"""Over-approximate time-to-react on a CommonRoad 2018b set-based scenario.

Follows Soentges, Koschi & Althoff (IV 2018): point-mass ego, |u_x|,|u_y| <= a_max,
velocity box, circular ego of radius r. The reachable set is over-approximated on a
position grid, each cell carrying a velocity box; a cell is removed only if every
position in it collides (so the result stays an over-approximation). TTR_max is the
first start time along the intended trajectory from which the set is empty at t_f.

Usage:
    python commonroad_ttr.py DEU_Ffb-2_2_S-1.xml --h 0.125 [--road] [--keep 2:0] [--from-step K]

Needs Python 3.11+, NumPy and Shapely (versions in environment/requirements-grid.txt).
With --from-step, the narrowed prediction applies only to start times from step K on, so the
empty/nonempty pattern need not be monotone. Every start time is then logged, and TTR_max is the
first start time with no evasive trajectory under the information available at that time.
"""
import sys, math, argparse
import xml.etree.ElementTree as ET
import numpy as np
from shapely.geometry import Polygon, Point, box
from shapely.ops import unary_union
from shapely.prepared import prep

DT, T_F = 0.1, 3.0
NSTEP = int(round(T_F / DT))


def pts(el):
    return [(float(p.find('x').text), float(p.find('y').text)) for p in el.iter('point')]


def load(path):
    r = ET.parse(path).getroot()
    lanelets = {}
    for ll in r.findall('lanelet'):
        L, R = pts(ll.find('leftBound')), pts(ll.find('rightBound'))
        succ = [s.get('ref') for s in ll.findall('successor')]
        lanelets[ll.get('id')] = dict(L=np.array(L), R=np.array(R), succ=succ)
    obs = {}
    for ob in r.findall('obstacle'):
        occ = {}
        for o in ob.find('occupancySet').iter('occupancy'):
            t = o.find('time')
            k = int(t.find('intervalStart').text) if t.find('intervalStart') is not None else int(t.find('exact').text)
            polys = [Polygon(pts(p)) for p in o.find('shape')]
            occ[k] = polys
        obs[ob.get('id')] = occ
    pp = r.find('planningProblem').find('initialState')
    g = lambda tag: float(pp.find(tag).find('exact').text)
    p0 = np.array(pts(pp.find('position'))[0])
    ego = dict(p=p0, v=g('velocity'), th=g('orientation'))
    return lanelets, obs, ego


def centerline(ll):
    n = min(len(ll['L']), len(ll['R']))
    return (ll['L'][:n] + ll['R'][:n]) / 2


def intended_path(lanelets, ego, length):
    # start lanelet: the one whose polygon contains the ego position
    start = None
    for i, ll in lanelets.items():
        poly = Polygon(np.vstack([ll['L'], ll['R'][::-1]]))
        if poly.buffer(1e-6).contains(Point(ego['p'])):
            start = i
            break
    assert start is not None
    path = [ego['p']]
    cl = centerline(lanelets[start])
    # keep the part of the centreline ahead of the ego
    d = cl - ego['p']
    ahead = d @ np.array([math.cos(ego['th']), math.sin(ego['th'])]) > 0
    path += list(cl[ahead])
    cur = start
    while True:
        P = np.array(path)
        if np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1)) > length or not lanelets[cur]['succ']:
            break
        h = P[-1] - P[-2]
        best, bestc = None, -2
        for s in lanelets[cur]['succ']:
            c2 = centerline(lanelets[s])
            dd = c2[-1] - c2[0]  # straight on: successor whose overall direction best matches the heading
            c = (h @ dd) / (np.linalg.norm(h) * np.linalg.norm(dd) + 1e-9)
            if c > bestc:
                best, bestc = s, c
        cur = best
        path += list(centerline(lanelets[cur])[1:])
    return np.array(path)


def along(path, s):
    path = path[np.r_[True, np.linalg.norm(np.diff(path, axis=0), axis=1) > 1e-12]]  # drop zero-length segments
    seg = np.diff(path, axis=0)
    L = np.linalg.norm(seg, axis=1)
    cum = np.concatenate([[0], np.cumsum(L)])
    i = min(np.searchsorted(cum, s, side='right') - 1, len(seg) - 1)
    u = seg[i] / L[i]
    return path[i] + u * (s - cum[i]), u


def intended_states(path, speed, n, amax, vmax):
    """States on the intended trajectory at steps 0..n. The point mass follows the centreline velocity with constant
    acceleration within each step, so each state is reachable from the previous one under the model (as Proposition 2
    of Soentges, Koschi & Althoff requires); positions stay within a few cm of the centreline."""
    vref = [along(path, speed * k * DT)[1] * speed for k in range(n + 1)]
    P, V = [along(path, 0.0)[0]], [vref[0]]
    for k in range(n):
        u = (vref[k + 1] - V[k]) / DT
        assert np.all(np.abs(u) <= amax) and np.all(np.abs(vref[k + 1]) <= vmax)
        P.append(P[k] + DT * (V[k] + vref[k + 1]) / 2)
        V.append(vref[k + 1])
    return np.array(P), np.array(V)


def obstacle_union(obs, k, r, keep=None):
    polys = []
    for oid, occ in obs.items():
        ps = occ.get(k, [])
        if keep and oid in keep:
            ps = [ps[j] for j in keep[oid] if j < len(ps)]
        polys += ps
    return prep(unary_union(polys).buffer(r)) if polys else None


def reach_empty(p0, v0, k0, args, colliders, road):
    """Propagate from exact state (p0, v0) at step k0 to the final step; True if empty."""
    h, a, vmax = args.h, args.amax, args.vmax
    # cells: arrays of (ix, iy, vx0, vx1, vy0, vy1)
    ix, iy = math.floor(p0[0] / h), math.floor(p0[1] / h)
    C = np.array([[ix, iy, v0[0], v0[0], v0[1], v0[1]]], float)
    for k in range(k0 + 1, NSTEP + 1):
        x0, y0 = C[:, 0] * h, C[:, 1] * h
        out = []
        for ax, s0, (lo, hi) in ((0, x0, (2, 3)), (1, y0, (4, 5))):
            s1 = s0 + h
            vl, vh = C[:, lo], C[:, hi]
            nlo = s0 + vl * DT - 0.5 * a * DT ** 2
            nhi = s1 + vh * DT + 0.5 * a * DT ** 2
            out.append((np.floor(nlo / h).astype(int), np.floor(nhi / h).astype(int), s0, s1, vl, vh))
        (ax0, ax1, xs0, xs1, vxl, vxh), (ay0, ay1, ys0, ys1, vyl, vyh) = out
        nx, ny = ax1 - ax0 + 1, ay1 - ay0 + 1
        cnt = nx * ny
        idx = np.repeat(np.arange(len(C)), cnt)
        off = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        tx = ax0[idx] + off % nx[idx]
        ty = ay0[idx] + off // nx[idx]

        def vrange(t, s0, s1, vl, vh):
            c0, c1 = t * h, t * h + h
            # s_new = s + v_new*DT - u*DT^2/2  =>  v_new in [(c0 - s1 - a DT^2/2)/DT, (c1 - s0 + a DT^2/2)/DT]
            lo = np.maximum.reduce([vl - a * DT, (c0 - s1 - 0.5 * a * DT ** 2) / DT, np.full_like(vl, -vmax)])
            hi = np.minimum.reduce([vh + a * DT, (c1 - s0 + 0.5 * a * DT ** 2) / DT, np.full_like(vh, vmax)])
            return lo, hi

        vx0n, vx1n = vrange(tx, xs0[idx], xs1[idx], vxl[idx], vxh[idx])
        vy0n, vy1n = vrange(ty, ys0[idx], ys1[idx], vyl[idx], vyh[idx])
        ok = (vx0n <= vx1n) & (vy0n <= vy1n)
        tx, ty, vx0n, vx1n, vy0n, vy1n = tx[ok], ty[ok], vx0n[ok], vx1n[ok], vy0n[ok], vy1n[ok]
        if len(tx) == 0:
            return True
        # merge per cell: hull of velocity boxes
        order = np.lexsort((ty, tx))
        tx, ty = tx[order], ty[order]
        vx0n, vx1n, vy0n, vy1n = vx0n[order], vx1n[order], vy0n[order], vy1n[order]
        new = np.ones(len(tx), bool)
        new[1:] = (tx[1:] != tx[:-1]) | (ty[1:] != ty[:-1])
        start = np.flatnonzero(new)
        C = np.column_stack([tx[start], ty[start],
                             np.minimum.reduceat(vx0n, start), np.maximum.reduceat(vx1n, start),
                             np.minimum.reduceat(vy0n, start), np.maximum.reduceat(vy1n, start)])
        # remove cells that collide everywhere (cell inside buffered obstacle) or lie wholly off-road
        coll = colliders[k]
        keep = np.ones(len(C), bool)
        for j, (cx, cy) in enumerate(C[:, :2]):
            b = box(cx * h, cy * h, cx * h + h, cy * h + h)
            if coll is not None and coll.contains(b):
                keep[j] = False
            elif road is not None and road.disjoint(b):
                keep[j] = False
        C = C[keep]
        if len(C) == 0:
            return True
    return False


def parse_keep(spec, obs):
    """Read --keep ("2:0" or "2:0,1;3:1") and check that each obstacle and polygon exists."""
    keep = {}
    for part in filter(None, spec.split(';')):
        try:
            oid, js = part.split(':')
            keep[oid] = [int(j) for j in js.split(',')]
        except ValueError:
            sys.exit(f'--keep: cannot read "{part}"; expected OBSTACLE:POLYGON[,POLYGON...], e.g. 2:0')
        if oid not in obs:
            sys.exit(f'--keep: obstacle {oid} is not in the scenario (obstacles: {", ".join(sorted(obs))})')
        n_poly = min(len(polys) for polys in obs[oid].values())
        if any(j < 0 or j >= n_poly for j in keep[oid]):
            sys.exit(f'--keep: obstacle {oid} has polygons 0 to {n_poly - 1} at every step')
    return keep


def main():
    ap = argparse.ArgumentParser(description='Grid over-approximation of the time-to-react (TTR_max).')
    ap.add_argument('scenario', help='CommonRoad 2018b scenario file (XML) with set-based predictions')
    ap.add_argument('--h', type=float, default=0.5, help='grid cell size (m)')
    ap.add_argument('--amax', type=float, default=10.0, help='acceleration bound per axis (m/s^2)')
    ap.add_argument('--vmax', type=float, default=14.0, help='speed bound per axis (m/s)')
    ap.add_argument('--r', type=float, default=0.9, help='radius of the disc modelling the ego vehicle (m)')
    ap.add_argument('--road', action='store_true', help='also remove cells wholly off the road')
    ap.add_argument('--keep', default='', help='narrowed prediction, e.g. "2:0" keeps polygon 0 of obstacle 2')
    ap.add_argument('--from-step', type=int, default=0,
                    help='narrowed prediction applies to start times >= this step (steps of 0.1 s)')
    a = ap.parse_args()
    lanelets, obs, ego = load(a.scenario)
    keep = parse_keep(a.keep, obs)
    road = None
    if a.road:
        road = prep(unary_union([Polygon(np.vstack([l['L'], l['R'][::-1]])) for l in lanelets.values()]))
    path = intended_path(lanelets, ego, ego['v'] * T_F + 5)
    full = {k: obstacle_union(obs, k - 1, a.r) for k in range(1, NSTEP + 1)}
    narrow = {k: obstacle_union(obs, k - 1, a.r, keep) for k in range(1, NSTEP + 1)} if keep else full
    # TTC along the intended trajectory (with the prediction in force at t0)
    P_int, V_int = intended_states(path, ego['v'], NSTEP, a.amax, a.vmax)
    dev = max(np.linalg.norm(P_int[k] - along(path, ego['v'] * k * DT)[0]) for k in range(NSTEP + 1))
    print(f'intended trajectory: constant acceleration per step, max deviation from centreline {dev:.3f} m')
    def state(k):
        return P_int[k], V_int[k]
    def ttc(cols):
        for k in range(1, NSTEP + 1):
            p, _ = state(k)
            if cols[k] is not None and cols[k].contains(Point(p)):
                return k
        return None
    print('TTC steps (full prediction):', ttc(full), ' (narrowed):', ttc(narrow) if keep else '-')
    res = []
    for k0 in range(0, NSTEP):
        cols = narrow if (keep and k0 >= a.from_step) else full
        p, v = state(k0)
        if cols.get(k0) is not None and k0 > 0 and cols[k0].contains(Point(p)):
            res.append((k0, 'collided'))
            break
        e = reach_empty(p, v, k0, a, cols, road)
        res.append((k0, 'empty' if e else 'nonempty'))
        print(f't*={k0*DT:.1f}s reachable set at t_f: {"EMPTY" if e else "nonempty"}', flush=True)
        if e and not (keep and a.from_step > 0):
            break  # fixed prediction: no evasive trajectory from any later start (Soentges et al., Prop. 2)
    first_empty = next((k for k, s in res if s == 'empty'), None)
    if first_empty is None and res and res[-1][1] == 'collided':
        print(f't*={res[-1][0]*DT:.1f}s ego on intended trajectory collides')
        print('TTR_max = None (no empty set found before the collision)')
    else:
        print('TTR_max =', None if first_empty is None else f'{first_empty*DT:.1f} s')


if __name__ == '__main__':
    main()
