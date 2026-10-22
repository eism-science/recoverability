"""Time-to-react (TTR) with CommonRoad-Reach for the scenario DEU_Ffb-2_2_S-1.

Parses the CommonRoad 2018b scenario file and builds the scenario with commonroad-io, keeping the set-based
occupancy predictions unchanged. TTR_max is the first start time along the ego vehicle's intended trajectory from
which its reachable set at t_f = 3.0 s is empty (Soentges, Koschi and Althoff, IEEE IV 2018).

Model: point mass in Cartesian coordinates, |a_x|, |a_y| <= 10 m/s^2, |v_x|, |v_y| <= 14 m/s, ego vehicle as a disc
of radius 0.9 m, road boundaries enforced; at step k the occupancy predicted for the interval [k-1, k] is used. The
same model is implemented independently, as a grid over-approximation, in commonroad_ttr.py.

Requires Linux and Python 3.10 or 3.11 (CommonRoad-Reach 2025.2.1 with commonroad-io 2024.3); see
commonroad_reach_crosscheck.ipynb.

With --from-step K, the narrowed prediction applies only to start times at or after step K, and the full
prediction before it: the other car's intended direction becomes known at K * 0.1 s. The empty/nonempty pattern
then need not be monotone, so every start time up to a collision is logged, and TTR_max is the first start time
with no evasive trajectory under the information available at that time.

CommonRoad-Reach removes a colliding rectangle whole once its diagonal is below --split, so its values are close
estimates, not certified upper bounds; commonroad_ttr.py gives the certified over-approximation.

Usage: python commonroad_reach_crosscheck.py SCENARIO.xml [--keep 2:0] [--from-step K] [--mode 2] [--grid 0.2]
       [--split 0.7] [--uncertainty U] [--build-only]
"""
import os
os.environ['MPLBACKEND'] = 'Agg'  # non-interactive plotting backend
import argparse
import math
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np

from commonroad.common.util import Interval
from commonroad.geometry.shape import Polygon, Rectangle, ShapeGroup
from commonroad.planning.goal import GoalRegion
from commonroad.planning.planning_problem import PlanningProblem
from commonroad.prediction.prediction import Occupancy, SetBasedPrediction
from commonroad.scenario.lanelet import Lanelet, LaneletNetwork
from commonroad.scenario.obstacle import DynamicObstacle, ObstacleType
from commonroad.scenario.scenario import Scenario, ScenarioID
from commonroad.scenario.state import CustomState, InitialState

DT, T_F = 0.1, 3.0
NSTEP = int(round(T_F / DT))
A_MAX, V_MAX, R_EGO = 10.0, 14.0, 0.9
OBS_ID_OFFSET = 900000  # obstacle ids must not collide with lanelet ids
UNCERTAINTY = ('uncertainty_p_lon', 'uncertainty_p_lat', 'uncertainty_v_lon', 'uncertainty_v_lat')


def pts(el):
    return [(float(p.find('x').text), float(p.find('y').text)) for p in el.iter('point')]


def load_2018b(path):
    """Parse the 2018b file: lanelets, set-based occupancies (step -> list of polygons) and the ego initial state."""
    r = ET.parse(path).getroot()
    lanelets = {}
    for ll in r.findall('lanelet'):
        L, R = np.array(pts(ll.find('leftBound'))), np.array(pts(ll.find('rightBound')))
        assert len(L) == len(R), f"lanelet {ll.get('id')}: bound lengths differ"
        lanelets[int(ll.get('id'))] = dict(L=L, R=R, succ=[int(s.get('ref')) for s in ll.findall('successor')])
    obs = {}
    for ob in r.findall('obstacle'):
        occ = {}
        for o in ob.find('occupancySet').iter('occupancy'):
            t = o.find('time')
            k = int(t.find('intervalStart').text)
            assert int(t.find('intervalEnd').text) == k + 1
            occ[k] = [np.array(pts(p)) for p in o.find('shape')]
        obs[int(ob.get('id'))] = occ
    init = r.find('planningProblem').find('initialState')
    g = lambda tag: float(init.find(tag).find('exact').text)
    ego = dict(p=np.array(pts(init.find('position'))[0]), v=g('velocity'), th=g('orientation'))
    return lanelets, obs, ego


# --- intended trajectory (same construction as commonroad_ttr.py) ---------------------------------------------
def centerline(ll):
    return (ll['L'] + ll['R']) / 2


def intended_path(lanelets, ego, length):
    from shapely.geometry import Point, Polygon as SPolygon
    start = None
    for i, ll in lanelets.items():
        if SPolygon(np.vstack([ll['L'], ll['R'][::-1]])).buffer(1e-6).contains(Point(ego['p'])):
            start = i
            break
    assert start is not None
    path = [ego['p']]
    cl = centerline(lanelets[start])
    ahead = (cl - ego['p']) @ np.array([math.cos(ego['th']), math.sin(ego['th'])]) > 0
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
            dd = c2[-1] - c2[0]
            c = (h @ dd) / (np.linalg.norm(h) * np.linalg.norm(dd) + 1e-9)
            if c > bestc:
                best, bestc = s, c
        cur = best
        path += list(centerline(lanelets[cur])[1:])
    return np.array(path)


def along(path, s):
    seg = np.diff(path, axis=0)
    L = np.linalg.norm(seg, axis=1)
    cum = np.concatenate([[0], np.cumsum(L)])
    i = min(np.searchsorted(cum, s, side='right') - 1, len(seg) - 1)
    u = seg[i] / L[i]
    return path[i] + u * (s - cum[i]), u


def intended_states(path, speed, n, amax=10.0, vmax=14.0):
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


# --- commonroad-io scenario ---------------------------------------------------------------------------------------
def build_scenario(lanelets, obs, keep):
    sc = Scenario(dt=DT, scenario_id=ScenarioID(country_id="DEU", map_name="Ffb", map_id=2))
    net = LaneletNetwork()
    for i, ll in lanelets.items():
        net.add_lanelet(Lanelet(left_vertices=ll['L'], center_vertices=centerline(ll), right_vertices=ll['R'],
                                lanelet_id=i, successor=ll['succ'] or None))
    sc.add_objects(net)
    for oid, occ in obs.items():
        occset = []
        for k in range(1, NSTEP + 1):  # step k uses the occupancy of interval [k-1, k], as in commonroad_ttr.py
            polys = occ.get(k - 1, [])
            if keep and oid in keep:
                polys = [polys[j] for j in keep[oid] if j < len(polys)]
            if polys:
                occset.append(Occupancy(k, ShapeGroup([Polygon(p) for p in polys])))
        c0 = np.mean(occ[0][0], axis=0)
        init = InitialState(position=c0, orientation=0.0, velocity=0.0, acceleration=0.0, yaw_rate=0.0,
                            slip_angle=0.0, time_step=0)
        sc.add_objects(DynamicObstacle(OBS_ID_OFFSET + oid, ObstacleType.CAR, Rectangle(0.1, 0.1), init,
                                       SetBasedPrediction(1, occset)))
    return sc


def planning_problem(ego):
    init = InitialState(position=ego['p'], orientation=ego['th'], velocity=ego['v'], acceleration=0.0,
                        yaw_rate=0.0, slip_angle=0.0, time_step=0)
    goal = GoalRegion([CustomState(time_step=Interval(0, NSTEP))])
    return PlanningProblem(1, init, goal)


# --- CommonRoad-Reach -------------------------------------------------------------------------------------------
def reach_empty(scenario, pp, k0, p, v, th, args):
    """True if Reach's reachable set at step NSTEP, started from the exact state (p, v) at step k0, is empty."""
    from commonroad_reach.data_structure.configuration_builder import ConfigurationBuilder
    from commonroad_reach.data_structure.reach.reach_interface import ReachableSetInterface

    config = ConfigurationBuilder(path_root=args.workdir).build_configuration("DEU_Ffb-2_2_S-1")
    config.planning.coordinate_system = "CART"
    config.planning.steps_computation = NSTEP - k0
    config.vehicle.ego.v_max = V_MAX
    config.vehicle.ego.a_max = A_MAX
    config.reachable_set.mode_computation = args.mode
    config.reachable_set.size_grid = args.grid
    config.reachable_set.radius_terminal_split = args.split
    config.reachable_set.consider_traffic = True
    config.reachable_set.prune_nodes_not_reaching_final_step = False
    if hasattr(config, 'debug'):
        config.debug.save_config = False
    config.update(scenario=scenario, planning_problem=pp)
    # override: start from the state on the intended trajectory at step k0, exact (tiny uncertainty as in defaults)
    config.planning.set_initial_states(k0, (float(p[0]), float(p[1])), (float(v[0]), float(v[1])), float(th))
    config.vehicle.ego.radius_inflation = R_EGO
    missing = [n for n in UNCERTAINTY if not hasattr(config.planning, n)]
    assert not missing, f'Reach configuration has no {missing}'
    default = ', '.join(f'{n}={getattr(config.planning, n)}' for n in UNCERTAINTY)
    if args.uncertainty is not None:  # initial position (m) and velocity (m/s) uncertainty, replacing Reach's defaults
        for name in UNCERTAINTY:
            setattr(config.planning, name, args.uncertainty)
    if k0 == 0:
        print('initial uncertainty: default', default, '| used',
              ', '.join(f'{n}={getattr(config.planning, n)}' for n in UNCERTAINTY), flush=True)
    rs = ReachableSetInterface(config)
    rs.compute_reachable_sets(verbose=False)
    final = rs.reachable_set_at_step(NSTEP)
    return final is None or len(final) == 0


def parse_keep(spec, obs):
    """Read --keep ("2:0" or "2:0,1;3:1") and check that each obstacle and polygon exists."""
    keep = {}
    for part in filter(None, spec.split(';')):
        try:
            oid, js = part.split(':')
            oid = int(oid)
            keep[oid] = [int(j) for j in js.split(',')]
        except ValueError:
            sys.exit(f'--keep: cannot read "{part}"; expected OBSTACLE:POLYGON[,POLYGON...], e.g. 2:0')
        if oid not in obs:
            sys.exit(f'--keep: obstacle {oid} is not in the scenario (obstacles: {", ".join(map(str, sorted(obs)))})')
        n_poly = min(len(polys) for polys in obs[oid].values())
        if any(j < 0 or j >= n_poly for j in keep[oid]):
            sys.exit(f'--keep: obstacle {oid} has polygons 0 to {n_poly - 1} at every step')
    return keep


def main():
    ap = argparse.ArgumentParser(description='Time-to-react (TTR_max) with CommonRoad-Reach.')
    ap.add_argument('scenario', help='CommonRoad 2018b scenario file (XML) with set-based predictions')
    ap.add_argument('--keep', default='', help='narrowed prediction, e.g. "2:0" keeps polygon 0 of obstacle 2')
    ap.add_argument('--from-step', type=int, default=0,
                    help='the narrowed prediction applies to start times at or after this step (steps of 0.1 s)')
    ap.add_argument('--mode', type=int, default=2, help='Reach mode_computation: 1 Python, 2 C++')
    ap.add_argument('--grid', type=float, default=0.2, help="Reach size_grid, drivable-area discretisation (m)")
    ap.add_argument('--split', type=float, default=0.7,
                    help='Reach radius_terminal_split (m): a colliding rectangle whose diagonal is below this is '
                         'dropped whole')
    ap.add_argument('--uncertainty', type=float, default=None,
                    help="initial position (m) and velocity (m/s) uncertainty; default: Reach's own values")
    ap.add_argument('--workdir', default='.',
                    help="root folder for Reach's ConfigurationBuilder, where it looks for a configurations/ folder; "
                         "the logged runs used a folder without one")
    ap.add_argument('--build-only', action='store_true', help='build and check the scenario, skip Reach')
    a = ap.parse_args()
    print('CMD:', ' '.join(sys.argv), flush=True)

    lanelets, obs, ego = load_2018b(a.scenario)
    keep = parse_keep(a.keep, obs)
    sc_full = build_scenario(lanelets, obs, {})
    sc = build_scenario(lanelets, obs, keep) if keep else sc_full
    pp = planning_problem(ego)
    path = intended_path(lanelets, ego, ego['v'] * T_F + 5)

    NS = int(round(T_F / DT))
    P_int, V_int = intended_states(path, ego['v'], NS)
    dev = max(np.linalg.norm(P_int[k] - along(path, ego['v'] * k * DT)[0]) for k in range(NS + 1))
    print(f'intended trajectory: constant acceleration per step, max deviation from centreline {dev:.3f} m')
    def state(k):
        p, v = P_int[k], V_int[k]
        return p, v, math.atan2(v[1], v[0])

    # TTC along the intended trajectory with the (possibly narrowed) prediction, centre point against obstacle + r
    from shapely.geometry import Point, Polygon as SPolygon
    from shapely.ops import unary_union
    def collides(k, p, narrowed=True):
        polys = []
        for oid, occ in obs.items():
            ps = occ.get(k - 1, [])
            if narrowed and keep and oid in keep:
                ps = [ps[j] for j in keep[oid] if j < len(ps)]
            polys += [SPolygon(x) for x in ps]
        return bool(polys) and unary_union(polys).buffer(R_EGO).contains(Point(p))
    ttc = next((k for k in range(1, NSTEP + 1) if collides(k, state(k)[0])), None)
    print('scenario built:', len(sc.lanelet_network.lanelets), 'lanelets,', len(sc.dynamic_obstacles),
          'obstacles; TTC steps:', ttc, flush=True)
    if a.build_only:
        for o in sc.dynamic_obstacles:
            print('  obstacle', o.obstacle_id, 'occupancy at step 10:', o.occupancy_at_time(10) is not None)
        return

    scan_all = bool(keep) and a.from_step > 0  # prediction changes partway: log every start time
    first_empty = None
    for k0 in range(0, NSTEP):
        narrowed = bool(keep) and k0 >= a.from_step
        p, v, th = state(k0)
        if k0 > 0 and collides(k0, p, narrowed):
            print(f't*={k0*DT:.1f}s ego on intended trajectory collides')
            break
        t = time.time()
        e = reach_empty(sc if narrowed else sc_full, pp, k0, p, v, th, a)
        print(f't*={k0*DT:.1f}s reachable set at t_f: {"EMPTY" if e else "nonempty"}  ({time.time()-t:.1f} s)',
              flush=True)
        if e and first_empty is None:
            first_empty = k0
        if e and not scan_all:
            break  # fixed prediction: no evasive trajectory from any later start (Soentges et al., Prop. 2)
    if first_empty is not None:
        print(f'TTR_max = {first_empty*DT:.1f} s')
    else:
        print('TTR_max = None (no empty set found before the collision or the end of the horizon)')


if __name__ == '__main__':
    main()
