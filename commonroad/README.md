# CommonRoad time-to-react: code, data and results

Code, data and results for the intersection example of the paper *Additional Information and Time to Use It* (Appendix B).

The logs were produced in the authors' development layout, so the command line at the top of each log names
the older paths `notebooks/recoverability/` (scripts) and `data/recoverability/commonroad/` (scenario); the files
are otherwise unchanged.

## Question

At a T-intersection, the ego vehicle does not know which way the car on the side road will turn. The time-to-react
(TTR) is the latest start time along the ego's intended trajectory from which an evasive trajectory still exists,
that is, from which the ego's reachable set at the end of a 3 s horizon is non-empty (Söntges, Koschi and Althoff,
"Worst-case Analysis of the Time-To-React Using Reachable Sets", IEEE IV 2018). Both tools report, as that paper's
Algorithm 1 does, TTR_max: the first start time, in steps of 0.1 s, at which the reachable set at the end of the
horizon is empty. For the grid over-approximation (`commonroad_ttr.py`) this is an upper bound for the model stated
below. CommonRoad-Reach can discard small regions that are only partly blocked (see Model), so its values are close
estimates rather than certified upper bounds. The question is how the TTR changes when the other car's intent is known.

## Files

| Path | Content |
|---|---|
| `scenario/DEU_Ffb-2_2_S-1.xml` | Scenario, downloaded 2026-10-10 from `gitlab.lrz.de/tum-cps/commonroad-scenarios` (branch `master`, `scenarios/hand-crafted/`). Format 2018b, authors Bui, Koschi, Lutz, set-based predictions included |
| `scenario/scene.png` | Plot of the scenario, used to read the two branches of obstacle 2 |
| `code/commonroad_ttr.py` | Our implementation: grid over-approximation of the reachable set |
| `code/commonroad_reach_crosscheck.py` | The same computation with CommonRoad-Reach 2025.2.1 (Linux, Python 3.10–3.11) |
| `code/commonroad_reach_crosscheck.ipynb` | Colab notebook that installs CommonRoad-Reach and runs the cross-check |
| `logs/grid_early/` | Runs of `commonroad_ttr.py` at 0.5 and 0.25 m grids without the road constraint (`full_h05`, `full_h025` and the `keep` cases; the 0.5 m logs carry no suffix), command in each log |
| `logs/grid_runs/` | Our runs at 0.25 and 0.125 m with and without the road constraint, command in each log |
| `logs/later_intent/` | Same-way case with the intent known only from a later time: `commonroad_ttr.py` at 0.125 m (`keep20_h0125_from<K>.log`) and CommonRoad-Reach at grid 0.2 m, split radius 0.2 m (`reach_keep20_from<K>_grid02_split02.log`), command in each log |
| `logs/reach_checks/` | CommonRoad-Reach at a finer grid with a smaller split radius (`grid010_split005`) and with initial uncertainty 10⁻⁶ instead of the default 0.01 (`_exact`), settings in hundredths of a metre, and `versions.txt` |
| `logs/reach_runs/` | CommonRoad-Reach runs, named `reach_<case>_grid<g>_split<s>.log` (grid and split radius in tenths of a metre), and `versions.txt` |

## Scenario and cases

`DEU_Ffb-2_2_S-1` is taken to be the scenario `GER_Ffb_2b` of Söntges, Koschi and Althoff (2018): ego speed
14 m/s driving straight, three other vehicles, 3 s of set-based predictions. Obstacle 2, the car on the side road,
has two occupancy polygons per step, in the same order at every step (checked from the centroids): polygon 0 moves
west, from (14.6, 10.6) to (6.0, 5.0) at the last interval, the ego's direction; polygon 1 moves east and south, to (19.8, 1.2), across
the ego's path into the oncoming lane.

| Case | Prediction of obstacle 2 | Option |
|---|---|---|
| `full` | both branches (intent unknown) | none |
| `keep20` | branch 0 only (turns the same way as the ego) | `--keep 2:0` |
| `keep21` | branch 1 only (turns across) | `--keep 2:1` |

## Model

Shared by both tools (Söntges, Koschi and Althoff 2018, Table I, Scenario III): point mass in Cartesian coordinates,
|a_x|, |a_y| ≤ 10 m/s² and |v_x|, |v_y| ≤ 14 m/s, each axis bounded separately (so diagonal speeds up to about
19.8 m/s are allowed); ego as a disc of radius 0.9 m; Δt = 0.1 s; horizon 3 s. At step k the ego is checked against
the occupancy of the interval [k−1, k]; positions between steps are not checked.

- `commonroad_ttr.py` keeps a position grid with a velocity box per cell and removes a cell only if it collides
  everywhere with the obstacle model above (at step k, the occupancy of the interval [k−1, k], buffered by the ego
  radius). It therefore over-approximates the reachable set of this discretised model, and its TTR is an upper bound
  for that model, the stated point-mass dynamics, the reconstructed intended trajectory, the road geometry and the
  3 s horizon. It is not a certified bound for the obstacle's unknown position at every instant. `--road` also
  removes cells wholly off the road.
- The intended trajectory is a trajectory of the model itself, so every state on it is reachable from the previous
  one, which Söntges, Koschi and Althoff's Proposition 2 (a later start never leaves more evasive options) needs. Both
  scripts build it the same way (`intended_states`): at each 0.1 s step the point mass takes the velocity of the lane
  centreline at 14 m/s, with constant acceleration within the step (at most 5.4 m/s² per axis, against the 10 m/s²
  bound), and its position follows from the dynamics. It stays within 2 cm of the centreline over the first second
  and within 6 cm over the 3 s horizon; each log prints the deviation. Sampling positions directly on the piecewise-
  linear centreline, as earlier versions did, cut its corners slightly too sharply (by up to 9 mm at one step) for
  every sampled state to be reachable from the one before. The scenario file supplies the starting position and speed;
  the starting heading is the lane's (2.930 rad), 0.6° from the file's orientation of 2.94 rad, as in every version
  of both scripts. The computation therefore does not reproduce the file's planning problem exactly.
- CommonRoad-Reach (Cartesian mode, C++ backend) propagates polytopes and removes any colliding rectangle whose
  diagonal is below `radius_terminal_split` (default 0.7 m). The default trades precision for speed; smaller values
  resolve more of the free space next to obstacles. Road boundaries are always enforced.

## Results

TTR in seconds.

| Run | `full` | `keep20` | `keep21` |
|---|---|---|---|
| `commonroad_ttr.py`, 0.5 m grid | 0.6 | none* | 0.6 |
| `commonroad_ttr.py`, 0.25 m grid | 0.5 | 1.0 | 0.5 |
| `commonroad_ttr.py`, 0.125 m grid | 0.5 | 0.9 | 0.5 |
| Reach, grid 0.2 m, split 0.7 m (defaults) | 0.4 | 0.5 | 0.4 |
| Reach, grid 0.1 m, split 0.7 m | 0.4 | 0.5 | 0.4 |
| Reach, grid 0.2 m, split 0.2 m | 0.4 | 0.8 | 0.4 |
| Reach, grid 0.2 m, split 0.1 m | 0.4 | 0.8 | 0.4 |
| Söntges, Koschi and Althoff 2018 (full prediction, `GER_Ffb_2b`) | 0.5 upper bound, 0.3 achieved | – | – |

\* The set never empties before the ego collides on its intended trajectory at 1.2 s (12 steps), so the grid is too
coarse to give a TTR.

The road constraint changes no result at 0.25 or 0.125 m. Reach is unchanged between split radii 0.2 and 0.1 m, and
also at grid 0.1 m with split radius 0.05 m (`full` 0.4 s, `keep20` 0.8 s) and at grid 0.2 m, split radius 0.2 m with
the initial position and velocity uncertainty set to 10⁻⁶ instead of Reach's default 0.01 (`full` 0.4 s, `keep20`
0.8 s); logs in `logs/reach_checks/`.

Reading: knowing that the side-road car turns the same way doubles the time left to act, from 0.4 to 0.8 s in
CommonRoad-Reach (0.5 to 0.9 s in the grid over-approximation, which sits 0.1 s later in every case); knowing that it
turns across adds nothing, because that branch is the binding one. The default split radius of 0.7 m trades precision
for speed. At that setting the same-way case reads 0.5 s, and a radius of 0.2 m or less is needed to resolve the free
space beside the side-road car.

### Intent known later

The same-way prediction applies only from the time the intent becomes known; earlier start times use both branches
(option `--keep 2:0 --from-step K` in both tools; `commonroad_ttr.py` at 0.125 m with the road constraint,
CommonRoad-Reach at grid 0.2 m and split radius 0.2 m). Because the prediction changes partway through, a later start
time can again have an evasive trajectory once the intent is known, so with `--from-step` both tools log every start
time up to a collision or the end of the horizon, and report as TTR_max the first start time with no evasive
trajectory under the information then available. The rule behind that choice is stated under Assumption.

| Intent known at | 0.2 s | 0.3 s | 0.4 s | 0.5 s | 0.6 s |
|---|---|---|---|---|---|
| CommonRoad-Reach | 0.8 s | 0.8 s | 0.8 s | 0.4 s | 0.4 s |
| `commonroad_ttr.py` | 0.9 s | – | 0.9 s | 0.9 s | 0.5 s |

In both tools the gain is kept in full as long as the intent is known by the time-to-react under the full prediction
(0.4 s in CommonRoad-Reach, 0.5 s in the grid over-approximation) and the TTR is no longer extended once it is known
later: information that arrives after the last moment to act does not move that moment. Each evasive trajectory is
computed with the prediction that applies at its start time, so a manoeuvre begun before the intent is known does not
use that information later; what the information might still add to a manoeuvre already under way is not measured.
Söntges, Koschi and Althoff's result that a later start never leaves more evasive options holds for a fixed
prediction, and here the prediction changes partway through, so these values are the first start times from which no
evasive trajectory is found under the information then available, not upper bounds in their sense. The logs show the
change directly. With the intent known at 0.5 s, CommonRoad-Reach finds no evasive trajectory from 0.4 s, under the
full prediction, and finds one again from 0.5 to 0.7 s, under the narrowed one; with it known at 0.6 s, the set is
empty at 0.5 s and non-empty from 0.6 to 0.7 s. The grid over-approximation shows the same pattern 0.1 s later.

## Reproducing

Grid over-approximation (any platform; NumPy and Shapely), from the repository root:

```
python commonroad/code/commonroad_ttr.py commonroad/scenario/DEU_Ffb-2_2_S-1.xml --h 0.125 [--road] [--keep 2:0]
```

CommonRoad-Reach: open `code/commonroad_reach_crosscheck.ipynb` in Colab (CPU runtime) and
run all cells. The notebook downloads the script and the scenario, and asks for them if a download fails.
`commonroad-io` is pinned to 2024.3, because 2026.1 is incompatible with CommonRoad-Reach 2025.2.1.

## Assumption

The later-intent results assume that the ego vehicle must keep an evasive trajectory available at every moment. Under
that rule, once no evasive trajectory exists under the information the vehicle has, it must already have started one,
so information that would restore an escape at a later start time comes too late.

All logs were produced on 2026-10-22 with the current scripts: the 18 `commonroad_ttr.py` logs on the machine whose
versions are in `environment/requirements-grid.txt`, and the 21 CommonRoad-Reach logs in Colab, with the full package list in
`environment/versions_reach.txt`. Earlier changes (an exact cell sort, and an intended trajectory that is a trajectory of the model) and
the changes of 2026-10-22 (checks on `--keep`, help text, and logging every start time with `--from-step`) were each
followed by a rerun of all logs, and none changed any value. The later-intent logs now continue past the first empty
set; every line they share with the earlier runs is unchanged.
