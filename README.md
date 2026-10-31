# Additional Information and Time to Use It: code and data

Code, logs and results for the paper *Additional Information and Time to Use It: Two-Boundary Recoverability* (L. Razo, 2026). The paper shows
that, as long as the task stays the same, a second source of information can improve a decision only if it is new
and arrives before the last moment a correction can still be started. This repository reproduces every number,
table and figure the paper computes.

The paper is under review at *Theory and Decision*. A copy is available from the author on request.

## Contents

| Folder | What it holds |
|---|---|
| `deadline_model/` | The deadline model of Section 3: the value of checking, Figure 1 and the numerical checks of the proposition |
| `commonroad/` | The intersection example of Appendix A: two time-to-react computations on a CommonRoad scenario, with every log; details in `commonroad/README.md` |
| `environment/` | Pinned package versions for each part |

## Claim-to-file map

| Paper item | Claim | Produced by | Evidence |
|---|---|---|---|
| Figure 1 | Region where an independent check pays, fixed deadline | `deadline_model/recoverability_deadline_model.py` | `deadline_model/results/fig1_where_checking_pays.png` |
| Section 3, "Where checking pays" | A sweep over 2,187 combinations of the assumed values, with a fixed deadline and three moving ones, agrees with the condition, with a search over every decision rule and with each part of the proposition | same | `deadline_model/results/deadline_model_results.json`, `sweep` |
| same | 3,000 random settings agree; in 0 of 3,000 further settings does scaling both chances of a move up in proportion lower the value | same | same file, `random_settings` |
| Section 3, "Postponing the step" | Postponing never worsens the best decision, and the best delay is at most e₀L(0)/h | same | same file, `postponement` |
| Table 2, CommonRoad-Reach row | Time-to-react 0.4, 0.8 and 0.4 s (turn unknown, same way, across) at grid 0.2 m and split radius 0.2 and 0.1 m | `commonroad/code/commonroad_reach_crosscheck.py`, run by the notebook | `commonroad/logs/reach_runs/` |
| Table 2 caption | The default split radius of 0.7 m gives 0.5 s for the same-way case | same | `commonroad/logs/reach_runs/reach_keep20_grid0{1,2}_split07.log` |
| Table 2, grid row | 0.5, 0.9 and 0.5 s at 0.125 m, with and without the road constraint | `commonroad/code/commonroad_ttr.py` | `commonroad/logs/grid_runs/` |
| Appendix A, "Computed deadlines" | The same-way gain is kept if the turn is known by 0.4 s and lost at 0.5 s (CommonRoad-Reach); the same pattern 0.1 s later in the grid computation | both tools, `--from-step` | `commonroad/logs/later_intent/` |

The grid computation is an upper bound for the stated model. CommonRoad-Reach can discard small regions that are only
partly blocked, so its values are close estimates rather than certified bounds. The scenario `DEU_Ffb-2_2_S-1` is
taken to be the T-intersection studied by Söntges, Koschi and Althoff (2018).

## Reproducing

From the repository root:

```
pip install -r environment/requirements-deadline.txt
python deadline_model/recoverability_deadline_model.py
```

This writes the figure and the results file to `deadline_model/output/`, for comparison with `deadline_model/results/`,
and exits with an error if any check fails. It takes about a minute.

```
pip install -r environment/requirements-grid.txt
python commonroad/code/commonroad_ttr.py commonroad/scenario/DEU_Ffb-2_2_S-1.xml --h 0.125 [--road] [--keep 2:0] [--from-step K]
```

The grid computation runs on any platform; a run at 0.125 m takes a few minutes to a quarter of an hour.
CommonRoad-Reach needs Linux and Python 3.10 or 3.11: open `commonroad/code/commonroad_reach_crosscheck.ipynb` in
Colab (CPU runtime) and run all cells. Every log begins with the command that produced it. Python 3.11 or later is
assumed throughout.

## Licences

- Code: MIT (`LICENSE`).
- Logs, results and figures: CC BY 4.0 (`LICENSE-data`).
- The scenario `commonroad/scenario/DEU_Ffb-2_2_S-1.xml` is third-party material from the CommonRoad scenario
  repository (`gitlab.lrz.de/tum-cps/commonroad-scenarios`), authors Vanessa Bui, Markus Koschi and Sebastian Lutz,
  Technical University of Munich, distributed unmodified under its BSD 3-Clause licence
  (`commonroad/scenario/LICENSE-scenario.txt`). Its inclusion does not imply endorsement by the Technical University
  of Munich.

## References

- M. Althoff, M. Koschi and S. Manzinger, "CommonRoad: Composable benchmarks for motion planning on roads", IEEE
  Intelligent Vehicles Symposium, 2017.
- S. Söntges, M. Koschi and M. Althoff, "Worst-case Analysis of the Time-To-React Using Reachable Sets", IEEE
  Intelligent Vehicles Symposium, 2018.
- E. I. Liu, G. Würsching, M. Klischat and M. Althoff, "CommonRoad-Reach: A Toolbox for Reachability Analysis of
  Automated Vehicles", IEEE International Conference on Intelligent Transportation Systems, 2022.

## Funding

This work was co-funded by the European Union under the Horizon Europe project EuropeanCity² (grant agreement No
101178170). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those
of the European Union or the European Research Executive Agency (REA). Neither the European Union nor the granting
authority can be held responsible for them.

## Citation

See `CITATION.cff`.
