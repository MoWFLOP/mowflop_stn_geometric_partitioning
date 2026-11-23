# A Novel Geometric Partitioning for Search Trajectory Networks in Multi-Objective Wind Farm Layout Optimization

This repository contains the implementation, benchmark instances, wind data and experimental results of the paper above, submitted to the IEEE Congress on Evolutionary Computation (CEC 2027). It includes the MOEA/D and NSGA-II implementations with the trajectory logger, the code that builds Search Trajectory Networks (STNs) with the Shannon-entropy scheme and with the proposed occupancy grid, and the scripts that reproduce the tables and figures of the paper.

## Repository structure

```
.
├── instance_map.csv              # Paper instance number -> internal name (Table I)
├── STN_MoWFLOP/                  # Metaheuristics and STN logger (C++17)
│   ├── source_code/              # MOEA/D, NSGA-II, STN logger, Makefile, run scripts
│   ├── instances/                # Wind roses, power curves, generator of instances 1-4
│   └── wflop_instances/          # The 14 instances used in the paper
├── STNs-MOCO-MoWFLOP/            # Partitioning, STN construction, metrics and plots
│   ├── scripts/                  # R scripts and Python package (mowflop)
│   ├── instances/site/           # Site geometry used by the occupancy grid
│   ├── raw_results/              # Trajectories, wind scenarios, runs of Silva et al.
│   └── sharing/                  # STNs, reference sets and metrics of the paper
├── paper/                        # Scripts for the tables and figures
└── README.md
```

## Requirements

- C++17 compiler (g++) and make.
- Python 3.10+ with the packages in `STNs-MOCO-MoWFLOP/requirements.txt` (pandas, matplotlib, scipy).
- R 4.5 with igraph, ggplot2 (>= 4.0), ggraph, tidygraph, dplyr, tidyr, RColorBrewer and reshape2. The figures of the paper were drawn with R 4.5.2, ggplot2 4.0.2, ggraph 2.2.2 and igraph 2.1.4.

Python modules under `STNs-MOCO-MoWFLOP/scripts/mowflop` are run from `STNs-MOCO-MoWFLOP/scripts`.

## Instances

| Paper | Internal | τ | \|P\| | Zones |
|---|---|---|---|---|
| 1 | ns506 (`506_e-02`) | 5 | 60 | 1 |
| 2 | ns507 (`507_e-03`) | 5 | 600 | 1 |
| 3 | ns508 (`508_e-04`) | 5 | 6,000 | 1 |
| 4 | ns509 (`509_e-05`) | 5 | 60,000 | 1 |
| 5 | ns48 | 21 | 4,229 | 1 |
| 6 | ns178 | 23 | 3,963 | 1 |
| 7 | ns192 | 31 | 4,113 | 1 |
| 8 | ns488 | 40 | 6,610 | 1 |
| 9 | ns203 | 41 | 6,734 | 2 |
| 10 | ns202 | 47 | 8,797 | 1 |
| 11 | ns101 | 63 | 4,727 | 1 |
| 12 | ns465 | 79 | 14,709 | 1 |
| 13 | ns41 | 123 | 13,248 | 3 |
| 14 | ns440 | 140 | 28,050 | 2 |

Instances 1-4 are in `STN_MoWFLOP/wflop_instances/sparse_instances/single_site_sweep/` and instances 5-14 in `STN_MoWFLOP/wflop_instances/New Sites/`. `instance_map.csv` holds the same table (Table I).

## Running the metaheuristics

```bash
cd STN_MoWFLOP/source_code/meta_heuristics && make && cd ..
# <instance> <algorithm> <run> <angle> <wind speed> <evaluations> <weight vectors p> <recording interval>
bash meta_heuristics/scripts/run_one_windcorrected.sh ns192 nsga2 0 330 13 1000000 100 50
```

The angle and wind speed of each run of the paper are in `STNs-MOCO-MoWFLOP/raw_results/wind_corrected/` (`cec_wind_map.csv` for instances 5-14, `sparse_wind_map_506-513_corrected.csv` for instances 1-4). The trajectory is written to `<instance>_<algorithm>_stn.csv` with one row per weight vector every 50 generations (`algorithm,instance,run_id,vector_id,generation,iteration,f_cost,f_power,weight1,weight2,occupied`).

## Reproducing the paper

| Item | Script |
|---|---|
| Table I | `instance_map.csv` |
| Tables III and IV, numbers in Section V | `python paper/make_tables.py` |
| Fig. 3 | `mowflop.schemes.shannon_entropy.diagnose_entropy` → `paper/entropy_csv_to_npy.py` → `paper/entropy_curve.py` |
| Fig. 4 | `paper/gerar_rq1_disc_inv.py` |
| Rank stability (Section V-A) | `mowflop.schemes.shannon_entropy.rank_stability` |
| Figs. 5–8 | `mowflop.run_scripts.run_plot_paper_r` on the STNs in `sharing/campanhas/*/stns/` (force-directed layouts are not deterministic: node positions change, the graphs do not) |

Python modules are run from `STNs-MOCO-MoWFLOP/scripts` (`python -m <module> --help`).

The STNs, reference sets and metrics can be rebuilt from `raw_results/` with `STNs-MOCO-MoWFLOP/scripts/run_entropy_campaign.sh`. The reference set of each wind scenario also uses the runs of Silva et al. with that scenario; for instances 1–4, which have none, the final archives of p = 10 runs 10–29 take their place (`raw_results/wflopcec26/<algo>/ns50x/`).

## Credits

The metaheuristics build on the code of [WFLOP_meta_learning](https://github.com/wflopbracis/WFLOP_meta_learning), which also provides the real-world instances and the MOEA/D and NSGA-II runs of Silva et al. The STN construction, metrics and plots build on the code of [STNs-MOCO](https://github.com/gabro8a/STNs-MOCO) by Ochoa, Liefooghe, Lavinas and Aranha.
