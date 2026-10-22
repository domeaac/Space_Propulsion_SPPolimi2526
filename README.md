# Space_Propulsion_SPPolimi2526

My codes for the Space Propulsion exam @Polimi MSc Space Engineering, AA2025-2026.

The repository contains the two course assignments, both written in MATLAB:

| Project | Folder | Report |
|---|---|---|
| **BARIA**: experimental characterisation of a BATES solid rocket motor (burning-rate law, c*, ballistic simulation) | [BARIA Solid Rocket Engine Analysis/](BARIA%20Solid%20Rocket%20Engine%20Analysis/) | [BARIA-SRM-Report.pdf](BARIA-SRM-Report.pdf) |
| **Stark100k**: design of a 100 kN-class solid rocket engine (grain, nozzle, structure, Monte Carlo, thermal protection, hydraulics) | [Stak100k Rocket Engine Development/](Stak100k%20Rocket%20Engine%20Development/) | [Stark100kReport.pdf](Stark100kReport.pdf) |

The assignment brief for Stark100k is in [Stark100k Assignment-V03.pdf](Stark100k%20Assignment-V03.pdf).

## Repository layout

```
Space_Propulsion_SPPolimi2526/
├── README.md                              this file
├── BARIA-SRM-Report.pdf                   BARIA final report
├── Stark100k Assignment-V03.pdf           Stark100k assignment brief
├── Stark100kReport.pdf                    Stark100k final report
├── BARIA Solid Rocket Engine Analysis/    BARIA code + experimental data (see its README.md)
└── Stak100k Rocket Engine Development/    Stark100k code, data and libraries (see its README.txt)
```

## Requirements

- MATLAB (developed on a recent release)
- Statistics and Machine Learning Toolbox: BARIA summary plots (`ksdensity`, `normpdf`, `quantile`)
- Optimization Toolbox: only for the optional `solverCoolantChannelFminCon.m` variant in Stark100k
- NASA CEA for MATLAB: Stark100k thermochemistry. **For security reasons the CEA files are not provided in this repository.** Get your own copy of `CEA.m`, `CEA.p` and `thermo_lib.mat` and put them in all of these folders under `Stak100k Rocket Engine Development/`:
  - `src/CEA/`: `CEA.m`, `CEA.p`, `thermo_lib.mat`
  - `src/complete_cooling/+CEAThermo/`: `CEA.m`, `CEA.p`, `thermo_lib.mat`
  - `src/complete_cooling/data/`: `thermo_lib.mat`

## Quick start

- **BARIA**: open MATLAB in `BARIA Solid Rocket Engine Analysis/` and run `main3.m`.
- **Stark100k**: open MATLAB in `Stak100k Rocket Engine Development/` and run `main.m` (whole or by section). Do not run `engine_montecarlo.m` unless you need to regenerate the Monte Carlo results; it takes a long time. `engine_resultsMontecarlo.m` loads the saved results instead.

Each project folder has its own README with details.
