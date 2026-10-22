# BARIA Solid Rocket Engine Analysis

Space Propulsion assignment, Politecnico di Milano, AA2025-2026.

This folder processes the static-firing pressure traces of the **BARIA** BATES motor. It:

1. extracts the action and burning time intervals from each trace,
2. fits the Vieille burning-rate law `r_b = a·p^n` with uncertainties,
3. estimates the experimental characteristic velocity `c*` and its spread,
4. simulates the mid-pressure motor's ballistics with a hand-written RK4 integrator.

The full discussion of the results is in [../BARIA-SRM-Report.pdf](../BARIA-SRM-Report.pdf).

## How to run

1. Set MATLAB's working directory to this folder.
2. Run `main3.m`, either all at once or section by section.

The helper files are scripts, not functions, so they share `main3.m`'s workspace. Run them only through `main3.m` and in its order.

## Motor and propellant data

These values are set at the top of `main3.m`.

| Quantity | Value |
|---|---|
| Grain | BATES (tubular, inhibited outer surface, both ends burning) |
| Inner radius `Ri0` | 50 mm |
| Web thickness `wt` | 30 mm (`Re0` = 80 mm) |
| Length `L0` | 290 mm |
| Throat diameters (LP / MP / HP) | 28.80 / 25.26 / 21.81 mm |
| Composition (mass fraction) | AP 68 % · Al 18 % · HTPB binder 14 % |
| Component densities | 1950 / 2700 / 920 kg/m³ |

## Experimental data

`tracesbar1.mat` holds **9 batches**. Each batch is a matrix of 3 pressure traces in bar, one per motor firing at low, mid and high pressure (LP, MP, HP), sampled at 1 kHz. That gives 27 firings in total.

> In batch 7 the MP and HP columns are swapped. `timeIntervals.m` and `FinalPlots.m` correct for this.

## Workflow (`main3.m`)

| Step | File | What it does |
|---|---|---|
| 0 | `main3.m` | Defines the grain geometry, throat areas and propellant mass, then loads `tracesbar1.mat`. |
| 1 | `plotLoadedData.m` | Quick-look plots of batch 1's traces, plus a zoom on the ignition transient. |
| 2 | `timeIntervals.m` | For each of the 27 traces it finds `p_max` (skipping the first 200 samples of transient) and the **action time** `[tA, tD]` at 5 % of `p_max`. It sets `p_ref` to half the mean action-time pressure and finds the **burning time** `[tB, tC]` at `p_ref`. It then computes `p_eff` (mean pressure over the burning time), `r_b = wt / t_b` and `c* = ∫p dt · A_t / M_p`. |
| 3 | `main3.m` | Takes the mean `c*`, the sample standard deviation `Sigma_c` (Bessel correction) and `Is = c*/g0`. |
| 4 | `UncertaintyScript.m` | Least-squares fit of `ln r_b = n ln p + ln a` over all 27 points. It returns `a`, `n`, their uncertainties `Inc_a` and `Inc_n`, and `R²`. |
| 5 | `MotorGiorgio.m` | Ballistic simulation of the **MP motor** (see below). |
| 6 | `FinalPlots.m` | Report-quality summary figures (see below). |

### Motor simulation (`MotorGiorgio.m`)

The simulation uses the auxiliary variable `z(t) = ∫ p^n dt`. Burned web distance is then `e = a·z`, and the quasi-steady mass balance becomes a single ODE:

```
dz/dt = p^n = ( ρ_p · a · c* · A_b(z) / A_t )^( n / (1-n) )
A_b(z) = 2π(Ri0 + a z)(L0 - 2 a z) + 2π(Re0² - (Ri0 + a z)²)     while a z < wt
```

It is integrated with a fixed-step RK4 (`dt = 1 ms`, 0 to 10 s). `p(t)`, `A_b(t)` and `r_b(t)` are recovered from the solution, and the results are compared against the experimental MP averages.

## Figures produced by `FinalPlots.m`

1. **Pressure traces**: 3×3 grid, one tile per batch, LP/MP/HP overlaid.
2. **Vieille law**: true log-log and linear views with the fit and ±1σ band. The linear view also shows the per-level mean ± std.
3. **c\* statistics**: per-test scatter, a violin + IQR + strip plot per pressure level, and a histogram with a Gaussian overlay.
4. **Motor simulation (MP)**: `p(t)`, `A_b(t)`, `r_b(t)`, web burn distance `e(t)`, and the `A_b`–`p` phase portrait coloured by time.

## File index

| File | Type | Role |
|---|---|---|
| `main3.m` | script | Entry point: inputs, data loading, calls every step |
| `plotLoadedData.m` | script | Raw-trace quick-look plots |
| `timeIntervals.m` | script | Time intervals, `p_eff`, `r_b`, `c*` per test |
| `UncertaintyScript.m` | script | Vieille-law regression and uncertainties (used by `main3.m`) |
| `Uncertainty.m` | function | Standalone version of the same regression: `[a, Inc_a, n, Inc_n, R2] = Uncertainty(p, rb)` |
| `MotorGiorgio.m` | script | RK4 ballistic simulation of the MP motor |
| `FinalPlots.m` | script | Final summary figures |
| `PlotCstar.m` | script | Standalone `c*`-per-test figure (run after `main3.m`) |
| `bonus.mlx` | live script | Bonus part of the assignment |
| `tracesbar1.mat` | data | Experimental pressure traces (9 batches × 3 firings) |

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox, used by `FinalPlots.m` (`ksdensity`, `quantile`, `normpdf`)
