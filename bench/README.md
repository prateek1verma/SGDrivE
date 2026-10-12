# bench

Performance benchmarks (excluded from the built package). Added from
milestone M3 onwards using the `bench` package, comparing the C++ solvers
with the pure-R reference model and, for the mosquito-like limit, with
MGDrivE2fast. Results are recorded per version so speed regressions are
visible.

## Scripts

| Script | Measures |
|---|---|
| `m2_reference_baseline.R` | Pure-R reference ODE solver (M2): seconds per simulated year, steps per day and time per rate evaluation for five scenarios (1 to 10 genotypes). Results in `results/m2_reference_baseline.csv`. |
| `m3_cpp_vs_reference.R` | C++ backend (M3) against the R backend on the same scenarios and tolerances: run times (C++ Dormand-Prince day loop, C++ rates with lsoda, R), speed-ups, time per rate evaluation and the largest difference between backends, plus the full 24-year gp2014 run. Results in `results/m3_cpp_vs_reference.csv`. |

Run benchmarks with the computer plugged in and awake.

## Latest results (M3, Apple silicon, R 4.6.1)

At rtol 1e-8, atol 1e-10:

| Scenario | States | R (s) | C++ (s) | Speed-up |
|---|---|---|---|---|
| wild, gp2014 preset, 1 year | 150 | 7.83 | 0.21 | 38x |
| wild, defaults, 1 year | 49 | 3.15 | 0.04 | 73x |
| SIT W/S, weekly releases, 1 year | 118 | 2.48 | 0.06 | 41x |
| Mendelian (3 genotypes), 1 year | 207 | 4.03 | 0.13 | 30x |
| homing drive (10 genotypes), 90 days | 1390 | 2.23 | 0.18 | 12x |
| wild, gp2014 preset, 24 years | 150 | about 480 | 4.50 | about 100x |

The backends agree to about 1e-8 (relative). The C++ kernel costs about the
same per state in every scenario (17-18 ns per state and evaluation); the
speed-up is smaller for large models because the vectorised R code spreads
its fixed per-call overhead over more states.
