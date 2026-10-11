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

Run benchmarks with the computer plugged in and awake.
