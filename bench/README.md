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
| `m4_cb_speed.R` | Chain-binomial solver (M4): seconds per simulated year and replicate of `simulate_cb()` from a one-year burn-in, for five scenarios, two node areas and two time steps, and the speed-up of the C++ step over the R reference step. Results in `results/m4_cb_speed.csv`. |
| `m4_dt_convergence.R` | Time step of the chain-binomial solver (M4): deterministic map of the reference step against the ODE for dt from 1 to 0.05 day, six scenarios over a year. Results in `results/m4_dt_convergence.csv`. |
| `phase2_cb_nodes.R` | Chain-binomial solver with several nodes (v0.2 phase 2): seconds per simulated year and replicate for 1, 10 and 50 nodes in a row, with and without movement, starting in every node or in the first only. Results in `results/phase2_cb_nodes.csv`. |
| `m3_cpp_vs_reference.R` | C++ backend (M3) against the R backend on the same scenarios and tolerances: run times (C++ Dormand-Prince day loop, C++ rates with lsoda, R), speed-ups, time per rate evaluation and the largest difference between backends, plus the full 24-year gp2014 run. Results in `results/m3_cpp_vs_reference.csv`. |

Run benchmarks with the computer plugged in and awake.

## Latest results (M3, Apple silicon, R 4.6.1)

At rtol 1e-8, atol 1e-10:

| Scenario | States | R (s) | C++ (s) | Speed-up |
|---|---|---|---|---|
| wild, gp2014 preset, 1 year | 150 | 8.07 | 0.21 | 38x |
| wild, defaults, 1 year | 49 | 3.23 | 0.04 | 74x |
| SIT W/S, weekly releases, 1 year | 118 | 2.56 | 0.06 | 42x |
| Mendelian (3 genotypes), 1 year | 207 | 4.18 | 0.14 | 30x |
| homing drive (10 genotypes), 90 days | 1390 | 2.36 | 0.18 | 13x |
| wild, gp2014 preset, 24 years | 150 | about 480 | 4.77 | about 100x |
| SIT W/S, weekly releases, 24 years | 118 | not run | 1.20 | |

The backends agree to about 1e-8 (relative). The C++ kernel costs about the
same per state in every scenario (17-18 ns per state and evaluation); the
speed-up is smaller for large models because the vectorised R code spreads
its fixed per-call overhead over more states.

## Chain-binomial solver (M4, Apple silicon, R 4.6.1)

Seconds per simulated year and replicate (`simulate_cb()`, eggs in batches,
one-year burn-in, area 1000 km^2):

| Scenario | States | dt = 0.1 | dt = 0.05 | C++ vs R reference step |
|---|---|---|---|---|
| wild, defaults | 49 | 0.074 | 0.102 | 55x |
| wild, gp2014 preset | 150 | 0.34 | 0.41 | 27x |
| SIT W/S, weekly releases | 118 | 0.063 | 0.084 | 30x |
| Mendelian (3 genotypes) | 207 | 0.099 | 0.132 | 31x |
| homing drive (10 genotypes) | 1390 | 0.33 | | 23x |

With 100 times more individuals (area 1e5 km^2, up to 11 million adults) the
time per year grows by only about 10-35% (0.081 s and 0.092 s for the wild
and SIT models): the cost depends on the number of compartments in use, not
on the number of individuals.

## Chain-binomial solver with several nodes (v0.2 phase 2, Apple silicon, R 4.6.2)

Seconds per simulated year and replicate (`simulate_cb()`, Mendelian, nodes
of 50 km^2 in a row on one forcing series, `record = "totals"`):

| Nodes | No movement, all occupied | Movement, all occupied | Movement, flies start in node 1 |
|---|---|---|---|
| 1 | 0.095 | | |
| 10 | 0.40 | 0.45 | 0.42 |
| 50 | 1.71 | 2.11 | 1.33 |

Time grows about linearly with the number of occupied nodes; nodes on one
forcing series share their chains, and empty nodes cost about a third of an
occupied one. Larval crowding rebuilds the larval chain per node and step
(+80% for ten nodes; D-032).
