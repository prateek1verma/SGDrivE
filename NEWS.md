# SGDrivE 0.0.0.9000

## Milestone M4: chain-binomial stochastic solver (in progress)

* Design of the chain-binomial solver (D-022; model specification S15-S19):
  exact moves along sub-stage chains within a step, mating before ageing,
  a canonical order of random draws shared by an R reference and the C++
  step, a deterministic map for testing, per-replicate random-number
  streams.
* The C++ rate kernel is split into day rates and state aggregates, shared
  by the ODE and chain-binomial solvers; ODE results are unchanged.
* Exact advancement probabilities along a chain over a time step, with
  death in the chain, and their time average over the step for arrivals
  (`src/chains.h`, uniformization; an independent R version by matrix
  exponential), tested against closed forms.
* Reference chain-binomial step in pure R (stochastic and deterministic
  modes; internal, the test oracle for the C++ step). The deterministic map
  matches the ODE derivatives for short steps, and its mean equals the
  stochastic step. Time-step study (`bench/m4_dt_convergence.R`): placing
  arrivals at their expected position within the step brings the error at
  dt = 0.1 from 3-5% to below 0.6% in every block; default dt = 0.1.
* Egg laying in batches (D-024): new parameters `oviposition` ("batches",
  default, or "eggs") and `clutch_size` (default R * phi_ov, the eggs of one
  day). In the chain-binomial solver, oviposition events are drawn first
  and eggs per event after, which keeps the mean (and the ODE) but gives
  the variance of real egg masses; the `egg_batches` counter counts events.
* C++ chain-binomial step and day loop (`src/cb_step.h`, `src/cb_step.cpp`):
  identical, draw for draw, to the R reference with the same seed, and its
  deterministic map to 1e-12. Releases at the start of their day, exact
  integer flow accounting, early stop once a population is extinct with no
  releases left. About 0.06 s per simulated year for a two-genotype release
  model, nearly independent of population size.

## Milestone M3: C++ backend

* `simulate_ode()` now uses the C++ backend by default (`backend = "cpp"`,
  `method = "dopri5"`; D-021). The pure-R backend stays available as the
  readable reference (`backend = "R"`). The C++ backend needs 'deSolve' only
  if a day falls back to `lsoda` (or for a deSolve `method`).
* Documentation: solver description in `?simulate_ode` and the model
  specification (section 10), C++ coding conventions, README, benchmarks
  (24 years of a two-genotype release model in about 1 s).
* New manually triggered `sanitizers` workflow: tests with AddressSanitizer
  and UndefinedBehaviorSanitizer for the C++ code.
* Internal model view (`R/model-view.R`, `src/model_view.h`): everything the
  C++ backend needs, flattened once per simulation with 0-based indices,
  column-major matrices and inverse adult sub-stage widths, and validated on
  the C++ side (D-021).
* C++ rate kernel (`src/rates.h`, `src/rates.cpp`): a block-by-block
  translation of the R reference rate function, allocation-free per
  evaluation. Derivatives, counter rates and the largest-rate diagnostic
  equal the R reference to round-off (about 1e-16; tested to 1e-12) for every
  option, preset and cube type, including sex-linked cubes.
* Dormand-Prince 5(4) integrator (`src/dopri5.h`): Hairer's PI step-size
  control, starting step and stiffness test; first-same-as-last (6
  evaluations per step); tested for accuracy on problems with known
  solutions, 5th-order convergence and stiffness detection.
* `simulate_ode(backend = "cpp")`: the C++ backend. The whole day loop
  (releases, Dormand-Prince integration with daily counters, diagnostics)
  runs in C++ (`src/simulate.cpp`); a day on which stiffness is detected is
  integrated again with `lsoda` using the C++ rates, and the C++ loop resumes
  on the next day (`stiff_fallback`, default `TRUE`). Any deSolve method can
  also be driven by the C++ rates (`method = "lsoda"`, `"radau"`, ...).
  About 25-30 times faster than the R backend on a year of McAllen weather,
  agreeing to the integration tolerance (about 1e-11 at `rtol = 1e-10`).
  `method` now defaults to `"dopri5"` (C++) or `"lsoda"` (R).
* Solver diagnostics gain `rejected` (Dormand-Prince rejected steps) and
  `fallback` (day integrated by the `lsoda` fallback); `meta` records
  `stiff_fallback` and the number of fallback days, and `print()` reports
  them.
* Backend equivalence: every M2 simulation, biology and stress test now runs
  on both backends, and a new test compares them directly over ten scenarios
  (options, releases of every kind, lifespan costs, the gp2014 preset,
  MGDrivE homing-drive and sex-linked cubes) to 1e-8.
* The C++ backend reproduces the legacy Gutierrez & Ponti (2014) ODE over the
  full 24 years (1960-1983) to 6e-9 in a routine test (4.5 s; new fixture
  of totals every 7th day), without needing the stiff fallback on any day.
* `bench/m3_cpp_vs_reference.R`: the C++ backend is 12-73 times faster than
  the R backend at the same tolerances (24 years of McAllen weather: 4.5 s
  instead of about 8 minutes).

## Milestone M2: reference model

* `sg_init()` builds initial states: `"constant"` counts per sub-stage split
  by genotype weights (defaults reproduce the Gutierrez & Ponti starting
  state), `"hardy_weinberg"` from allele frequencies, or a full `"state"`
  vector. Adults go only to genotypes of the matching sex and to adult classes
  they reach; mated females are split across mates by random mating (D-020).
* Internal reference rate function (pure R, model specification S7-S14):
  all state derivatives and cumulative-counter rates (eggs laid by genotype,
  egg batches, emergences, matings by female x mate genotype, deaths by
  block, deaths by old age) for one day's forcing; the test oracle for the
  C++ kernel.
* `simulate_ode(model, init, backend = "R")` solves the model with
  'deSolve', one day at a time (exact daily breakpoints; releases added as
  jumps at the start of their day). Returns an `sg_sim` object with daily
  states (or every `output_every` days), daily flow totals, releases per
  block and solver diagnostics (steps, evaluations, stiff-method use,
  largest rate, smallest state value).
* Milestone checks of the reference solver: with one wild genotype and the
  `"gp2014"` preset it reproduces the legacy Gutierrez & Ponti ODE to a
  relative error below 1e-6 (1960 in routine tests; 1960-1983 as a slow test,
  `SGDRIVE_SLOW_TESTS=true`, also available as a manual CI workflow); a
  neutral Mendelian cube stays at Hardy-Weinberg with constant allele
  frequency; rate mating approaches instant mating as nu grows; sterile-male
  releases divert matings and suppress the wild population; genotype
  lifespan costs act as specified; sex-linked cubes keep sexes consistent.
* Stress tests: very large sterile-male releases (10x to 1000x the wild
  males), larval crowding with a huge fertile egg release, a very large
  mating rate and many sub-stages; lsoda and radau agree within 1e-6 and
  states stay non-negative. Release size adds no stiffness.
* `bench/m2_reference_baseline.R`: runtime baseline of the reference solver
  (results in `bench/results/`), the reference for the C++ backend (M3).

## Milestone M1: inputs

* `sg_params()` creates validated model parameters, with presets
  `"default"` and `"gp2014"` (the Gutierrez & Ponti 2014 code
  configuration); `update()`, `print()` and `as.data.frame()` methods.
* `sg_param_table()` returns the definition of every parameter: default,
  allowed range, units, source and model equation.
* `sg_cube()` reads and validates MGDrivE inheritance cubes (including
  sex-linked cubes, whose impossible crosses have zero inheritance slices),
  converts omega to lifespan fractions and precomputes the egg tensor
  s * ih * tau.
* `sg_omega_to_ell()` inverts `MGDrivE::calcOmega()`.
* Temperature functions of Gutierrez & Ponti (2014): `sg_degree_days()`,
  `sg_temp_mortality()`, `sg_egg_viability()`, `sg_ovip_scalar()`
  (identical to the legacy code under the `"gp2014"` preset).
* `sg_forcing()` checks daily weather and precomputes the daily forcing
  terms (degree-days, mortality, egg viability, oviposition scalar, moisture
  index).
* Weather tools: `sg_example_weather()` (McAllen, Texas, 1960-1983),
  `sg_read_ghcn()`, `sg_fill_weather()` (Gutierrez et al. 2019 gap-filling
  procedure) and `sg_cold_index()`, which reproduces the published McAllen
  cold-stress index (1971 = 6.67, 1979 = 9.80).
* `sg_releases()` defines release schedules of any genotype as eggs, pupae
  or unmated adults (male, female or both), on day numbers or dates, per
  node or per km2; schedules combine with `c()`.
* `sg_model()` assembles parameters, cube, forcing and releases into a model:
  segmented adult chains with exact maturity thresholds for every genotype
  (D-018), state-vector layout, counts per node, release events mapped to
  state positions, and an MD5 fingerprint of all inputs. `sg_state_index()`
  and `sg_state_names()` locate and name states.
* New parameter `X_ovi` (age at first oviposition); the class overrides
  `i0`, `i_m` and `i_mM` are removed (D-018).

## Milestone M0: foundations

* Package skeleton: DESCRIPTION, GPL-3 licence, Rcpp/C++17 build setup.
* `sg_backend_info()`: toolchain check reporting how the C++ backend was
  compiled. No model code yet.
* Testing (testthat 3rd edition), linting (lintr, styler, clang-format) and
  GitHub Actions workflows (R CMD check on macOS, Windows and Ubuntu; lint).
* Design documents in `dev/`: decisions log, v0.1 model specification,
  roadmap, conventions.
* Reference outputs of the legacy Gutierrez & Ponti (2014) R translations in
  `data-raw/reference/` for later regression tests.
