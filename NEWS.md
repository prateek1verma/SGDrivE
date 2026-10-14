# SGDrivE 0.0.0.9000

## Scenario scripts

* New folder `inst/scenarios/` (installed; `system.file("scenarios",
  package = "SGDrivE")`) with short user-style scripts for scientific testing
  and future tutorials. First script: `01_wild_type_dynamics.R`, a wild-type
  population at McAllen, Texas (model check, burn-in, deterministic and
  stochastic runs, plots, R0 by temperature).
* `SGDrivE.Rproj` RStudio project file for local development.

## Milestone M5: checks, outputs, calibration and plots

* Documentation: README section on outputs, plots, checks and calibration
  (with a figure), plotting conventions, ROADMAP and model specification
  updated.
* Output tables for ODE and chain-binomial simulations (D-026):
  `sg_totals()` (by stage, sex, genotype and mate genotype, per node or per
  km^2), `as.data.frame()` methods (long tables), `sg_allele_freq()`,
  `sg_females()` (all, unmated, mated, wild type, egg-laying and
  fertile-equivalent females) and `sg_cases()` (the case proxy by day, week,
  month or year, with the daily fraction of wounds newly infested and cases
  per host). Simulations now store the model information these need.
* Cases are wounds newly infested with viable eggs: an oviposition event
  counts when its batch holds at least one viable egg, so females mated to
  sterile males cause no cases (D-026). New counter `fertile_batches` in
  both solvers; in the chain-binomial step each pair's batches are split
  into fertile ones and their eggs (zero-truncated Poisson), an exact
  refinement that keeps the egg distribution, so random draws differ from
  earlier versions for the same seed.
* Plots with ggplot2 (suggested package): `plot()` methods for ODE and
  chain-binomial simulations return ggplot objects of totals by stage or
  genotype, allele frequencies, female classes or cases, with stochastic
  replicates drawn as separate lines; `sg_plot_data()` returns the plotted
  table. Colours are colour-blind friendly: Okabe-Ito by default,
  ColorBrewer (`palette = "brewer"`) or viridis; ColorBrewer Paired for 9 to
  12 series.
* `check_model()`: model size and options, the fastest sub-stage rate and
  the mating rate, R0 over the forcing (days with R0 > 1, longest run),
  genotypes that can never lay, mate or mature, sterile releases, and the
  seasonal low of the wild population (ODE after a burn-in) with a warning
  when stochastic runs may go extinct by chance.
* `sg_r0()`: eggs per egg at low density of a wild-type population at
  constant temperatures. `sg_flow_balance()`: the flow balance of any
  simulation between output days.
* The solver diagnostic `max_rate` includes the mating rate of unmated
  females (O-007).
* `calibrate_H()`: the wound density at which a summary of the population
  (by default the mean adult density over the model's forcing period, per
  km^2) equals a target. Without an Allee effect or releases the population
  is exactly proportional to H, so two model runs suffice (closed form);
  otherwise a root search on log H. Returns the calibrated parameters and
  model.
* Runtime guards: a deSolve day that does not reach its end (lsoda can
  return its starting values with a success code when its step size
  underflows) or a non-finite state stops the ODE with the day; a
  non-finite or negative count or flow stops a chain-binomial run.

## Milestone M4: chain-binomial stochastic solver

* Documentation: README section on stochastic simulation, random-number
  conventions (`dev/CONVENTIONS.md`), `simulate_cb()` in the reference
  index.
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
* `simulate_cb()`: the stochastic chain-binomial solver for users.
  Replicates with one L'Ecuyer-CMRG random-number stream each (replicate
  r is the same whether run alone, among others, or in parallel with
  `parallel = TRUE` through 'future.apply'); `seed`, or the session's
  generator, which is otherwise left untouched; `keep = "state"` or
  `"totals"`, `output_every`, a memory check, extinction days and `print()`.
  Non-integer initial states and releases are rounded without bias in each
  replicate.
* `sg_init()` gains `method = "burn_in"` (the ODE over the years before the
  start, from the weather or the model's first year repeated, with a check
  that the starting state no longer matters) and `method = "equilibrium"`
  (closed-form equilibrium of a wild-type population at constant
  temperature, with R0 reported; D-025).
* `simulate_cb()` also reports `wild_extinct_day`, the first day without
  individuals of the wild type: the measure of eradication when released
  sterile insects remain (`extinct_day` requires every compartment empty).
* Statistical validation of the chain-binomial solver: one-step
  distributions match binomial and multinomial theory (means, variances,
  covariances over 40,000 steps); the replicate mean converges to the
  deterministic map like 1/sqrt(N) and the map matches the ODE; neutral
  allele frequencies keep their mean while drift grows in small
  populations; sterile-male releases eliminate the wild type; egg batches
  raise extinction in small populations (22% against 1% at 3 km^2 through a
  McAllen winter). Stress tests: counts above 1e9, step sizes 0.05-1 day,
  fast mating, 25 sub-stages, a ten-genotype homing drive.
* `bench/m4_cb_speed.R`: 0.06-0.34 s per simulated year and replicate,
  nearly independent of population size; the C++ step is 23-55 times
  faster than the R reference step.

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
