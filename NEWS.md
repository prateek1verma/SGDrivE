# SGDrivE 0.0.0.9000

## Milestone M2: reference model (in progress)

* `sg_init()` builds initial states: `"constant"` counts per sub-stage split
  by genotype weights (defaults reproduce the Gutierrez & Ponti starting
  state), `"hardy_weinberg"` from allele frequencies, or a full `"state"`
  vector. Adults go only to genotypes of the matching sex and to adult classes
  they reach; mated females are split across mates by random mating (D-020).

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
