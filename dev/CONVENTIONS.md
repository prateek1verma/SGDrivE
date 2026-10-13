# SGDrivE coding conventions

These conventions apply to all code in the package. `dev/` and `data-raw/`
are excluded from the built package.

## Naming

| Item | Convention | Example |
|---|---|---|
| Exported R functions | `sg_` prefix, snake_case | `sg_params()`, `sg_cube()` |
| Simulation entry points | verb first, no prefix | `simulate_ode()`, `simulate_cb()` |
| Internal R helpers | leading dot or no export | `.check_cube_dims()` |
| S3 classes | `sg_` prefix | `sg_params`, `sg_model`, `sg_sim` |
| C++ functions exported to R | `sg_..._cpp` (called only from R wrappers) | `sg_backend_info_cpp()` |
| C++ files | one responsibility per file | `rates.cpp`, `dopri5.h` |
| Model symbols in code | follow `dev/MODEL_SPEC.md`, ASCII | `q_A`, `mu_T`, `ell_g`, `i0` |

R code follows the tidyverse style guide (apply with `Rscript dev/style.R`), checked by
`lintr` with a 100-character line limit. C++ follows `.clang-format`
(Google-based, 2-space indent, 100 columns).

## Units

All quantities carry fixed units. Every parameter's units are documented in
the parameter table (milestone M1) and in its roxygen `@param` entry.

| Quantity | Unit |
|---|---|
| Time | calendar day (d) |
| Physiological time | degree-day (dd), as defined by Gutierrez & Ponti's Δx(T) |
| Temperature | °C |
| Population states | counts per node (individuals) |
| Rates and hazards | d⁻¹ |
| Node area | km² (default 1, so counts equal densities numerically) |
| Densities (output only) | individuals km⁻² = counts / area |

Daily probabilities from the literature are converted to hazards with
`mu = -log(survival)` before use. Never mix probabilities and rates.

## Array layout and index order

State blocks are stored in one contiguous numeric vector. Index order, from
slowest to fastest varying:

```
[node, host_class, genotype, (mate_genotype), sub_stage]
```

* `node` and `host_class` exist from v0.1 but have length 1.
* `mate_genotype` exists only for mated females `F`.
* `sub_stage` is the Erlang sub-stage (age class) within a life stage.

An index map built in R (`sg_model()`) gives the offset and dimensions of each
block. C++ reads blocks through that map and never hard-codes offsets.
R uses 1-based and C++ 0-based indices; conversion happens only at the R/C++
boundary.

## C++ code

* Files in `src/`, one responsibility each: `model_view.*` (the model
  flattened for C++), `rates.*` (rate kernel), `dopri5.h` (integrator),
  `simulate.cpp` (day loop and R entry points), `backend_info.cpp`.
  Test-only entry points live in `*_test.cpp`.
* Plain C++17, Rcpp and `std::vector`; no other dependencies. Code lives in
  `namespace sgdrive`; only `sg_..._cpp` functions are exported to R (run
  `Rcpp::compileAttributes()` after changing them) and they are called only
  from internal R wrappers.
* C++ reads the model only through the model view built in R
  (`.sg_model_view()`), validated once in the `ModelView` constructor; it
  never re-derives parameters or hard-codes offsets.
* The rate kernel allocates no memory per evaluation: work arrays are
  allocated once when it is constructed. Objects that live across R calls
  are held in an `Rcpp::XPtr` handle.
* Long loops call `Rcpp::checkUserInterrupt()` once per simulated day.
  Errors use `Rcpp::stop()` with a message that names the problem.
* Any change to the model goes into the R reference rate function and the
  C++ kernel together, with an equivalence test (`test-rates-cpp.R`,
  `test-backend-equivalence.R`); the R code is the readable specification.
* Format with `.clang-format`. After larger changes to `src/` run
  `bash dev/sanitize.sh` (AddressSanitizer and UBSan, macOS, a few minutes);
  before a release also run the manual `sanitizers` GitHub workflow (Linux,
  R-devel).

## Errors, warnings and messages

* Validate inputs in R, before any C++ call.
* Error messages say what is wrong, where, and how to fix it, e.g.
  `"cube$ih slices must sum to 1 over offspring genotypes; slice [AA, Aa, ] sums to 0.97."`
* Warnings are for valid but risky settings (e.g. a parameter that has no
  effect with the chosen options).
* Functions never modify their inputs and never use global state; the only
  exception is R's random-number state (next section).

## Random numbers

* All random draws use R's own generator: `stats::rbinom`/`rpois` in R and
  the same Rmath functions (`R::rbinom`, `R::rpois`) in C++, called only
  inside an Rcpp function (its `RNGScope` handles the generator state). No
  other generator, no C++ threads; parallel runs use R processes.
* The order of draws in a chain-binomial step is part of the specification
  (`dev/MODEL_SPEC.md`, section 10). Any change to the step goes into the R
  reference step (`R/cb-reference.R`) and the C++ step (`src/cb_step.h`)
  together and must keep the draw-for-draw test (`test-cb-cpp.R`) passing;
  multinomials are sequential conditional binomials in category order.
* Replicates run on their own L'Ecuyer-CMRG streams (`.sg_with_stream()`).
  A function that uses random numbers draws its seed from the session's
  generator only when the user gives none, and restores the session's
  generator kind and state afterwards.
* Stochastic tests fix the seed and test distributions (means, variances,
  frequencies) with bounds of about five standard errors, or compare with
  the deterministic map (draws replaced by expectations).

## Comments and documentation

* Every exported function has roxygen documentation: purpose, every argument
  with units and allowed values, return value, and a runnable example.
* Comments explain *why*, not *what*.
* Each block in the C++ rate kernel cites its equation in `dev/MODEL_SPEC.md`
  (later the model-description vignette), e.g. `// Eq. (S7): egg production`.
* Each optional model switch is tagged with its source:
  `[GP2014]`, `[G2019]`, `[MGDrivE]`, `[MGDrivE2]` or `[adapted]`.

## Testing

* testthat 3rd edition; one test file per R source file where practical.
* Shared tolerances are defined in `tests/testthat/helper-tolerance.R`:

| Name | Value | Use |
|---|---|---|
| `sg_tol_exact` | 1e-12 | identities (probabilities sum to 1, flow accounting) |
| `sg_tol_ode` | 1e-6 | two solvers of the same model (C++ vs pure-R reference), relative |
| `sg_tol_reference` | 1e-2 | agreement with the legacy Gutierrez & Ponti R translations |

* Stochastic tests fix the seed with `set.seed()` and test distributional
  properties, not exact draws, unless testing reproducibility or R/C++
  equality (see Random numbers).
* Slow tests (> 5 s) are skipped on CRAN with `skip_on_cran()`.

## Version control

* `main` is always releasable: `devtools::check()` gives 0 errors, 0
  warnings, 0 notes.
* Work on a branch, open a pull request, merge when CI passes.
* Commit messages: short imperative summary line, blank line, details.
