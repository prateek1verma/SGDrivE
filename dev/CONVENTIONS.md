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
| C++ files | one responsibility per file | `rates.cpp`, `ode_rk45.cpp` |
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

## Errors, warnings and messages

* Validate inputs in R, before any C++ call.
* Error messages say what is wrong, where, and how to fix it, e.g.
  `"cube$ih slices must sum to 1 over offspring genotypes; slice [AA, Aa, ] sums to 0.97."`
* Warnings are for valid but risky settings (e.g. a large chain-binomial step).
* Functions never modify their inputs and never use global state.

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

* Stochastic tests fix the seed with `withr::local_seed()` or `set.seed()`
  and test distributional properties, not exact draws, unless testing
  reproducibility.
* Slow tests (> 5 s) are skipped on CRAN with `skip_on_cran()`.

## Version control

* `main` is always releasable: `devtools::check()` gives 0 errors, 0
  warnings, 0 notes.
* Work on a branch, open a pull request, merge when CI passes.
* Commit messages: short imperative summary line, blank line, details.
