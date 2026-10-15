# SGDrivE v0.1 development plan

Updated version of `plan_SGDRIVE.rtf` (2026-10-08), revised to match the
decisions in `DECISIONS.md`. The model itself is specified in
`MODEL_SPEC.md`; milestones and completion criteria are in `ROADMAP.md`;
coding rules are in `CONVENTIONS.md`.

Main changes from the earlier plan:
* genotype lifespan effect scales physiological lifespan X_A, not the
  mortality hazard (D-002);
* default sub-stages at most 10 (D-003);
* wound-search competitors are searching mated females by default (D-006);
* all populations tracked as counts per node (D-009);
* name, licence, repository and tooling fixed (D-001, D-007, D-008).

## 1. Architecture

```
R layer (once per run)                     C++ layer (every step)
─────────────────────                      ──────────────────────
sg_params()   parameters + units + checks  model_view.h   flat read-only model
sg_cube()     cube → tensors, ω → ℓ_g      rates.cpp      ONE rate kernel (aggregates → flows)
sg_forcing()  weather → daily vectors      ode_rk45.cpp   Dormand–Prince, per-day integration
sg_releases() release table                cb_step.cpp    chain-binomial stepper
sg_model()    assemble + index maps        events.cpp     releases at day boundaries
sg_init()     initial state / burn-in      summary.cpp    on-the-fly aggregation
simulate_ode(), simulate_cb()              checks.cpp     invariants, flow accounting
summaries, plots, check_model(), calibrate_wound_density()
```

Everything done once per run (validation, tensors, forcing) is in R;
everything done every step (rates, integration, random draws, summaries) is
in C++. Random numbers use R's generator; replicates run in parallel R
processes with L'Ecuyer streams and stored seeds.

## 2. Built-in checks

* **Inputs:** cube dimensions and names, Σ_k ih = 1, entries ≥ 0, φ and ξ in
  [0, 1], ω ≥ 0; parameter ranges and units; weather completeness and
  plausibility; release genotype, stage, sex and timing.
* **Runtime (C++):** no NaN/Inf; non-negative counts (stochastic), no
  meaningful negatives (ODE); a cap on extreme rates; flow accounting
  (change in total = inflows − outflows); ODE step failures reported.
* **`check_model()`:** largest (q + μ)·dt, approximate R₀ at a reference
  temperature, number of compartments, memory estimate, options in use.
* **Tests:** unit tests, regression against the legacy Gutierrez & Ponti R
  translations (`data-raw/reference/`), limit cases (one genotype, neutral
  Mendelian, MGDrivE2-like limit), reproducibility, snapshot tests.

## 3. Documentation

1. roxygen2 for every exported function: purpose, every argument with units
   and allowed values, return structure, runnable example, links to
   vignettes.
2. Parameter reference dataset `sg_param_table`: symbol, R name,
   description, units, default, allowed range, source tag. Printed by
   `print(sg_params())`.
3. Code comments cite equations (S1, S2, … in `MODEL_SPEC.md`), headers
   describe array layout, switches carry source tags, comments explain why.
4. Vignettes: getting started; model description; parameters, units and
   sources; genetic control with MGDrivE cubes; temperature and seasonality;
   stochastic simulation; validation and benchmarks; developer guide.
5. pkgdown site, README, NEWS, CITATION.

## 4. Development practices

* Imports: Rcpp, stats, utils. Suggests: MGDrivE, deSolve, testthat, knitr,
  rmarkdown, ggplot2, bench, future.apply.
* S3 classes with print/summary/plot methods; no modification of inputs; no
  global state; informative error messages.
* styler, lintr, clang-format; coverage measured locally (`dev/coverage.R`);
  benchmarks with bench.
* Semantic versioning; outputs record package version, model hash and seed.
* Example data in `inst/extdata`: gap-filled McAllen weather.
