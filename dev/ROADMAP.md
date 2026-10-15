# SGDrivE roadmap to v0.1.0

| Milestone | Deliverables | Done when |
|---|---|---|
| **M0 Foundations** | Package skeleton (DESCRIPTION, GPL-3, Rcpp/C++17), C++ toolchain check, testthat, lintr/styler/clang-format, GitHub Actions (R CMD check on 3 OSes, lint), README/NEWS/CITATION, `dev/` design documents, reference data | `devtools::check()` 0/0/0; CI passing; documents approved |
| **M1 Inputs** | `sg_params()` with units table (`sg_param_table`), `sg_cube()` (validation, tensors, ω → ℓ_g), `sg_forcing()` (Δx, μ_T, μ_E, φ_ov; NOAA reader; gap filling), `sg_releases()`, `sg_model()` (index maps) | Input tests pass; bad inputs give clear errors; McAllen cold index reproduced (1971 = 6.67, 1979 = 9.80) |
| **M2 Reference model in R** | Pure-R genotype-resolved ODE (deSolve), the "truth" for later tests | With one genotype and n = 25, i0 = 8, matches the legacy `gp_ode_simulate()` (`data-raw/reference/`) to 1e-6; neutral Mendelian cube stays at Hardy–Weinberg |
| **M3 C++ ODE** | Rate kernel, RK45 solver, `simulate_ode()`, releases | Matches M2 to 1e-6; ≥ 10× faster; 24 years, 2 genotypes in about 1 s |
| **M4 C++ chain-binomial** | `simulate_cb()`, replicates, parallel runs, seeds | Stochastic mean converges to the ODE at large sizes; same seed gives the same result; counts never negative |
| **M5 Checks and outputs** | `check_model()`, `sg_r0()`, flow accounting (`sg_flow_balance()`), runtime guards, `calibrate_wound_density()`, output tables (stage × sex × genotype, allele frequencies, female classes, case proxy from fertile oviposition), ggplot2 plots; extinction and replicate statistics are left to users' own analyses of the raw replicate outputs (D-026) | All runtime checks tested; calibration hits target adult numbers |
| **M6 Scientific validation** | Mosquito-like limit vs MGDrivE2fast; SIT, pgSIT, fsRIDL on McAllen; Knipling and Gutierrez 2019 release checks; lumped model for R₀ and release thresholds | Validation vignette reproducible |
| **M7 Documentation and release** | Vignettes, pkgdown site, README, NEWS, CITATION, benchmarks | Coverage ≥ 85%; 0/0/0 check; tag v0.1.0 |

## After v0.1

| Version | Adds |
|---|---|
| v0.2 | Multiple nodes and fly movement |
| v0.3 | Host and wound model; modelled case counts |
| v0.4 | Host movement, larval transport, quarantine and treatment |
| later | Explicit gonotrophic cycle; male mating capacity |
