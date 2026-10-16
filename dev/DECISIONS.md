# SGDrivE design decisions

Each entry records a decision, its date, the reasoning, and its effect on the
model. Decisions are changed only by adding a new entry that supersedes an
old one.

Status: **Accepted**, **Superseded** or **Open**.

---

## D-001 Package name: `SGDrivE` (Accepted, 2026-10-06)

"Screwworm Genetic-control Dynamics Explorer". Echoes MGDrivE ("Mosquito Gene
Drive Explorer") because the package reads MGDrivE cubes and uses MGDrivE2's
mating structure, but "genetic control" rather than "gene drive", since the
main use cases (SIT, pgSIT, fsRIDL) are self-limiting, not drives.
Name confirmed free on CRAN, Bioconductor and GitHub (2026-10-08).

## D-002 Genotype lifespan effect scales physiological lifespan only (Accepted, 2026-10-06)

*Supersedes an earlier draft that divided the mortality hazard by ℓ_g.*

Each genotype g has a lifespan fraction ℓ_g ∈ (0, 1+] that scales its adult
physiological lifespan:

    X_{A,g} = ℓ_g · X_A,   q_{A,g}(t) = n_A · Δx(T(t)) / (ℓ_g · X_A)

* Temperature mortality μ_A(T) and any background mortality are identical for
  all genotypes and are **not** multiplied by ℓ_g.
* Maturity thresholds (pre-oviposition, female receptivity, male maturity)
  stay in absolute degree-days, so each genotype has its own class threshold
  (see D-004 and MODEL_SPEC §5). A shorter lifespan does not mean earlier
  maturity. If a genotype's lifespan is shorter than the pre-oviposition
  period it never lays eggs.
* ℓ_g comes from the cube's MGDrivE v1 `omega` or is given directly. For a
  cube built with `MGDrivE::calcOmega(mu_ref, L)`, the exact inverse is
  ℓ_g = μ_ref / (1 − ω_g (1 − μ_ref)). ω = 0 maps to ℓ_g = 0.001.

Reasoning: the user wants a fixed genotype cost to total lifespan that is
independent of temperature. In this model adults also die of old age (ageing
chain), so scaling the hazard alone barely changes lifespan when senescence
dominates. Expected behaviour: realised lifespan ratio ≈ ℓ_g when senescence
dominates, moving towards 1 when cold or field mortality dominates.

## D-003 Default sub-stage numbers tuned for speed, at most 10 (Accepted, 2026-10-06)

Long stages get more sub-stages, short stages fewer; **no default exceeds 10**.
Proposed defaults (MODEL_SPEC §3): eggs 3, larvae 6, pupae 10, adults 10.
Users may set any number, e.g. 25 for every stage to reproduce Gutierrez &
Ponti (2014), which is used in regression tests.

## D-004 Mating defaults (Accepted, 2026-10-06)

* ν = ln 2 d⁻¹ (Gutierrez et al. 2019: half of virgins mate per day).
* h_M = 0 (no mate-finding Allee effect by default).
* Female receptivity and male maturity thresholds are user-specifiable;
  otherwise derived from the number of adult sub-stages using the
  nearest-class rule of MODEL_SPEC §5. Defaults correspond to mating from
  emergence (class 1), as in Gutierrez et al. (2019).
* A switch `mating = "instant"` reproduces MGDrivE2 (mating at emergence).

## D-005 Releases (Accepted, 2026-10-06)

MGDrivE2-like release events:
* life stage: eggs, pupae, or unmated adults (males and/or females);
* released individuals enter the **first sub-stage** of that stage;
* sex composition for adults: `"male"`, `"female"` or `"both"`.

## D-006 Wound-search denominator: searching mated females (Accepted, 2026-10-06)

The competitors in the wound-search term Ψ are **searching mated females**
(age ≥ pre-oviposition threshold), including females mated to sterile males,
who also search and lay non-viable eggs. "All adults" (Gutierrez & Ponti 2014
code) is an option, for benchmarking only.

Reasoning: with "all adults", released sterile males would compete for wounds
and suppress wild egg laying, overstating SIT. Males and immature females do
not seek wounds.

## D-007 Licence: GPL (>= 3) (Accepted, 2026-10-06)

Compatible with MGDrivE and MGDrivE2.

## D-008 Repository and tooling (Accepted, 2026-10-08)

* Location `~/Documents/GitHub/SGDrivE`, outside Box (cloud sync can corrupt
  git folders and compiled objects).
* Private GitHub repository `prateek1verma/SGDrivE`.
* Author and maintainer: Prateek Verma <prateekverma@berkeley.edu>.
* Minimum R version 4.1.
* Test coverage measured locally only (`dev/coverage.R`); no Codecov.
* C++17 via Rcpp.

## D-009 States are counts per node (Accepted, 2026-10-06)

All populations are counts per node of area A (km², default 1). Biological
quantities that depend on density (wound density, Allee half-saturation,
release rates per area) are given per km² and converted to per-node counts
when the model is assembled. Counts are required by the chain-binomial solver
and make eradication ("< 1 individual") well defined.

## D-010 No Petri-net layer; one C++ rate kernel (Accepted, 2026-10-06)

The model is coded directly on structured arrays. One C++ rate kernel feeds
both the ODE solver (RK45, integrated one day at a time) and the
chain-binomial solver. MGDrivE2fast showed that speed comes from compiled
rate evaluation; the Petri net adds no benefit here and would make NWS
features (adult ageing chains, age kept on mating, wound-limited egg laying)
harder to express.

## D-011 Offspring viability τ applies at egg production (Accepted, 2026-10-08; was O-001)

τ[g,m,k] multiplies egg production (MODEL_SPEC S9), as in MGDrivE2. Non-viable
eggs are never created.

## D-012 Converting cube ω needs μ_ref (Accepted, 2026-10-08; was O-002)

`sg_cube()` requires `mu_ref` (the daily adult mortality used when the cube's ω
values were computed with `MGDrivE::calcOmega()`) whenever any ω ≠ 1, unless
lifespan fractions `ell` are given directly. The error message explains why.

## D-013 "both"-sex releases split 0.5 by default (Accepted, 2026-10-08; was O-003)

`sg_releases(sex = "both")` splits each release equally between males and
females; the `split` argument (fraction female) changes this.

## D-014 Egg development temperature: ambient by default, host option (Accepted, 2026-10-08; was O-005)

Default `egg_temperature = "ambient"` ages eggs with ambient degree-days, as
in the Gutierrez & Ponti code. Option `"host"` uses a constant host body
temperature `T_host` = 38.5 °C (Gutierrez & Ponti 2014, citing Rubink 1987),
so q_E = n_E Δx(T_host) / X_E is constant.

## D-015 Simulation time origin (Accepted, 2026-10-08)

Day 1 of a simulation is the first date of the forcing. Release times may be
given as day numbers or as `Date`s (converted using the forcing dates).

## D-016 Base-R error handling (Accepted, 2026-10-08)

Errors and warnings use base R (`stop(call. = FALSE)`, `warning()`), with
messages that say what is wrong and how to fix it. No `cli`/`rlang`
dependency.

## D-017 Parameter table defined in R code, returned by `sg_param_table()` (Accepted, 2026-10-08)

The parameter table (names, defaults, ranges, units, sources, spec equations)
is defined once, in `R/param-table.R`, and returned by the exported function
`sg_param_table()`, instead of being shipped as a lazily loaded dataset. It is
the single source for defaults, validation, `print()` and documentation, and
changes to it are visible in code review.

## D-018 Segmented adult chains with exact maturity thresholds (Accepted, 2026-10-08)

*Replaces the nearest-class rounding rule for maturity thresholds (former
MODEL_SPEC S6) and the class overrides `i0`, `i_m`, `i_mM`.*

Each adult chain is split at its maturity thresholds into segments whose mean
durations are exact:

* females (U and F share one chain because age is kept on mating): boundaries
  at X_rec (receptivity) and X_ovi (first oviposition); males: boundary at
  X_mat (maturity). Zero-length segments are dropped; `mating = "instant"`
  drops X_rec.
* The total number of adult sub-stages stays n_A for each sex (D-003). It is
  allocated to segments in proportion to the wild-type segment lengths, at
  least one each (largest-remainder rounding); `sg_model(adult_substages = )`
  overrides the allocation per segment.
* Sub-stage counts are the same for all genotypes. Genotype g has lifespan
  L_g = ℓ_g X_A (D-002): boundaries stay at the same absolute degree-days and
  only the segment in which L_g falls is shortened (or the last segment
  lengthened when ℓ_g > 1). Segments starting at or after L_g are never
  reached: the chain of that genotype ends there (death), and it never
  reaches the corresponding threshold (flagged by print()).
* Each adult sub-stage i of genotype g has its own width w_{g,i} (dd), so the
  transfer rate is q_{g,i}(t) = Δx(T(t)) / w_{g,i}.
* New parameter `X_ovi` (age at first oviposition, dd, default 61) separates
  the egg-laying threshold from X_pre, which keeps its role of scaling the
  degree-day increment (S1). The "gp2014" preset sets X_ovi = 7/25 × X_A =
  56.784 dd with n_A = 25, which gives 7 + 18 equal sub-stages, identical to
  the Gutierrez & Ponti chain (egg laying from class 8).

Reasoning: exact mean threshold ages for every genotype and every n_A, and an
exact version of D-002 (maturity fixed in degree-days, lifespan cost acting on
the last segment). The spread of threshold ages is still set by the number of
sub-stages before the threshold (CV = 1/sqrt(n)); about 25 sub-stages before
oviposition would match the observed 61 ± 12.2 dd.

## D-019 Reference solver design for M2 (Accepted, 2026-10-09)

* API: `simulate_ode(model, init, backend = "R")`. The pure-R reference
  backend comes first (M2); the C++ backend becomes the default in M3 with the
  same interface, so the two are compared directly in tests.
* Integration day by day with `deSolve`, holding each day's forcing constant,
  so midnight changes are exact breakpoints; releases are added between days,
  at the start of their release day (D-005).
* Output: daily states (end of each day) by default, with an `output_every`
  option; daily increments of cumulative counters (eggs laid, egg batches,
  emergences, matings, deaths by block, releases) for flow accounting.
* Regression target: new tight-tolerance legacy runs (rtol 1e-10, atol 1e-12,
  no floors) are added to `data-raw/reference/` as `*_tight.rds`; the existing
  reference files are kept.
* The 24-year legacy regression is a slow test, run only when the environment
  variable `SGDRIVE_SLOW_TESTS` is `"true"` (and in a manually triggered CI
  job); routine tests use the first 3 years.

## D-020 Initial conditions (Accepted, 2026-10-09)

`sg_init()` provides, in M2:
* `"constant"`: a count per sub-stage for each block (eggs, larvae, pupae,
  males, unmated females, mated females), split between genotypes by given
  weights (default: equally between the wild types). Adults go only to
  genotypes of the matching sex (males only where phi < 1, females only where
  phi > 0) and only to adult classes the genotype reaches (D-018). Mated
  females are split across mate genotypes in proportion to eta-weighted
  starting males (random mating).
* `"hardy_weinberg"`: genotype weights from allele frequencies for
  single-locus two-allele-letter genotype names (e.g. "AA", "Aa"), then as
  `"constant"`.
* `"state"`: a full state vector supplied by the user.
* (added by D-025) `"burn_in"` and `"equilibrium"`.

The defaults (0.25 per sub-stage for E, L, P; 0.125 for M and F; 0 for U)
reproduce the starting state of the Gutierrez & Ponti runs (0.25 per class,
adults split equally between the sexes, all females mated). Starting at
equilibrium or after a burn-in is deferred to M5 (calibration).

## D-021 C++ backend design for M3 (Accepted, 2026-10-11)

1. Integrator: our own Dormand-Prince 5(4) (Hairer's step-size control,
   componentwise error norm atol + rtol |y|), no Boost/odeint dependency.
2. Stiffness: Hairer's stiffness test flags a day; that day is then redone
   with deSolve's lsoda driven by the C++ rate kernel, recorded in the
   diagnostics (automatic per-day fallback rather than only a warning).
3. The whole day loop (releases, integration, recording) runs in C++; R and
   C++ meet once per simulation.
4. After M3 the default backend is "cpp"; "R" stays as the readable reference
   and test oracle.
5. Default tolerances are the same for both backends (rtol 1e-8,
   atol 1e-10).
6. The 24-year legacy regression becomes a routine test with the C++
   backend; with the R backend it stays slow-gated (D-019).
7. A manually triggered CI job runs the tests with sanitizers
   (AddressSanitizer, UBSan) for the C++ code.

The C++ code uses plain C++17 with Rcpp and std::vector, allocates no memory
inside the rate kernel, checks for user interrupts once per simulated day,
and reads the model through a flat "model view" built once in R.

Implemented in M3 (WP1-WP6). Outcome: no stiff day over 24 years of McAllen
weather (defaults and gp2014 preset); the fallback triggers only in
artificial settings (nu >= 1e4 per day). The C++ backend is 12-73 times
faster than the R backend at the same tolerances (bench/m3_cpp_vs_reference.R)
and needs 'deSolve' only when a fallback day occurs.

---

## D-022 Chain-binomial solver design for M4 (Accepted, 2026-10-12)

1. **Exact chain moves within a step.** For an individual in sub-stage i of
   a chain (hazard mu, the same in every sub-stage), the probabilities of
   being alive in sub-stage j >= i at the end of a step of length dt, of
   having left the chain alive, or of having died in it come from the
   matrix exponential of the chain's generator with advancement and death
   (computed once per day by uniformization; S16). An individual that
   leaves the chain during the step no longer dies at this chain's rate.
   Chain transit is then exact for any dt, unlike the classic
   one-sub-stage-per-step chain-binomial, whose sub-stage durations are too
   long by about (q + mu) dt / 2.
2. **Arrivals during a step** (eggs laid, individuals leaving the previous
   stage, emerged adults) enter at times spread over the step. They face
   the hazard of their new chain for the rest of the step and are placed at
   their expected position: the time average of the chain's transition
   matrix over the step (S16). Emerged females can also mate during the rest
   of the step (S18). (Amended 2026-10-12 in WP2: first drafted as "arrivals
   wait in sub-stage 1 until the end of the step", which the time-step study
   showed to be the dominant error, 3-5% at dt = 0.1 against below 0.6%
   now.) The remaining error comes from quantities held fixed within a step
   (wound search, male availability, larval crowding).
3. **Mating before ageing.** Unmated and mated females share their age chain
   and hazard, so mating can be drawn first (probability
   1 - exp(-nu_{g,i} dt), mate genotype from pi) and all females then age
   and die together; only male numbers and the receptivity flag are held at
   their start-of-step values. Large nu causes no stiffness. With
   `mating = "instant"` emerging females draw their mate at emergence.
   **Order of events and draws in a step** (from the state at the start of
   the step): aggregates; eggs (Poisson per female x mate genotype,
   multinomial split over offspring genotypes); mating; chain moves of E, L,
   P, M, U, F in state order (binomial death, then multinomial destination);
   emergence (binomial sex, then binomial survival, by genotype); arrivals in
   state order (binomial death, then multinomial placement; emerged females
   also binomial mating and multinomial mate genotype). Multinomials are
   drawn as sequential conditional binomials in category order. This
   canonical order lets the R reference step and the C++ step consume the
   same random numbers.
4. **Deterministic map.** Replacing every draw by its expectation gives a
   deterministic map whose value is the exact conditional mean of one
   stochastic step. Used in tests: replicate mean vs map (large sizes), map
   vs ODE (dt -> 0, first order).
5. **Shared kernel.** The C++ rate kernel is split into day rates and state
   aggregates used by both solvers (D-010); ODE results are unchanged
   (bitwise).
6. **Random numbers.** R's generator through Rmath (`set.seed()` works).
   Replicate r uses its own L'Ecuyer-CMRG stream derived from `seed`, so it
   is identical whether run alone, among others, sequentially or in
   parallel. Parallel runs use R processes (`future.apply`), not C++
   threads. The user's RNG state is restored; seeds are stored in `meta`.
7. **Integers.** Non-integer initial states and releases are rounded
   without bias in each replicate (floor plus a Bernoulli draw on the
   fraction). Counts are stored as doubles (exact integers to 2^53). Flow
   accounting is an exact integer identity. Egg laying: see D-024 (batches
   by default; the original Poisson single eggs remain an option).
8. **Output.** `simulate_cb()` returns an `sg_sim_cb` object: states
   [output day, state, replicate] or block x genotype totals (`keep`),
   daily counters per replicate, releases, extinction day per replicate
   (first day end with every compartment zero), the day the wild type is
   eliminated (first day end with no individual of the wild-type
   genotype(s); added in WP5, since released sterile insects keep the total
   above zero under SIT), and `meta`. A memory
   estimate is checked before running. A replicate stops early when extinct
   with no releases left.
9. **Time step.** dt divides the day (1/dt a whole number); default 0.1
   day. Confirmed by the time-step study (bench/m4_dt_convergence.R, six
   scenarios over a year of McAllen weather): at dt = 0.1 the deterministic
   map is within 0.13-0.55% of the ODE in every block (relative to the
   block's maximum), within 0.08% in total adults and 0.35% in eggs laid per
   year; at dt = 0.25 up to 1.9%, 0.27% and 0.9%.
10. **R reference step** (`.sg_cb_step_r`), pure R, used to test the C++ step
    draw for draw and the deterministic map to round-off.
11. **Population scale.** The defaults (1 km² node, H = 5 wounds per km²)
    give tens of adults at the seasonal peak and fewer than one individual
    in most McAllen winters, so stochastic runs at 1 km² go extinct almost
    every winter. Defaults are kept (ODE regressions); stochastic tests and
    examples use node areas of at least 1000 km²; realistic wound densities
    are calibrated in M5 (`calibrate_wound_density()`).

---

## D-023 Meaning of the cube's mating parameter eta (Accepted, 2026-10-12)

1. **eta is relative mating competitiveness** (MGDrivE v1 and MGDrivE2
   convention). It sets which male genotypes females mate with,
   pi_{g,m} = eta[g,m] M_bar_m / sum_h eta[g,h] M_bar_h, but not how fast
   females mate: a female mates at rate nu g_find whatever the eta of the
   males present, and only eta = 0 for every male present prevents mating.
   Same in the ODE and the chain-binomial solver.
2. **g_find counts all mature males equally**, not weighted by eta:
   g_find = M_bar / (h_M A + M_bar).

---

## D-024 Oviposition in batches in the chain-binomial solver (Accepted, 2026-10-12)

Screwworm females lay egg masses at wounds, not single eggs. In the
chain-binomial solver, egg laying is by default (`oviposition = "batches"`)
a two-level draw for each female x mate genotype pair in a step:
oviposition events ~ Poisson(Psi f W dt), then viable eggs ~ Poisson(events
x clutch x sum_k T) split over offspring genotypes. With `clutch_size = NA`
(default) a batch is R phi_ov eggs, the eggs a female lays in a day at the
day's temperature, and f = 1: Gutierrez & Ponti's Psi is the probability of
finding a wound per day, as in their egg-batch output. With `clutch_size`
set, f = R phi_ov / clutch_size. The mean egg production, and so the ODE and
the deterministic map, are unchanged; the variance of egg numbers grows from
the mean (Poisson single eggs, `oviposition = "eggs"`, the original D-022
choice, kept for comparison with MGDrivE-style models) to about
(1 + clutch) x mean. This matters in small populations (extinction
probability, time to eradication, genetic drift through sibships) and not
for large ones. The `egg_batches` counter becomes a count of events
(including those of sterile-mated females, whose eggs are not viable); in
the ODE its rate is Psi f W. Not modelled: the regular laying cycle of each
female, and aggregation of several females in one wound. Amends D-022 item 7.

---

## D-025 Initial states from burn-in and constant-temperature equilibrium (Accepted, 2026-10-12)

Extends D-020 (added in M4-WP4 at the user's request):
1. `sg_init(method = "burn_in", years, weather)`: run the ODE (C++
   backend, no releases) over the `years` before the start date, using the
   given weather or, without it, the model's first 365 days of forcing
   repeated, starting from a `"constant"` (or `"hardy_weinberg"`) state,
   and return the state at the start date. The result records the weather
   source and `start_effect`, the relative difference to the same burn-in
   from a start ten times larger (a warning above 1e-3). The ODE forgets
   its starting state within months at McAllen, so a one-year burn-in is
   the default.
2. `sg_init(method = "equilibrium", temperature)`: the equilibrium of a
   wild-type population at constant temperature, in closed form (model
   specification section 11), or by a one-dimensional root search with
   larval crowding or an Allee effect. Only for a single wild-type genotype
   that produces both sexes and breeds true; an empty state (with a message)
   when R0 <= 1. Useful for constant-temperature scenarios; at McAllen the
   seasonal burn-in is the appropriate start.

---

## D-026 Outputs, checks, calibration and plots for M5 (Accepted, 2026-10-12)

User decisions for milestone M5:
1. **Plots** use ggplot2 throughout, as a suggested package: `plot()`
   methods return ggplot objects (an error with an install hint when ggplot2
   is missing); the plotted data are also available as tables. Colours are
   colour-blind friendly (user request, WP4): Okabe-Ito by default,
   ColorBrewer Dark2 or Paired (9 to 12 series), viridis beyond.
2. **Case proxy**: wounds newly infested with viable eggs, per km^2 and
   period, with cases per host when the user gives a host density, and the
   daily fraction of available wounds newly infested. (Amended in WP2 at the
   user's request: a female mated to a sterile male lays eggs that cannot
   hatch and causes no myiasis case.) An oviposition event is a case when its
   batch holds at least one viable egg: probability 1 - exp(-C v_gm) for a
   batch of C eggs (clutch size, D-024) and cross viability v_gm = sum_k s ih
   tau. New counter `fertile_batches` (ODE rate Psi f sum W_gm (1 -
   exp(-C v_gm))); in the chain-binomial step the batches of each pair are
   split into fertile ones ~ Binomial(batches, 1 - exp(-C v)) whose eggs are
   zero-truncated Poisson sums, an exact refinement that leaves the egg
   distribution unchanged. Cases are counted at oviposition; viable eggs
   that later die of temperature are not removed. Planned replicate summary
   bands, extinction confidence intervals and time-to-eradication
   statistics were dropped (user decision): users analyse the raw
   replicate outputs.
3. **Calibration** (`calibrate_wound_density()`): by default the wound density H is set
   so that the mean adult density over a reference year equals a target
   (adults per km^2, a required argument). With h_M = 0 and no releases all
   equations are homogeneous of degree one in (state, H), so after a burn-in
   populations are exactly proportional to H (verified: 10.000000 times the
   adults at 10 times H) and the calibration is closed form; otherwise a
   one-dimensional root search.
4. **O-007 resolved**: the largest-rate diagnostic includes the mating rate
   nu g_find, and `check_model()` reports it separately.

Output definitions (WP1): long tables by stage (eggs, larvae, pupae, males,
unmated females, mated females), sex (juveniles "unsexed"), genotype and
mate genotype; allele frequencies from single-letter alleles in genotype
names by default (or a user allele table); female classes (all, unmated,
mated, wild type, egg-laying, and fertile-equivalent: mated females weighted
by the viability of their cross relative to wild type x wild type).
Simulations store the information these need (`meta$info`: area, H,
genotypes, wild type, laying classes, cross viabilities).

---

## D-027 Descriptive parameter and argument names (Accepted, 2026-10-14)

User-facing parameter names were the symbols of Gutierrez & Ponti (2014)
(`X_P`, `c_A`, `p_E0`, `a_s`, ...), Greek spellings (`nu`, `sigma0`,
`gamma_L`) or ad hoc words. At the user's request they are renamed to
descriptive snake_case names (rules in CONVENTIONS.md), in a clean rename
before the first release: the old names give an error naming the new one
(`.sg_param_renames()`), and the symbols remain in the `symbol` column of
`sg_param_table()`, in the model specification (which tabulates symbol and
name, section 13) and in internal code. Results are unchanged.

| Former name | New name |
|---|---|
| `T_L` | `dev_temp_min` |
| `T_U` | `dev_temp_decline` |
| `a_dev` | `dev_rate_slope` |
| `b_dev` | `dev_decline_steepness` |
| `dx_floor` | `dev_rate_floor` |
| `X_pre` | `dd_scale_preoviposition` |
| `T_opt` | `optimal_survival_temp` |
| `TU_ovip` | `oviposition_temp_max` |
| `egg_temperature` | `egg_temp_source` |
| `T_host` | `host_body_temp` |
| `X_E` | `egg_duration_dd` |
| `tau_L` | `larva_duration_days` |
| `X_P` | `pupa_duration_dd` |
| `X_A` | `adult_lifespan_dd` |
| `n_E` | `egg_substages` |
| `n_L` | `larva_substages` |
| `n_P` | `pupa_substages` |
| `n_A` | `adult_substages` |
| `X_ovi` | `oviposition_age_dd` |
| `X_rec` | `female_receptive_age_dd` |
| `X_mat` | `male_mature_age_dd` |
| `p_A0` | `adult_mortality_min` |
| `c_A` | `adult_mortality_curvature` |
| `p_E0` | `egg_failure_min` |
| `c_E` | `egg_failure_curvature` |
| `s_L` | `larva_survival_daily` |
| `mu_b_F` | `female_background_mortality` |
| `mu_b_M` | `male_background_mortality` |
| `larval_dd` | `larval_crowding_mode` |
| `gamma_L` | `larval_crowding_coef` |
| `mu_W` | `pupa_drought_mortality` |
| `rain_threshold` | `drought_free_rainfall` |
| `eps_s` | `survival_floor` |
| `R` | `eggs_per_day_max` |
| `H` | `wound_density` |
| `a_s` | `wound_search_rate` |
| `sigma0` | `wound_finding_max` |
| `search_denominator` | `wound_competitors` |
| `oviposition` | `oviposition_mode` |
| `nu` | `mating_rate` |
| `h_M` | `male_density_half_mating` |
| `mating` | `mating_mode` |
| `ext_threshold` | `extinction_threshold` |

Function arguments: `sg_cube(ell, mu_ref)` became `lifespan_fraction`,
`reference_mortality` (also in `sg_omega_to_ell()`); `p` became `params`
in `sg_forcing()`, `sg_degree_days()`, `sg_temp_mortality()`,
`sg_egg_viability()`, `sg_ovip_scalar()` and `sg_cold_index()`;
`sg_model(adult_substages)` became `adult_segment_substages`;
`calibrate_H()` became `calibrate_wound_density()`, its result field `H`
became `wound_density`; `simulate_cb()` metadata fields follow the new
names. Decisions above this one use the former names.

---

## D-028 Population size from a target number of adult females (Accepted, 2026-10-15)

For scenarios at a constant temperature the user wants a fixed population
size (e.g. 1000 adult females) without scaling by hand in the script.
`sg_model(adult_females = n)` sets the node area A = n / d, where d is the
equilibrium density of adult females (U + F, per km^2) of the wild type at
the forcing's constant temperature (D-025). The area is scaled rather than
the wound density H because (1) the densities in the parameters, H in
particular, are the calibrated biology (GP2014) and stay unchanged, so only
the size of the region (and with it the number of wounds, H A) changes; and
(2) every rate depends on counts only through densities (wounds, half-
saturation h_M and larval crowding per wound are all scaled by A), so the
equilibrium is exactly proportional to A for every model option, whereas
scaling H is exact only without an Allee effect. The forcing must have one
constant temperature (constant `tmean` and 5-day mean, and no drought stress
the equilibrium would miss); for weather-driven runs, `area` or
`calibrate_wound_density()` remain the tools. `sg_init(method =
"equilibrium")` defaults to that temperature. The target is kept in
`model$meta$adult_females` for `print()` and is dropped when the wound
density is changed by calibration.

---

## D-029 Allee threshold as an initial state (Accepted, 2026-10-15)

With an Allee effect (h_M > 0, mating "rate") the wild-type equilibrium
equation at constant temperature has two positive roots in the egg input
B: the stable upper root (D-025) and a lower, unstable root, the Allee
threshold. `sg_init(method = "threshold")` returns the lower root, found by
the same one-dimensional search (stepping down from a point between the
roots, in half-decade steps so a narrow positive range is not missed).
Without an Allee effect it returns an empty state with a message. At the
threshold almost all females are unmated (e.g. 1% mated at h_M = 20 males
km^-2, 20 degC), so the threshold is a property of a sparse population, not
of a number of adult females alone: a population of the same size taken
from the upper equilibrium (two thirds mated, with juveniles in the
pipeline) grows. The demonstration (scenario 03) therefore starts from
multiples of the threshold state. Because wounds do not limit at such low
density (Ψ = σ₀), the threshold is proportional to h_M. Equilibrium and
threshold states record `adult_females` (U + F) in their attribute.

---

## D-030 Host-wound-case model, stage 1 (Accepted, 2026-10-16)

The static model counts every egg batch with a viable egg as a case, with no
limit from the number of animals, and keeps a fixed number of wounds. At the
user's request an optional host-wound model (`wound_model = "hosts"`) adds
hosts without an attractive wound (HS), wounded (HW) and infested (HI, an
Erlang chain), coupled to the flies through the attractive sites
A_s = HW + alpha HI that replace H in wound search; a case is a new
infestation of a host, so at most N_H hosts are infested at once while hosts
can be infested again over time. Equations: MODEL_SPEC section 7a (H1-H3);
design and plan: dev/HOST_WOUND_DESIGN.md. Decisions (user, 2026-10-16):

1. Wounding set by `wounded_fraction` (wounded share without flies) and
   `wound_healing_days`; host_density x wounded_fraction takes the role of
   `wound_density`.
2. After an infestation the host returns to HS (wound healed).
3. Attractiveness of infested wounds alpha = 1 by default.
4. Infestation duration: Erlang chain with its own `infestation_days` (6)
   and `infestation_substages` (6), independent of the egg and larval
   chains.
5. `treatment_rate` (default 0) ends infestations and kills their larvae.
6. `"static"` stays the default; static results are unchanged.
7. `host_density` and `wounded_fraction` have no defaults; host mode
   requires them.
8. Names: blocks HS, HW, HI; counters hosts_wounded, wounds_healed,
   infestations, infestations_ended, infestations_treated,
   superinfestations; `sg_hosts()` columns unwounded, wounded, infested,
   hosts, prevalence.
9. Calibration: one `calibrate_population(parameter = )` (wound_density,
   wounded_fraction - the host default, since host density usually comes
   from a census - or host_density); `calibrate_wound_density()` kept.
10. Equilibria: a grid scan with stability classification for all
    nonlinear cases (`sg_equilibria()`); `sg_init()` uses the highest stable
    equilibrium with a message when there are several; `"threshold"` is the
    unstable one just below it.

Findings during implementation: (a) a wound receives at most
sigma0 a_s dx oviposition events per day; with the GP2014 a_s (0.001, set
for a relative H) that is about 0.01 per day at 25 degC, so with real wound
counts `wound_search_rate` needs calibrating (scenario 04 uses 0.02); (b)
at equilibrium egg balance fixes the strike rate per wound, so host
prevalence does not depend on alpha and alpha does not create alternative
stable states (it scales the flies); (c) two existing bugs fixed:
`sg_cases()` with `output_every > 1`, and HS/HW totals with
`keep = "totals"`. Later stages: superinfestation extending infestations,
seasonal wounding, several host classes (the `n_host` dimension),
myiasis mortality, detection and reporting of cases.

---

## Open questions

| ID | Question | Current assumption |
|---|---|---|
| O-006 (resolved by D-021) | Integrator for the C++ ODE backend (M3) | M2 diagnostics (1 year, McAllen): `lsoda` never switched to its stiff method; largest rate 13 per day with defaults and 40 per day under "gp2014", about 40-95 steps per day. 24-year gp2014 run (M2-WP5): stiff method on 2 of 8766 days, largest rate 43.3 per day, 93 steps per day. An explicit adaptive Runge-Kutta (Dormand-Prince) is expected to suffice; artificial stiffness arises only for extreme settings such as nu = 1000. Stress tests (M2-WP6): weekly sterile releases of 10x-1000x the wild males leave the largest rate unchanged with no stiff days; larval crowding with 1e6 released eggs raises it from 4.9 to 30.7 per day, still without stiff days. M3 should still include a stiffness check (Hairer's Dormand-Prince test) with a warning or fallback to a stiff method, and reuse tests/testthat/test-stress.R against the C++ solver |
| O-007 (resolved by D-026) | Should the largest-rate diagnostic (`max_rate`) include the effective mating rate ν·g_find? It currently counts transfer and mortality only, so stiffness caused by fast mating (found in M3-WP4: ν = 1e5, h_M = 1e5, large male release) is not visible in it | Unchanged (transfer + mortality), pending the user's decision |
| O-004 | Suspected `calcOmega` argument misuse in the legacy parameter table (0.4767 and 0.0740 look like reductions; `calcOmega` expects the fraction retained) | To be checked before building SIT cubes in M6 |
