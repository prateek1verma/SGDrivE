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

---

## Open questions

| ID | Question | Current assumption |
|---|---|---|
| O-004 | Suspected `calcOmega` argument misuse in the legacy parameter table (0.4767 and 0.0740 look like reductions; `calcOmega` expects the fraction retained) | To be checked before building SIT cubes in M6 |
