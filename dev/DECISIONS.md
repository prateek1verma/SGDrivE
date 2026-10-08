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

---

## Open questions

| ID | Question | Current assumption |
|---|---|---|
| O-001 | Where offspring viability τ applies | At egg production, as in MGDrivE2 |
| O-002 | μ_ref needed to convert cube ω to ℓ_g | User supplies μ_ref, or gives ℓ_g directly |
| O-003 | "both" sexes release: how the count is split | Equal split unless the user gives separate numbers |
| O-004 | Suspected `calcOmega` argument misuse in the legacy parameter table (0.4767 and 0.0740 look like reductions; `calcOmega` expects the fraction retained) | To be checked before converting legacy ω values |
| O-005 | Egg development temperature: ambient degree-days (Gutierrez code) or host temperature (biology) | Ambient degree-days, as in the Gutierrez code |
