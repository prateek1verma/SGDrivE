# SGDrivE v0.1 model specification

This document is the agreed specification that milestones M1–M3 implement.
Equation labels (S1, S2, …) are cited in code comments. Changes to the model
require a change here and, where relevant, a new entry in `DECISIONS.md`.

Source tags: **[GP2014]** Gutierrez & Ponti (2014) paper or Pascal code;
**[G2019]** Gutierrez, Ponti & Arias (2019); **[MGDrivE]** MGDrivE v1 cubes;
**[MGDrivE2]** MGDrivE2 lifecycle model; **[adapted]** new for SGDrivE.

---

## 1. Scope of v0.1

* One node (area A km²); one host class. Array dimensions for node and host
  class exist but have length 1 (D-009, CONVENTIONS).
* Any MGDrivE v1 inheritance cube, used unchanged.
* Deterministic ODE and stochastic chain-binomial solvers sharing one C++ rate
  kernel (D-010).
* Daily weather input; temperature-dependent terms precomputed per day and
  held constant within each day.
* Not in v0.1: movement, host/wound dynamics, larval transport, explicit
  gonotrophic cycle, male mating capacity.

## 2. Indices and states

Genotypes g, m, k ∈ G (|G| = n_G): g = own genotype (female genotype in F),
m = genotype of the male a female mated, k = offspring genotype. Sub-stages
i = 1…n_X for stage X.

| State | Meaning | Size | Source |
|---|---|---|---|
| E_{g,i} | eggs | n_G · n_E | [GP2014] |
| L_{g,i} | larvae | n_G · n_L | [GP2014] |
| P_{g,i} | pupae | n_G · n_P | [GP2014] |
| M_{g,i} | adult males, by physiological age | n_G · n_A | [GP2014] + [MGDrivE2] |
| U_{g,i} | unmated (virgin) adult females, by age | n_G · n_A | [MGDrivE2], [adapted] |
| F_{g,m,i} | mated females of genotype g carrying sperm of genotype m, by age | n_G² · n_A | [MGDrivE2], [adapted] |

All states are **counts per node** (D-009). Cumulative counters (eggs laid,
egg batches, emergences, matings, deaths by cause, releases) are tracked for
flow accounting and outputs.

## 3. Default sub-stage numbers (D-003)

| Stage | Default n | Mean duration | Reason |
|---|---|---|---|
| Eggs | 3 | 10 dd (< 1 day) | short stage |
| Larvae | 6 | 6 days | medium |
| Pupae | 10 | 95 dd (≈ 8 days at 27 °C) | long |
| Adults | 10 | 202.8 dd (≈ 17 days at 27 °C) | long; carries maturity thresholds |

No default exceeds 10. Setting every n to 25 (with i0 = 8) gives the
Gutierrez & Ponti configuration used in regression tests.

## 4. Temperature-driven functions (precomputed per day)

Daily mean temperature T(t) (°C); T̄₅(t) = mean of the last 5 days (fewer at
the start of a run).

**(S1) Physiological time** [GP2014 Eq. 1, code]

    r(T)  = a (T − T_L) / (1 + b^(T − T_U))            [h⁻¹]
    Δx(T) = X_pre · max(dx_floor, 24 · r(T))           [dd d⁻¹]

with a = 0.00065, b = 1.75, T_L = 14.5 °C, T_U = 40.5 °C (code; 40 in paper),
X_pre = 61 dd. `dx_floor` = 0 by default; 0.01 reproduces the Pascal floor of
0.61 dd/d.

**(S2) Pupal/adult temperature mortality** [GP2014 Eq. 2]

    p_A(T) = c_A (T − T_opt)² + p_A0,      μ_T(T) = −ln(max(ε_s, 1 − p_A(T)))

c_A = 0.00036, p_A0 = 0.0035, T_opt = 27.2 °C, ε_s = 1e−5.

**(S3) Egg viability** [GP2014 Eq. 3]

    h(T̄₅) = 1 − min(1, c_E (T̄₅ − T_opt)² + p_E0),     μ_E(t) = −ln(max(ε_s, h(T̄₅)))

c_E = 0.00375, p_E0 = 0.3.

**(S4) Oviposition temperature scalar** [GP2014 Eq. 6]

    φ_ov(T) = max(0, 1 − ((T − T_L − ψ)/ψ)²),    ψ = (43.5 − 14.5)/2 = 14.5

## 5. Development rates and maturity thresholds

**(S5) Transfer rates** (per sub-stage, d⁻¹)

    q_E(t)   = n_E Δx(T) / X_E            X_E = 10 dd       [GP2014 code; O-005]
    q_L      = n_L / τ_L                  τ_L = 6 d         [GP2014 code]
    q_P(t)   = n_P Δx(T) / X_P            X_P = 95 dd       [GP2014]
    q_A,g(t) = n_A Δx(T) / (ℓ_g X_A)      X_A = 202.8 dd    [GP2014] + [adapted, D-002]

**(S6) Maturity thresholds in absolute degree-days** [adapted, D-002, D-004]

For a threshold of X_thr dd and genotype g, the first adult class at or past
the threshold is

    i_thr,g = round(n_A · X_thr / (ℓ_g · X_A)) + 1

If i_thr,g > n_A the threshold is never reached (that genotype never matures).
Users may instead set class indices directly.

| Threshold | Symbol | Default | Applies to |
|---|---|---|---|
| Pre-oviposition | X_pre → i0_g | 61 dd | egg laying (F) |
| Female receptivity | X_rec → i_m,g | 0 dd (class 1) | mating of U |
| Male maturity | X_mat → i_mM,g | 0 dd (class 1) | male availability |

With n_A = 10 and ℓ_g = 1, i0 = round(3.008) + 1 = 4 (laying after ≈ 61 dd).
Gutierrez & Ponti used classes 8–25 of 25 (i0 = 8); the regression preset
sets i0 = 8 directly.

## 6. Mortality hazards

**(S7)**

    μ_E(t)          from (S3)
    μ_L(t)          = −ln(s_L) [+ optional density dependence]      s_L = 0.95 d⁻¹ survival [GP2014]
    μ_P(t)          = μ_T(T) [+ optional moisture term]
    μ_M(t) = μ_U(t) = μ_F(t) = μ_T(T) + μ_b

μ_b = 0 by default (optional sex-specific background field mortality).
Genotype does **not** change these hazards (D-002); genotype lifespan effects
act through q_A,g only.

## 7. Reproduction

**(S8) Egg-laying females and wound search** [GP2014 Eq. 4] + [adapted, D-006]

    W_{g,m} = Σ_{i ≥ i0_g} F_{g,m,i}                (egg-laying mated females)
    N_s     = Σ_{g,m} W_{g,m}                        (default: searching mated females, incl. sterile-mated)
    N_s     = Σ M + Σ U + Σ F                        (option "all_adults", [GP2014 code])
    Ψ(t)    = σ₀ [1 − exp(−a_s Δx(T) H / N_s)]       (Ψ = σ₀ when N_s = 0)

H = number of wounds available in the node (= wound density × A), a_s =
search constant (dd⁻¹), σ₀ = baseline search efficiency (default 1).

**(S9) Egg production by offspring genotype** [MGDrivE] + [GP2014 Eq. 5]

    B_k(t) = R φ_ov(T) Ψ(t) Σ_g Σ_m s_g · ih[g,m,k] · τ[g,m,k] · W_{g,m}

R = maximum eggs per egg-laying female per day (90 [GP2014]). τ applied at
egg production (O-001). Implementation: precompute T[g,m,k] = s_g · ih · τ,
then B = R φ_ov Ψ · (vec(W)ᵀ · matrix(T, n_G², n_G)).

## 8. Life-cycle equations

Chains: for a stage with rate q and hazard μ,
dX_1/dt = inflow − (q + μ) X_1, dX_i/dt = q X_{i−1} − (q + μ) X_i (i ≥ 2).

**(S10) Eggs, larvae, pupae** (for each genotype k)

    dE_{k,1}/dt = B_k − (q_E + μ_E) E_{k,1}
    dL_{k,1}/dt = q_E E_{k,n_E} − (q_L + μ_L) L_{k,1}
    dP_{k,1}/dt = q_L L_{k,n_L} − (q_P + μ_P) P_{k,1}
    ε_k         = q_P P_{k,n_P}                         (emergence flux)

**(S11) Emergence** [MGDrivE2]

    into M_{k,1}: (1 − φ_k) ξ_M,k ε_k        into U_{k,1}: φ_k ξ_F,k ε_k
    the remainder (1 − ξ) dies at emergence.

**(S12) Male availability and mating** [MGDrivE2] + [adapted, D-004]

    M̄_m     = Σ_{i ≥ i_mM,m} M_{m,i}              (mature males of genotype m)
    M̄       = Σ_m M̄_m
    g_find  = M̄ / (h_M A + M̄)                    (= 1 if h_M = 0 and M̄ > 0; 0 if M̄ = 0)
    π_{g,m} = η[g,m] M̄_m / Σ_h η[g,h] M̄_h        (no mating if the denominator is 0)
    ν_{g,i} = ν · 1[i ≥ i_m,g] · g_find

h_M is a density (males km⁻²) converted to counts by the node area A.

**(S13) Adults** (for each genotype g, sub-stage i; X_{·,0} terms are the inflows)

    dM_{g,i}/dt   = q_A,g (M_{g,i−1} − M_{g,i}) − μ_M M_{g,i}
    dU_{g,i}/dt   = q_A,g (U_{g,i−1} − U_{g,i}) − (μ_U + ν_{g,i}) U_{g,i}
    dF_{g,m,i}/dt = q_A,g (F_{g,m,i−1} − F_{g,m,i}) + ν_{g,i} π_{g,m} U_{g,i} − μ_F F_{g,m,i}

with M_{g,0} q_A,g ≡ (1 − φ_g) ξ_M,g ε_g, U_{g,0} q_A,g ≡ φ_g ξ_F,g ε_g and
F_{g,m,0} ≡ 0. Outflow q_A,g X_{n_A} from the last class is death by old age.
A mating female keeps her age class i.

**(S14) Instant-mating switch** `mating = "instant"` [MGDrivE2]: emerging
females go directly to F_{g,m,1} with weights π_{g,m}, or to U_{g,1} only when
M̄ = 0; U then mates at rate ν once males appear. Equivalent to the limit
ν → ∞ with i_m = 1, without stiffness.

## 9. Releases (D-005)

Release events at day boundaries add individuals of a release genotype r:

| Stage | Added to |
|---|---|
| eggs | E_{r,1} |
| pupae | P_{r,1} |
| adults, sex "male" | M_{r,1} |
| adults, sex "female" | U_{r,1} (unmated) |
| adults, sex "both" | M_{r,1} and U_{r,1} (split per O-003) |

Release sizes may be given per km² and are converted with A.

## 10. Solvers

**ODE.** Dormand–Prince RK45 with adaptive steps, integrated one day at a
time so the daily forcing is constant within each call; releases applied
between days. A pure-R `deSolve` path calling the same rate function is kept
for cross-checks.

**Chain-binomial (time step dt, rates fixed within the step).**

| Process | Draw |
|---|---|
| Chain compartment (advance or die) | leaving ~ Binomial(N, 1 − e^{−(q+μ)dt}); advancing ~ Binomial(leaving, q/(q+μ)) |
| Virgins (advance, die, mate) | Multinomial over the three outcomes with total probability 1 − e^{−(q+μ+ν)dt}; mate genotype ~ Multinomial(π_{g,·}) |
| Eggs from (g, m) | viable eggs ~ Poisson(R φ_ov Ψ s_g W_{g,m} Σ_k ih τ · dt); split ~ Multinomial(ih τ normalised) |
| Emergence | survival ~ Binomial(ξ); sex ~ Binomial(φ) |

Counts never become negative. A warning is issued when max (q + μ) dt > 0.1.
Random numbers come from R's generator, so `set.seed()` reproduces runs.

**Extinction.** No population floors. Eradication is reached when all wild
stages are below an extinction threshold (default 1 individual).

## 11. Initial conditions

`sg_init()` options: constant counts per sub-stage, user-supplied state, or a
burn-in run under repeated weather.

## 12. Parameters (defaults)

| Symbol | R name (planned) | Description | Units | Default | Source |
|---|---|---|---|---|---|
| T_L | `T_L` | lower developmental threshold | °C | 14.5 | [GP2014] |
| T_U | `T_U` | upper parameter of development curve | °C | 40.5 | [GP2014] code |
| a | `a_dev` | slope of r(T) | h⁻¹ °C⁻¹ | 0.00065 | [GP2014] |
| b | `b_dev` | high-temperature decline | – | 1.75 | [GP2014] |
| X_pre | `X_pre` | pre-oviposition period; scales Δx | dd | 61 | [GP2014] |
| dx_floor | `dx_floor` | floor on 24·r(T) | – | 0 (0.01 = GP) | [GP2014] code |
| X_E | `X_E` | mean egg duration | dd | 10 | [GP2014] code |
| τ_L | `tau_L` | mean larval duration | d | 6 | [GP2014] code |
| X_P | `X_P` | mean pupal duration | dd | 95 | [GP2014] |
| X_A | `X_A` | mean adult physiological lifespan | dd | 202.8 | [GP2014] |
| n_E, n_L, n_P, n_A | `n_E` … | Erlang sub-stages | – | 3, 6, 10, 10 | [adapted] D-003 |
| T_opt | `T_opt` | temperature of minimum mortality | °C | 27.2 | [GP2014] |
| c_A, p_A0 | `c_A`, `p_A0` | pupal/adult mortality curve | °C⁻²; – | 0.00036; 0.0035 | [GP2014] |
| c_E, p_E0 | `c_E`, `p_E0` | egg hatch-failure curve | °C⁻²; – | 0.00375; 0.3 | [GP2014] |
| s_L | `s_L` | daily larval survival | – | 0.95 | [GP2014] |
| μ_b | `mu_b_F`, `mu_b_M` | background field mortality | d⁻¹ | 0 | [adapted] |
| R | `R` | max eggs per egg-laying female per day | eggs d⁻¹ | 90 | [GP2014] |
| T_ov, ψ | – | oviposition optimum, half-range | °C | 29, 14.5 | [GP2014] |
| a_s | `a_s` | search constant | dd⁻¹ | 0.001 | [GP2014] |
| H | `H` | wounds available in the node | count | calibrated (`calibrate_H()`) | [GP2014] + D-009 |
| σ₀ | `sigma0` | baseline search efficiency | – | 1 | [adapted] |
| search denominator | `search_denominator` | competitors for wounds | – | `"searching_females"` | D-006 |
| ν | `nu` | mating rate of receptive virgins | d⁻¹ | ln 2 | [G2019] |
| h_M | `h_M` | mate-finding half-saturation | males km⁻² | 0 | [adapted] |
| X_rec, X_mat | `X_rec`, `X_mat` | female receptivity, male maturity | dd | 0, 0 | [adapted] D-004 |
| ℓ_g | from cube | genotype lifespan fraction | – | 1 | [MGDrivE] → D-002 |
| φ_g, ξ_F,g, ξ_M,g, s_g, η, ih, τ | from cube | sex ratio, emergence, fecundity, mating fitness, inheritance, viability | – | cube | [MGDrivE] |
| A | `area` | node area | km² | 1 | D-009 |
| ε_s | – | survival clamp | – | 1e−5 | [GP2014] code |

Values marked "calibrated" are set so the model matches target densities
(e.g. 7–54 females km⁻² in the field, Gutierrez et al. 2019).
