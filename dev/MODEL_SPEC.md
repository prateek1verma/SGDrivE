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

No default exceeds 10. Setting every n to 25 (with X_ovi = 56.784 dd) gives the
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

    q_E(t)   = n_E Δx(T) / X_E            X_E = 10 dd       [GP2014 code; D-014]
               (egg_temperature = "host": q_E = n_E Δx(T_host) / X_E, constant)
    q_L      = n_L / τ_L                  τ_L = 6 d         [GP2014 code]
    q_P(t)   = n_P Δx(T) / X_P            X_P = 95 dd       [GP2014]
    q_{g,i}(t) = Δx(T) / w^F_{g,i}  (females),  Δx(T) / w^M_{g,i}  (males)
                                          per adult sub-stage, see S6  [adapted, D-018]

**(S6) Segmented adult chains and maturity thresholds** [adapted, D-002, D-004, D-018]

Maturity thresholds are fixed physiological ages (dd), the same for every
genotype:

| Threshold | Symbol | Default | Applies to |
|---|---|---|---|
| First oviposition | X_ovi | 61 dd | egg laying and wound search (F) |
| Female receptivity | X_rec | 0 dd | mating of U (dropped when mating = "instant") |
| Male maturity | X_mat | 0 dd | male availability as mates |

Each adult chain is split at its thresholds into segments: females at
{X_rec, X_ovi}, males at {X_mat}; zero-length segments are dropped. With
wild-type boundaries 0 = b_0 < b_1 < … < b_K = X_A, segment j covers
[b_{j−1}, b_j) and has n_j sub-stages, with Σ_j n_j = n_A. The n_j are
allocated in proportion to the segment lengths, at least one each (largest
remainder), unless given by `sg_model(adult_substages = )`.

Genotype g has lifespan L_g = ℓ_g X_A (D-002). Its segment j covers
[b_{j−1}, min(b_j, L_g)) for j < K and [b_{K−1}, L_g) for the last segment;
each of its n_j sub-stages has width

    w_{g,i} = (segment end − segment start) / n_j      (dd)

Segments that start at or after L_g are never reached: the chain of
genotype g ends with the sub-stage in which L_g falls (its last class
l_g), and outflow from that class is death. Sub-stage i is *receptive*,
*laying* or *mature* when its segment starts at or after X_rec, X_ovi or
X_mat; these flags are the same for all genotypes, and a genotype whose last
class l_g comes before the first laying (receptive, mature) class never lays
(mates, counts as a mate).

Default females (n_A = 10, X_rec = 0): segments [0, 61) and [61, 202.8) with
3 and 7 sub-stages, so the mean age at the first laying class is exactly
61 dd. Preset "gp2014": X_ovi = 56.784 dd and n_A = 25 give 7 + 18 sub-stages
of equal width 8.112 dd, identical to the Gutierrez & Ponti chain (egg laying
from class 8).

## 6. Mortality hazards

**(S7)**

    μ_E(t)          from (S3)
    μ_L(t)          = −ln(s_L) [+ optional density dependence]      s_L = 0.95 d⁻¹ survival [GP2014]
    μ_P(t)          = μ_T(T) [+ optional moisture term]
    μ_M(t) = μ_U(t) = μ_F(t) = μ_T(T) + μ_b

μ_b = 0 by default (optional sex-specific background field mortality).
Moisture index I_W(t) = min(1, annual rainfall / rain_threshold) for the
calendar year of day t, with incomplete years scaled to a full year; I_W = 1
when no rainfall is supplied (μ_W = 0 by default, so no effect).
Genotype does **not** change these hazards (D-002); genotype lifespan effects
act through the sub-stage widths w_{g,i} only (S6).

## 7. Reproduction

**(S8) Egg-laying females and wound search** [GP2014 Eq. 4] + [adapted, D-006]

    W_{g,m} = Σ_{i laying} F_{g,m,i}                 (egg-laying mated females, S6)
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

    M̄_m     = Σ_{i mature} M_{m,i}                 (mature males of genotype m, S6)
    M̄       = Σ_m M̄_m
    g_find  = M̄ / (h_M A + M̄)                    (= 1 if h_M = 0 and M̄ > 0; 0 if M̄ = 0)
    π_{g,m} = η[g,m] M̄_m / Σ_h η[g,h] M̄_h        (no mating if the denominator is 0)
    ν_{g,i} = ν · 1[i receptive] · g_find

h_M is a density (males km⁻²) converted to counts by the node area A.

**(S13) Adults** (for each genotype g, sub-stage i; X_{·,0} terms are the inflows)

    dM_{g,i}/dt   = q^M_{g,i−1} M_{g,i−1} − (q^M_{g,i} + μ_M) M_{g,i}
    dU_{g,i}/dt   = q^F_{g,i−1} U_{g,i−1} − (q^F_{g,i} + μ_U + ν_{g,i}) U_{g,i}
    dF_{g,m,i}/dt = q^F_{g,i−1} F_{g,m,i−1} − (q^F_{g,i} + μ_F) F_{g,m,i} + ν_{g,i} π_{g,m} U_{g,i}

with q^M_{g,0} M_{g,0} ≡ (1 − φ_g) ξ_M,g ε_g, q^F_{g,0} U_{g,0} ≡ φ_g ξ_F,g ε_g and
F_{g,m,0} ≡ 0. Outflow from the last reached class l_g (S6) is death by old age;
sub-stages after l_g stay empty.
A mating female keeps her age class i.

**(S14) Instant-mating switch** `mating = "instant"` [MGDrivE2]: emerging
females go directly to F_{g,m,1} with weights π_{g,m}, or to U_{g,1} only when
M̄ = 0; U then mates at rate ν once males appear. Equivalent to the limit
ν → ∞ with X_rec = 0, without stiffness.

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

**ODE** (`simulate_ode()`, D-019, D-021). The model is integrated one day at
a time, from t = d − 1 to t = d, with that day's forcing held constant, so the
change of weather at midnight is an exact breakpoint. Releases are added to
the state at the start of their day, before it is integrated. The counters
(S-flows: eggs laid, matings, emergences, deaths) are integrated alongside the
states and reset each day, giving exact daily totals.

* **C++ backend** (default). The day loop runs in C++ with the package's own
  Dormand–Prince 5(4) integrator (Hairer's DOPRI5: PI step-size control with
  β = 0.04, starting-step estimate, error norm
  sqrt(mean((err_j / (atol + rtol·max(|y_j|, |y_new_j|)))²)), first same as
  last). Hairer's stiffness test estimates h|λ| from the last two stages at
  every accepted step; if it exceeds 3.25 on 15 steps of a day, that day is
  integrated again from its start with `lsoda` driven by the C++ rates, and
  the C++ loop resumes on the next day. Default tolerances rtol = 1e-8,
  atol = 1e-10.
* **R backend.** The same equations in pure R (the readable reference and
  test oracle), solved day by day with any `deSolve` method (default
  `lsoda`). Any `deSolve` method can also be driven by the C++ rates.

The two backends agree to within the integration tolerance (tested to 1e-8 at
rtol = 1e-10 over many scenarios), and both reproduce the legacy Gutierrez &
Ponti (2014) ODE to 1e-6. Over 24 years of McAllen weather the model was
never stiff with the defaults or the gp2014 preset; artificial stiffness
needs extreme settings such as near-instant mating (ν ≥ 10⁴ d⁻¹).

**Chain-binomial** (`simulate_cb()`, D-022). Counts are whole numbers.
Each day is divided into steps of length dt (1/dt a whole number, default
0.1 day); every step of day d uses day d's forcing. Releases are added at the
start of their day, as in the ODE.

**(S15) Aggregates.** At the start of each step, W_{g,m}, N_s, Ψ (S8), M̄_m,
g_find, π_{g,m} (S12) and the larval hazard μ_L (S7) are computed from the
current counts and held fixed for the step.

**(S16) Chain moves.** For a chain with sub-stage rates q_1, …, q_n, last
reached sub-stage l and hazard μ (the same in every sub-stage), let A(t) =
exp(G t), where G is the generator of advancement (sub-stage i → i + 1 at
rate q_i for i < l, and from l out of the chain, an absorbing state) and
death (rate μ in every sub-stage). A_{i,j}(dt) is the probability of being
alive in sub-stage j at the end of a step that started in sub-stage i,
A_{i,out}(dt) of having left the chain alive (death stops acting at exit),
and δ_i = 1 − Σ_j A_{i,j} − A_{i,out} of having died in it. A is computed
once per day for each chain (E, P per day; adults per genotype and sex;
larvae in every step, as μ_L includes crowding from S15). From sub-stage i
with X_i individuals:

    deaths   D_i ~ Binomial(X_i, δ_i)
    movers   (N_{i→i}, …, N_{i→l}, N_{i→out}) ~ Multinomial(X_i − D_i; A_{i,·} / (1 − δ_i))

Leaving the chain means: E → L, L → P, P → emergence (S19), adults → death
by old age. Transit through a chain is exact for any dt.

**Arrivals.** Individuals entering a chain during a step (eggs laid, those
leaving the previous stage, emerged adults) enter at times spread uniformly
over the step. With Ā = (1/dt) ∫_0^dt A(u) du, an arrival dies before the end
of the step with probability δ̄ = 1 − Σ_j Ā_{1,j} − Ā_{1,out}, and survivors
are placed over sub-stages by Multinomial(Ā_{1,·} / (1 − δ̄)), the small
probability Ā_{1,out} of leaving again being kept in sub-stage l. Arrivals of
a step are added at its end.

**(S17) Eggs.** For each female × mate genotype pair (g, m), with
T[g,m,k] = s_g · ih[g,m,k] · τ[g,m,k] (S9). With `oviposition = "batches"`
(default, D-024):

    B_{g,m}         ~ Poisson(Ψ f W_{g,m} dt)                     oviposition events
    N_{g,m}         ~ Poisson(B_{g,m} · C · Σ_k T[g,m,k])          viable eggs
    (N_{g,m,k})_k   ~ Multinomial(N_{g,m}; T[g,m,·] / Σ_k T[g,m,k])

with clutch size C = R φ_ov and f = 1 when `clutch_size` is NA (a batch is
the eggs of one day; Ψ is the daily probability of finding a wound), or
C = `clutch_size` and f = R φ_ov / C. With `oviposition = "eggs"`:

    N_{g,m}         ~ Poisson(R φ_ov Ψ W_{g,m} Σ_k T[g,m,k] dt)

then the same split. Both have mean R φ_ov Ψ W_{g,m} Σ_k T dt (S9); batches
raise the variance to about (1 + C) times the mean. The eggs arrive in the
egg chain of genotype k. The egg_batches counter counts events B (ODE rate
Ψ f Σ W).

**(S18) Mating** (before ageing; unmated and mated females share age chain
and hazard, so the order does not bias ageing or survival):

    Y_{g,i}             ~ Binomial(U_{g,i}, 1 − e^{−ν_{g,i} dt})      ν_{g,i} = ν · 1[i receptive] · g_find
    (Y_{g,m,i})_m       ~ Multinomial(Y_{g,i}; π_{g,·})

The mated females move to F_{g,m,i} and then age with all other females
(S16). ν_{g,i} = 0 when no compatible male is present. Emerged females that
survive their arrival mate during the rest of the step with probability
1 − (1 − e^{−ν_{g,1} dt}) / (ν_{g,1} dt) (a uniform arrival time), mates from
π_{g,·}. With `mating = "instant"` the receptivity flag is not used, and
surviving emerged females are split over mates by Multinomial(π_{g,·}) (to U
only when no compatible male is present).

**(S19) Emergence.** For the Z_g pupae of genotype g leaving the pupal chain:

    females Z^F_g ~ Binomial(Z_g, φ_g)
    surviving females: Binomial(Z^F_g, ξ_F,g)      surviving males: Binomial(Z_g − Z^F_g, ξ_M,g)

the rest die at emergence. Survivors arrive in the adult chains (S16, S18).

**Order of draws in a step.** Aggregates (S15); eggs (S17, pairs (g, m) with
m fastest: events, then eggs, then split); mating (S18, genotype then sub-stage); chain moves (S16) block by
block (E, L, P, M, U, F) and compartment by compartment in state order;
emergence (S19, by genotype); arrivals block by block in state order (eggs
by genotype, larvae, pupae, males; then females by genotype: death, mating,
mates, placement of unmated, placement of mated by mate genotype). All draws
use the counts at the start of the step except that mating moves females
before the chain moves. Multinomials are sequential conditional binomials in
category order. The R reference step and the C++ step follow this order, so
they use the same random numbers.

**Deterministic map.** Replacing every draw by its expectation gives a
deterministic map, the exact conditional mean of one step. As dt → 0 it
converges to the ODE, with an error from the S15 quantities held fixed
within a step; at the default dt = 0.1 it is within about 0.5% of the ODE
in every block (bench/m4_dt_convergence.R).

**Integers.** Non-integer initial states and releases are rounded without
bias in each replicate: x → ⌊x⌋ + Bernoulli(x − ⌊x⌋). Counts never become
negative. Random numbers come from R's generator; replicate r uses
L'Ecuyer-CMRG stream r derived from the seed (D-022).

**Extinction.** No population floors. A stochastic replicate is extinct on
the first day that ends with every compartment at zero. Eradication
summaries with an extinction threshold (default 1 individual) for the ODE
follow in M5.

## 11. Initial conditions

`sg_init()` (D-020, D-025):

* `"constant"`, `"hardy_weinberg"`, `"state"`: given counts per sub-stage
  (the default reproduces the Gutierrez & Ponti starting state), genotype
  weights from allele frequencies, or a full state vector.
* `"burn_in"`: the ODE without releases is run for a number of years before
  the start date, over the weather of those years or, without it, over the
  model's first 365 days of forcing repeated, and the state at the start
  date is returned. With the defaults at McAllen the ODE forgets its
  starting state within about six months (relative differences below 1e-6
  after 180 days, 1e-14 after a year); the reported `start_effect` is the
  difference to the same burn-in from a start ten times larger.
* `"equilibrium"`: the equilibrium of a wild-type population (one wild
  genotype that breeds true) at constant temperature T. With the day's
  quantities at T (S1-S7; no moisture term), every chain is linear, so for
  an egg input B all compartments are B times a fixed profile: Erlang
  chains E_i = B/(q+μ) (q/(q+μ))^{i−1}, segmented adult chains
  X_i = X_{i−1} q_{i−1}/(q_i + μ + ν_i), unmated females with the mating
  rate ν_i (S12) and mated females as the difference to all females. With
  W(1) and N_s(1) the egg-laying and wound-searching females per unit egg
  input and T_ww = s_w ih[w,w,w] τ[w,w,w], equilibrium requires
  R φ_ov T_ww Ψ* W(1) = 1, so

      Ψ*  = 1 / (R φ_ov T_ww W(1))            (needs R0 = σ₀ R φ_ov T_ww W(1) > 1)
      N_s* = a_s Δx H / (−ln(1 − Ψ*/σ₀))       (S8 solved for N_s)
      B*  = N_s* / N_s(1)

  R0 is the number of eggs per egg at low density; without a positive
  equilibrium (R0 ≤ 1 or no development) the state is empty. With larval
  crowding (μ_L depends on L, solved by a one-dimensional root) or an Allee
  effect (g_find depends on M̄) B* is found by a one-dimensional root
  search, taking the upper (stable) root.

## 12. Parameters

The authoritative list of parameters, with names, defaults, allowed ranges,
units, sources and the equations that use them, is the parameter table in
`R/param-table.R`, returned by `sg_param_table()` (D-017). `print(sg_params())`
shows the current values.

Groups: temperature (T_L, T_U, a_dev, b_dev, dx_floor, T_opt, TU_ovip,
egg_temperature, T_host); development (X_pre, X_E, tau_L, X_P, X_A); stages
(n_E = 3, n_L = 6, n_P = 10, n_A = 10); thresholds (X_rec, X_mat and the class
X_ovi; thresholds define the segmented adult chains, S6); mortality (c_A, p_A0, c_E, p_E0, s_L, mu_b_F, mu_b_M,
larval_dd, gamma_L, mu_W, rain_threshold, eps_s); reproduction (R, a_s, H,
sigma0, search_denominator); mating (nu, h_M, mating); numerical
(ext_threshold).

Cube-derived quantities (ℓ_g, φ_g, ξ_F,g, ξ_M,g, s_g, η, ih, τ) come from
`sg_cube()`; the node area A comes from `sg_model()`.

Preset `"gp2014"` (regression configuration of the Gutierrez & Ponti code):
n_E = n_L = n_P = n_A = 25, X_ovi = 56.784 dd, dx_floor = 0.01,
search_denominator = "all_adults", mating = "instant".

Values marked "calibrate" in the table (notably H) should be set so the model
matches target densities (e.g. 7–54 females km⁻² in the field, Gutierrez et
al. 2019).
