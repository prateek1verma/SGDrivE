# Host–wound–case model, stage 1 (design note)

Status: **draft**, 2026-10-15; decisions of section 9 made by the user on
2026-10-16. Nothing here is implemented.
Proposed decision number: D-030.

## 1. Why

In the current model (MODEL_SPEC S7–S9, S17, D-026):

* wounds are a fixed count H = `wound_density` × area, never depleted,
  created or healed;
* there are no hosts; `sg_cases(host_density = )` only divides by a number
  given by the user;
* a "case" is every oviposition event with at least one viable egg
  (`fertile_batches`), so ten batches on one animal are ten cases, and the
  number of cases depends on how eggs are grouped into batches
  (`clutch_size`).

So nothing limits cases by the number of animals: with 2 animals the model
can report many simultaneous cases. Stage 1 adds a minimal host population
in which **at most N_H animals can be infested at any time**, while the same
animal can be infested again and again over a period.

## 2. Scope of stage 1

In stage 1:

* one host class (the `n_host` dimension of the state layout stays 1);
* host states S (no attractive wound), W (attractive, uninfested wound),
  I (infested, an Erlang chain I_1 … I_n), counts per node;
* constant wounding and healing rates (no seasonality);
* Erlang infestation duration with its own mean (`infestation_days`) and
  number of sub-stages (`infestation_substages`), set by the user;
* optional treatment of infested hosts (default off);
* new case definition: a new infestation of a host (W → I);
* ODE and chain-binomial solvers, equilibrium start, outputs, one scenario;
* the current model stays the default (`wound_model = "static"`), so existing
  results and the GP2014 regression are unchanged.

Later stages (not here): superinfestation resetting or extending the
infestation age; seasonal wounding (calving, shearing);
several host classes (livestock, wildlife, newborns); host mortality from
myiasis; detection and reporting of cases.

## 3. Model

### 3.1 Host states

Per node, with N_H = `host_density` × A hosts (A the node area):

| State | Meaning |
|---|---|
| S | hosts without an attractive wound |
| W | hosts with an attractive wound that is not infested |
| I = Σ_j I_j | hosts with an active infestation (eggs or larvae on the wound), in n_I sub-stages |

S + W + I = N_H at all times (no host births or deaths in stage 1).

```
            wounding ω                new infestation (case)
    S  ───────────────────►  W  ───────────────────────────►  I
       ◄───────────────────     ◄─ ─ ─ (not in stage 1) ─ ─ ─
            healing γ_W
    ▲                                                           │
    └──────── end of infestation: last sub-stage (q_I) or treatment (κ) ────────┘
          (larvae leave or the host is treated; the wound heals)
```

### 3.2 Rates

| Symbol | Meaning | From parameters |
|---|---|---|
| γ_W | healing of an uninfested wound (per day) | 1 / `wound_healing_days` |
| ω | wounding of a host without a wound (per day) | γ_W p₀ / (1 − p₀), p₀ = `wounded_fraction` |
| q_I | advance through one infestation sub-stage (per day) | n_I / `infestation_days`, n_I = `infestation_substages` |
| κ | treatment of an infested host, from any sub-stage (per day) | `treatment_rate` |
| α | attractiveness of an infested wound relative to an uninfested one | `infested_wound_attractiveness` |

Parameterising wounding by p₀, the fraction of hosts carrying an attractive
wound **when there are no flies**, makes the link to the current model
direct: at low fly density W = p₀ N_H, so the number of attractive wounds is
host_density × p₀ × A, the role H plays now. A calibration of
`wound_density` carries over as host_density × wounded_fraction.

### 3.3 Coupling to the flies

Attractive sites replace the fixed wound count:

    A_s = W + α I

Wound search (S8) keeps its form with A_s in place of H:

    Ψ = σ₀ [1 − exp(−a_s Δx A_s / N_s)]

Oviposition events and fertile batches are as now (S9, S17):

    E = Ψ f Σ_{g,m} W_{g,m}                                  (all batches)
    F = Ψ f Σ_{g,m} W_{g,m} (1 − e^{−C v_{g,m}})              (batches with viable eggs)

Batches are shared among wounds in proportion to attractiveness: a share
W / A_s goes to uninfested wounds, α I / A_s to infested ones. A fertile
batch on an uninfested wound starts an infestation; on an infested wound it
adds larvae only (superinfestation). Batches of females mated to sterile
males lay eggs that do not hatch and infest nobody (as decided for D-026).

    new infestations (cases)   C = F · W / A_s
    superinfestations              F · α I / A_s

### 3.4 Host equations (ODE)

    dS/dt   = −ω S + γ_W W + q_I I_n + κ I
    dW/dt   =  ω S − γ_W W − C
    dI_1/dt =  C − (q_I + κ) I_1
    dI_j/dt =  q_I I_{j−1} − (q_I + κ) I_j          (j = 2, …, n_I)

with I = Σ_j I_j. Without treatment an infestation lasts on average
`infestation_days`, with coefficient of variation 1/√n_I; every infestation
ends by the host returning to S (wound healed, decision 2).

with C = F W / A_s (C = 0 when A_s = 0). The fly equations are unchanged
except that Ψ uses A_s, and with treatment the larval hazard gains
`treatment_rate` (treating an infested host kills its larvae; all larvae
live in infested hosts, so in the mean field every larva is removed at the
treatment rate):

    μ_L = −ln(s_L) + treatment_rate [+ crowding, 3.5]

I ≤ N_H always: C is at most the flux of wounded hosts and every host is in
exactly one state.

### 3.5 Larval crowding

The current per-wound crowding divides larvae by all wounds, empty ones
included (the main weakness found in the review of S7). With hosts, the
natural denominator is the number of infested hosts:

    μ_L = −ln(s_L) + treatment_rate + γ_L L / (I + ε N_H)

`larval_crowding_mode = "per_infested_host"`; ε N_H (ε = 1e-6) keeps the
term finite when larvae outlive the infestations in the mean field and keeps
the model homogeneous in the area. Crowding stays off by default.

### 3.6 Chain-binomial step

Host counts are whole numbers. Within a step (start-of-step counts, as for
the other S15 aggregates):

1. Fertile batches B⁺ are drawn as now (S17), per female × mate pair.
2. Batches on uninfested wounds: n_W ~ Binomial(B⁺, W / A_s).
3. Newly infested hosts: n_W batches fall at random on the W wounded hosts;
   the number of distinct hosts hit, k, follows the occupancy distribution
   (drawn batch by batch: a batch hits a host not yet hit with probability
   (W − k) / W). So **k ≤ W** exactly. These are the cases of the step.
4. Healing among the other W − k: Binomial(W − k, 1 − e^{−γ_W dt}).
5. Infestation chain: moves as in S16 (exact transitions over dt for the
   Erlang chain), with treatment in the role of the hazard; hosts leaving
   the last sub-stage or treated return to S. New infestations enter I_1.
6. Wounding: Binomial(S, 1 − e^{−ω dt}).
7. All moves are applied at the end of the step; a host infested in this
   step cannot recover in the same step (bias of order dt, as for other
   arrivals).

Order: after eggs (S17) and before mating (S18), so the case count of a step
uses that step's batches. The R reference step and the C++ step follow the
same order and random numbers.

### 3.7 Equilibrium at constant temperature

Given the fertile-batch rate F, the host block has a unique steady state.
With b = C / W, the infestation hazard of a wounded host:

    I = b W D,   S = (γ_W + b) W / ω,   W = N_H / (1 + (γ_W + b)/ω + b D)

where D = I / C is the mean time a host stays infested: D = `infestation_days`
without treatment and D = [1 − (q_I / (q_I + κ))^{n_I}] / κ with treatment
(the sub-stages I_j follow the usual Erlang steady state). b solves
b W(b) (1 + α b D) = F, whose left side increases in b, so a
bracketed root exists. The fly part is solved as now by a one-dimensional
root search in the egg input; inside it, A_s (and so Ψ) comes from the host
steady state. The equilibrium therefore becomes a nested root search (as
for crowding and the Allee effect). R0 is unchanged: at low fly density
I → 0 and Ψ → σ₀.

**Area scaling** (D-028) stays exact: N_H, S, W, I, F are proportional to the
area and all host rates are per host, so the model is homogeneous of degree
one in the state and the area. `adult_females =` and `"threshold"` work in
host mode too.

**Feedback.** With α > 1 more infestations make the region more attractive
(A_s grows as hosts move from W to I), a positive feedback that may make the
host–fly system bistable; the equilibrium search must check for several
roots and the tests must include α > 1.

## 4. Parameters

New rows in `sg_param_table()`, group "hosts". `host_density` and
`wounded_fraction` have no default (decision 7): with `wound_model = "hosts"`
`sg_params()` stops until both are given, and its message suggests
host_density × wounded_fraction = 5 per km², the current `wound_density`, as
a starting point. In `"static"` mode they are not used. Values marked
*assumption* have no source yet.

| Name | Symbol | Default | Units | Basis |
|---|---|---|---|---|
| `wound_model` | – | `"static"` | – | `"static"` (current model) or `"hosts"` |
| `host_density` | N_H/A | none (required) | animals km⁻² | from livestock census for the region |
| `wounded_fraction` | p₀ | none (required) | – | fraction of hosts with an attractive wound without flies; product with host_density calibrated like `wound_density` |
| `wound_healing_days` | 1/γ_W | 7 | d | *assumption*; no source found for how long an uninfested wound (navel, husbandry wound, tick bite) stays attractive |
| `infestation_days` | n_I / q_I | 6 | d | mean duration without treatment; eggs hatch within about 24 h, larvae feed 5–7 days after hatching, then leave the wound [WOAH; CFSPH; FAO] |
| `infestation_substages` | n_I | 6 | – | Erlang sub-stages of an infestation (whole number 1–10, D-003); 6 as for larvae (CV 0.41) |
| `infested_wound_attractiveness` | α | 1 | – | infested wounds attract further gravid females and multiple infestations are common [WOAH; CFSPH]; up to 3000 larvae in one wound, about 15 batches of 200 [FAO]; magnitude unknown, so neutral default, sensitivity 1–10 |
| `treatment_rate` | κ | 0 | d⁻¹ | management; topical treatments are highly effective in calves [Brazil field study], so κ reflects how fast infestations are found and treated |
| `larval_crowding_mode` | – | `"none"` | – | adds `"per_infested_host"` (3.5) |

Related existing values from the same sources, for later decisions:

* egg batch size 10–400, average 200, about every 3 days, about 4 batches
  per female [FAO]; 343 eggs per mass [WOAH chapter]. With host-based cases
  the case count no longer depends on `clutch_size`, but `clutch_size` = 200
  would be a better default than one day's eggs (a separate decision);
* navel infestation of newborns can reach 100% where flies are abundant
  [APHIS]; 26% of untreated control calves in a field trial [fluralaner
  trial]; about 62–65% of beef calves in the rainy season [Brazil field
  study]: targets for stage 2 (newborn host class, seasonal wounding);
* untreated animals can die within 1–2 weeks [APHIS; CFSPH]: host mortality
  for a later stage.

## 5. Outputs

* New counters: `infestations` (cases, W → I) and `superinfestations`;
  `egg_batches` and `fertile_batches` are kept.
* `sg_hosts(sim)`: S, W, I per day (per node or km²), prevalence I / N_H.
* `sg_cases()` in host mode: `cases` = infestations; `cases_per_host` uses
  the model's host density (the `host_density` argument is then an error);
  `wound_fraction` becomes infestations per wounded host and day.
* `plot(sim, "hosts")`.
* `check_model()` reports the host rates and the largest host rate.

## 6. What the change gives

| Aspect | Now (`"static"`) | Stage 1 (`"hosts"`) |
|---|---|---|
| Case | fertile oviposition event | new infestation of a host |
| Animals infested at once | unbounded | at most N_H |
| Repeated infestation of an animal | – | yes: I → S → W → I |
| Effect of `clutch_size` on cases | strong | none (cases count animals) |
| Wounds | fixed | created, healed, infested |
| Superinfestation | counted as cases | adds larvae, no case |
| Larval crowding | per wound, empty ones included | per infested host |
| Treatment | – | ends infestations, kills larvae |

Example to demonstrate (scenario 04): a farm with 2 animals. The number
infested at any time is 0, 1 or 2; over a year the same animals are infested
many times; with `"static"` the same flies would report many simultaneous
cases.

## 7. Tests

* `"static"` reproduces all current results (existing tests, GP2014
  regression, snapshots) bit for bit.
* S + W + I = N_H in every ODE step and every stochastic replicate.
* I ≤ N_H and new infestations ≤ W in every stochastic step (also N_H = 2).
* Without flies the host block relaxes to W = p₀ N_H, I = 0.
* At low fly density the host model gives the static model's Ψ with
  H = host_density × wounded_fraction.
* The equilibrium is a fixed point of the ODE; the deterministic map of the
  stochastic step converges to the ODE as dt → 0.
* Area scaling: `adult_females =` exact in host mode.
* α > 1: equilibrium search finds the stable root; bistability checked.
* Treatment: higher `treatment_rate` lowers prevalence and larvae; the
  steady-state infested time matches D = [1 − (q_I/(q_I + κ))^{n_I}] / κ.
* Infestation chain: without treatment the mean duration is
  `infestation_days` and its variance `infestation_days`² / n_I.
* `wound_model = "hosts"` without `host_density` or `wounded_fraction`
  stops with a message naming both.
* Flow balance (sg_flow_balance) still closes.

## 8. Implementation plan (work packages)

One branch (`host-wound`), one commit per work package, one PR at the end,
as for the milestones. Each WP ends with tests, `dev/style.R`, lint and the
full `dev/check.R`. The R reference code (the oracle) is written before the
C++ code that must reproduce it (as in M2–M4).

**Completion criteria for stage 1**

1. `wound_model = "static"` gives numerically identical results to the
   current package (all existing tests and the GP2014 regression pass;
   only parameter-print snapshots and model hashes change, because the
   parameter list grows).
2. Hosts are conserved (S + W + Σ I_j constant) in the ODE and in every
   stochastic replicate; I ≤ N_H and new infestations ≤ W in every
   stochastic step.
3. C++ = R reference: ODE to 1e-8; chain-binomial with the same draws.
4. The host-mode equilibrium is a fixed point of the ODE; `adult_females =`
   is exact in host mode; R0 is unchanged.
5. Scenario 04 shows a 2-animal farm with at most 2 animals infested at once
   and repeated infestations over a year.

### WP1. Parameters

* New rows in `R/param-table.R` (group "hosts"): `wound_model`
  (`"static"`/`"hosts"`), `host_density`, `wounded_fraction`,
  `wound_healing_days`, `infestation_days`, `infestation_substages`,
  `infested_wound_attractiveness`, `treatment_rate`; new choice
  `"per_infested_host"` for `larval_crowding_mode`.
* Required values (decision 7): `host_density` and `wounded_fraction` default
  to NA; the parameter system learns "NA allowed, required when
  `wound_model = "hosts"`". `sg_params(wound_model = "hosts")` without them
  stops with a message naming both and suggesting
  host_density × wounded_fraction = 5 per km².
* Cross-checks: `wounded_fraction` in (0, 1); `infestation_substages` a whole
  number 1–10 (D-003); `"per_infested_host"` only with `"hosts"`; host
  parameters given in `"static"` mode give a warning that they are unused.
* `gp2014` preset stays `"static"`.
* Tests: defaults, validation messages, required values, param table and
  print snapshots (updated).

### WP2. State layout and model assembly

* Host block appended after F in `.sg_layout()`: H_S, H_W, H_I1 … H_In per
  node and host class (`n_host` = 1); length 0 in `"static"` mode, so static
  state vectors are unchanged. State names e.g. `H_W[1]`, `H_I[1,3]`;
  `sg_state_index()` gains the host blocks.
* `sg_model()`: host counts per node (N_H = host_density × area), derived
  rates (ω, γ_W, q_I, κ) in `model$counts` / the setup; `print()` shows the
  host model; hash covers the new inputs.
* `sg_init()`: `"constant"` places hosts at the fly-free steady state
  (W = p₀ N_H, S = (1 − p₀) N_H, I = 0); `"state"` accepts the longer
  vector; state rounding in stochastic runs keeps each host compartment a
  whole number (the host total is then fixed per replicate).
* Model view (`R/model-view.R`, `src/model_view.*`): flag, host rates, n_I,
  block offsets.
* Tests: layout lengths and names in both modes, static layout unchanged,
  initial host state, view contents.

### WP3. R reference ODE (oracle)

* `R/reference-rates.R`: attractive sites A_s = W + α I replace H in Ψ (S8);
  new infestations C = F W / A_s; host derivatives (section 3.4, Erlang
  chain); larval hazard + `treatment_rate` and the `"per_infested_host"`
  crowding (3.5).
* New counters: `infestations`, `superinfestations`, `infestations_ended`,
  `treatments` (host flows), alongside the existing ones.
* `simulate_ode(backend = "R")` runs host mode.
* Tests: conservation; relaxation to W = p₀ N_H without flies; at low fly
  density Ψ equals the static Ψ with H = host_density × wounded_fraction;
  infestation chain mean and variance; treatment lowers prevalence and
  larvae; `"static"` unchanged.

### WP4. C++ ODE kernel

* `src/rates.cpp` / `rates.h`: the same terms in the aggregates and
  derivatives; host rates in `max_rate` and the stiffness diagnostics.
* Backend equivalence to 1e-8 (R vs C++) in both modes and with treatment,
  crowding and α > 1; speed check on the McAllen year (host mode should
  cost little extra).

### WP5. Chain-binomial step

* R reference step (`R/cb-reference.R`), then C++ (`src/cb_step.*`), in the
  order of section 3.6: fertile batches → batches on uninfested wounds
  (binomial) → distinct hosts hit (occupancy draw, k ≤ W) → healing →
  infestation chain moves (S16 with treatment as the hazard) → wounding;
  moves applied at the end of the step.
* Tests: R and C++ steps identical with the same seed; I ≤ N_H and k ≤ W in
  every step (including N_H = 2); deterministic map → ODE as dt → 0;
  statistical test of the occupancy draw; replicate means vs ODE.

### WP6. Initial states, equilibrium, scaling, calibration

* `"burn_in"` runs the host block with the flies.
* `"equilibrium"` and `"threshold"`: nested root search with the host steady
  state (section 3.7, mean infested time D); with α > 1, scan for several
  roots and return the stable upper one, with a message if there are more.
* `adult_females =` in host mode (exact by homogeneity).
* `calibrate_wound_density()` in host mode: scales `host_density` (exact
  without an Allee effect, as now for H) and says so; to be confirmed.
* Tests: equilibrium is an ODE fixed point; area scaling; R0 unchanged;
  α > 1 cases.

### WP7. Outputs and checks

* `sg_hosts(sim)`: S, W, I per day, per node or km², prevalence I / N_H
  (replicates for stochastic runs).
* `sg_cases()` in host mode: cases = infestations, `superinfestations`
  column, `cases_per_host` from the model's host density (`host_density`
  argument an error); `wound_fraction` = infestations per wounded host and
  day. Static mode unchanged.
* `plot(sim, "hosts")`; `check_model()` reports host rates and the largest
  host rate × dt; `sg_flow_balance()` adds a host conservation check.
* Tests and plot snapshots.

### WP8. Scenarios, documentation, PR

* `inst/scenarios/04_host_limited_cases.R`: a farm with 2 animals (fixed
  area and host density so N_H = 2), stochastic runs showing 0–2 animals
  infested at once and repeated infestations over a year, compared with the
  static model's case counts for the same flies.
* MODEL_SPEC: host section (equations of sections 3.3–3.7), S7/S8/S17
  amendments, parameter table; DECISIONS: D-030 (this note's decisions);
  NEWS; scenario README; pkgdown reference; roxygen for new functions.
* Full check, then the PR.

### Points to settle during implementation

* Names of the host blocks and counters (proposed above).
* Whether `calibrate_wound_density()` scales `host_density` in host mode or
  a separate `calibrate_host_density()` is clearer.
* How many roots to search for with α > 1 (bistability), and what to report.

## 9. Decisions (made by the user, 2026-10-16)

1. **Wounding parameterisation:** `wounded_fraction` + `wound_healing_days`;
   host_density × wounded_fraction takes the role of `wound_density`.
2. **After an infestation:** the host returns to S (wound healed).
3. **α default:** 1, neutral, until a value is found.
4. **Infestation duration:** Erlang chain with its own mean
   (`infestation_days`, default 6) and number of sub-stages
   (`infestation_substages`, default 6), both set by the user, not tied to
   the egg and larval chains.
5. **Treatment:** `treatment_rate` included in stage 1, default 0.
6. **Default model:** `"static"` stays the default.
7. **Host numbers:** `host_density` and `wounded_fraction` have no defaults;
   host mode requires the user to set them.

## References

* WOAH. New World screwworm (*Cochliomyia hominivorax*) and Old World
  screwworm. Terrestrial Manual chapter 3.1.14:
  <https://www.woah.org/fileadmin/Home/fr/Health_standards/tahm/3.01.14_SCREWW.pdf>;
  disease card:
  <https://www.woah.org/fileadmin/Home/eng/Animal_Health_in_the_World/docs/pdf/Disease_cards/SCREWWORM.pdf>
* WOAH Terrestrial Manual chapter 3.1.13 (egg masses):
  <https://www.woah.org/fileadmin/Home/eng/Health_standards/tahm/3.01.13_SCREWW.pdf>
* CFSPH (Iowa State University). Screwworm myiasis factsheet:
  <https://www.cfsph.iastate.edu/Factsheets/pdfs/screwworm_myiasis.pdf>
* FAO. Screwworm flies as agents of wound myiasis:
  <https://www.fao.org/4/u4220t/u4220t07.htm>
* USDA APHIS. New World screwworm:
  <https://www.aphis.usda.gov/livestock-poultry-disease/cattle/ticks/screwworm>;
  disease response strategy (2025):
  <https://www.aphis.usda.gov/sites/default/files/nws-myiasis-disease-strategy.pdf>
* Brazil field study: number of rainy days and navel myiasis in beef calves;
  efficacy of antiparasitics: <https://pubmed.ncbi.nlm.nih.gov/36252323/>
* Fluralaner field trial, umbilical myiasis in newborn calves:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC13520140/>

The values above were taken from these sources' summaries; the primary
literature behind them (e.g. on wound attractiveness and healing) still
needs to be read before the assumptions are replaced.
