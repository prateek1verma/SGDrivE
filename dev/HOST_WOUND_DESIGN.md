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

## 8. Implementation sketch (work packages)

1. Parameters, validation, presets (`"static"` default), param table.
2. State layout: a host block with S, W and I_1 … I_n per node and host
   class (`n_host` = 1), state names, `sg_init()` methods place hosts at
   W = p₀ N_H (constant start) or at the host steady state (burn-in,
   equilibrium).
3. R reference rates and ODE (S7, S8, host equations, counters).
4. C++ rates and ODE; backend equivalence tests.
5. Chain-binomial step (R reference and C++): host moves, occupancy draw.
6. Equilibrium, threshold and `adult_females` with the host block.
7. Outputs (`sg_hosts()`, `sg_cases()`, plots, `check_model()`).
8. Scenario 04 (small farm), MODEL_SPEC section, D-030, NEWS.

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
