# Multiple nodes and fly movement, v0.2 (design note)

Status: **accepted**, 2026-11-04 (phase 0 done); decisions of section 10
made by the user on 2026-11-02 and 2026-11-04; not implemented (phases
1–4 open). Equations in MODEL_SPEC section 8a (N1 movement terms, N2
losses, N3 clock, N4 chain-binomial step, N5 node extinction rule);
decision D-032 (DECISIONS.md), which also records the choices made in
writing the specification. Biology from the sources in section 2,
including Medical and Veterinary Entomology 23, Suppl. 1 (2009).

## 1. Why

ScrewDrivE v0.1 has one node (MODEL_SPEC section 1). The questions that matter
for New World screwworm (NWS) control are spatial: how far and how fast an
infestation spreads north each season, how wide a sterile-fly barrier must
be, whether a node far from the release area is reinvaded, and how county
case series (Texas 1962) are reproduced. ROADMAP lists v0.2 as "multiple
nodes and fly movement" and v0.4 as host movement and larval transport.

## 2. What the biology says

| Topic | Evidence | Consequence for the model |
|---|---|---|
| Local adult dispersal | Mayer & Atzeni (1993) pooled mark–release–recapture trials: displacement fits a **Cauchy** (fat-tailed) distribution; **females disperse farther than males**; typical female range about **3 km** in favourable habitat, up to **about 8×** farther where hosts are scarce or habitat poor. | Rates by class and by node; users can build fat-tailed movement with their own kernel. |
| Maximum range | Marked irradiated flies recovered up to **180 miles (≈ 290 km)** from release (Hightower, Adams & Alley 1965; Spradbery 1994), studied after the rapid Texas spread of spring 1963. | Small long-distance rates. |
| Emigration in the field | Matlock & Skoda (2009), eight mark–release–recapture studies in Costa Rica (19,573 marked flies, 4–7 bait stations per study): apparent daily survival φ = 0.798 (both sexes). A diffusion model fitted to the same data separates death from emigration: true daily survival **0.947 (♀)** and **0.886 (♂)**. So about **16% of females and 10% of males leave a station array per day**. | A leave rate at node scale, usable directly as δ (section 5.2); more female than male movement. |
| Growth rate | Same study: field λ = **1.005 d⁻¹** (doubling 139 days; range 0.942–1.011); λ from rearing-facility life history **1.25 d⁻¹** (doubling 3 days). | Bounds r for the spread-rate check below. |
| Population structure | Robinson et al. (2009), Lyra et al. (2009), Torres et al. (2007): low genetic differentiation among mainland South American populations (high dispersal) but **no isolation by distance** in Uruguay, explained by **local extinction and recolonisation**; Caribbean islands genetically isolated (water barriers); mtDNA vs nuclear differences suggest **female-biased dispersal**. | Stochastic metapopulation runs with local extinction; barriers as missing or weak edges; higher female rates. |
| Introductions | Virgin Islands reinfested in 1977, probably by **wind** (Rawlins 1985); Libya 1988 by **moved infested animals** (Lindquist et al. 1992), "the most common way" of introduction (Robinson et al. 2009). | Long-distance layers (phase 4) and host transport (v0.4). |
| Seasonal invasion | Overwintering in south Texas, Florida, Arizona, California and Mexico; spread north each spring at about **80–160 km per generation** (FAO); severe winters (1961–62) reduce the source. | Node-specific weather; source nodes outside the study area; northern nodes die out each winter and are recolonised. |
| Wind-borne migration | Gutierrez, Ponti & Arias (2019): southwestern US outbreaks in years when **North American monsoon winds** carried adults north from northern Mexico and south Texas. | A seasonal, directional layer on the calendar clock (phase 4). |
| Male behaviour | Krafsur (1978): males form mating aggregations on sunlit mesquite and thorny brush, territorial, independent of cattle. Males and virgin females stay in vegetation; **mated females orient to hosts**, about 3–4 days after mating. | Males move least; a hook for rates by adult segment (default 1). |
| Human-mediated spread | Nuevo León case (September 2025): a calf moved from infested southern Mexico. Livestock movement is the main explanation for the jumps through Central America. A 2026 PLOS ONE model (local diffusion ≈ 7.2 km per week) also needed host transport for long jumps. | Transport of infested hosts with their larvae (v0.4). |
| Barrier operations | Aerial release lanes about 1.6 km apart; 20 km infested zone in the 2026 Texas response; 40–50 million sterile flies per week in the Darién barrier (Robinson et al. 2009; McDonagh et al. 2009). | Releases per node; a barrier helper (phase 4). |

**Spread-rate check.** A Fisher–KPP front moves at c = 2√(rD). With a
lifetime displacement of about 3 km, D is of order 0.1–0.5 km² d⁻¹ and the
front speed well below 1 km d⁻¹: even with the optimal growth rate
r = ln 1.25 ≈ 0.22 d⁻¹ (Matlock & Skoda 2009) and D = 0.5 km² d⁻¹,
c ≈ 0.66 km d⁻¹, and with the field r ≈ 0.005 d⁻¹ it is about 0.1 km d⁻¹.
The observed 80–160 km per generation (3–7 km d⁻¹ at summer generation
times) needs long-distance movement (wind, long flights of gravid females,
moved livestock). The model therefore needs two scales of movement: local,
diffusion-like movement of all adults, and rare long-distance moves.

**Consistency check.** Matlock & Skoda's female leave rate,
−ln(0.84) ≈ 0.17 d⁻¹, on a 2-D lattice (δ = 4D/h²) with an array span of
about 3 km gives D ≈ 0.4 km² d⁻¹, inside the range above (rough: the
array geometry is approximate).

Unverified: only abstracts or summaries of Mayer & Atzeni (1993) and the
PLOS ONE (2026) paper were read. The Matlock & Skoda diffusion analysis was
"in preparation" in 2009; look for the published version, which would give
D by sex directly.

## 3. Design principle: raw rates, not a built-in kernel

The model's only movement input is **rates and destination probabilities**
(the generator of a continuous-time Markov chain over nodes), as in
MGDrivE2. The package does not choose a dispersal kernel.

Why (decision 3, section 10):

* What the user gives is exactly what runs; no hidden approximation.
* The processes that matter most for NWS are not distance kernels: water
  barriers, directed monsoon wind and livestock trade networks are networks
  or directed flows, which a matrix holds directly.
* A "lifetime displacement in km" is not one well-defined number: the
  growth of displacement with time (diffusive √t for Gaussian steps, linear
  t for Cauchy steps), the age mixture of recaptured flies and trap layout
  all enter, so a built-in kernel would look more rigorous than the
  evidence is.
* The best NWS field number (Matlock & Skoda 2009) is already a leave rate
  at node scale.
* Direct comparison with MGDrivE2 (two-node cross-check).

The weakness of raw rates is **scale dependence**: the same biology needs
different rates on different grids (D = 0.5 km² d⁻¹ gives δ = 0.01 d⁻¹ at
10 km spacing but 0.0004 d⁻¹ at 50 km), and a rate copied from one grid to
another is wrong without warning. Two additions remove this weakness
without committing the package to a kernel:

1. **Constructors with a user-supplied kernel** (section 6.1): the user
   passes any function of distance; the package does the bookkeeping. Exact
   helpers convert D, emigration fractions and lifetime probabilities to
   rates. All return ordinary objects that users can inspect and change, and
   a scalar such as D can drive sensitivity analyses.
2. **A diagnostic that translates rates back into biology**
   (`check_movement()`, section 6.2): implied lifetime leave probability,
   path length and D per node and class, so a 1 km-scale rate on a county
   grid shows up as an implausible lifetime displacement.

## 4. Scope of v0.2

In v0.2:

* N nodes, each with its own area, coordinates, weather (forcing), wound or
  host density and releases, and a role (study or source node); one node
  reproduces v0.1 exactly;
* adult movement (blocks M, U, F) between nodes, linear in the counts, in
  both solvers;
* movement given as a **list of layers**; phases 1–3 build and test one
  layer, phase 4 several (long distance, seasonal schedules);
* open or closed landscape edges through the destination probabilities;
* constructors, conversion helpers and `check_movement()`;
* outputs per node and spread metrics;
* first landscape: a transect.

Not in v0.2: movement of eggs, larvae or pupae; transport of infested hosts
(v0.4); host-dependent dispersal (later option); batch events of many flies
moved together (deferred until a scenario needs them).

## 5. Model

### 5.1 State layout

The index order of CONVENTIONS is unchanged:
`[node, host_class, genotype, (mate_genotype), sub_stage]`. The node is the
slowest index, so node n's blocks start at n × (states per node), and every
v0.1 offset is an offset within a node. Host blocks (D-030) are per node.

### 5.2 What moves and at what rate

* Adults only: M_{g,a}, U_{g,a}, F_{g,m,a}, classes s ∈ {M, U, F}. Unmated
  females move (decision 1). Eggs, larvae and pupae are in wounds or soil
  and move only with hosts (v0.4).
* A mover keeps its genotype g, its mate genotype m (F) and its
  physiological-age sub-stage a: X[i, g, (m), a] → X[j, g, (m), a]. The
  maturity thresholds of the segmented adult chains (S6, D-018) are
  unaffected, and a mated female carries her mate's sperm, so a single
  mated female can found a population in an empty node.

A **movement layer** ℓ has:

* leave rates δ^ℓ_{s,i} [d⁻¹] for each class s and node i, given **per day
  at the reference temperature T_ref** (one model parameter,
  `movement_t_ref`, default 27 °C, near T_opt; the same for all layers);
* destination probabilities π^ℓ_{ij}, sparse, π_ii = 0, **Σ_j π_ij ≤ 1**;
  one matrix for all classes. Class differences ("females disperse
  farther") come from the rates: a higher female rate in a layer, and with
  a local and a long-distance layer, class-specific rates in each give
  class-specific effective kernels. With the diffusion constructor D
  cancels from π anyway, so per-class matrices would add an option and
  nothing else;
* a clock (5.3);
* a schedule (phase 4): one daily multiplier m^ℓ(t) ≥ 0 for the whole
  layer, the same in every node (e.g. 1 in July–September and 0
  otherwise for a monsoon layer); default 1 every day.

**Time variation is kept minimal** (decision 8). In phases 1–3 rates vary
in time only through the clock (5.3), computed in C++ from the forcing it
already holds (Δx per node and day); no extra table is stored. Phase 4
adds only the daily multiplier above, one numeric vector per layer.
Deferred until a scenario needs them: multipliers that differ by node and
day, multipliers given as R functions, destination matrices that change
over time (e.g. monthly wind directions) and time series of transition
probabilities. Each can be added later on top of this without changing
the interface.

A hook for rates by adult segment (e.g. laying females, which search for
hosts), c_{s,a}, is in the specification with default 1 and is not exposed
in the interface until there is evidence for it.

### 5.3 Movement clock (decision 2)

    κ^ℓ(T) = min(κ_max, Δx(T) / Δx(T_ref))   clock = "physiological" (default)
    κ^ℓ(T) = 1                                clock = "calendar"

with T_ref = `movement_t_ref` (default 27 °C) and κ_max =
`movement_kappa_max`, one value for the model: **1 by default (capped)**,
`Inf` for uncapped degree-day scaling, values in between allowed. The cap
is a multiplier, not a temperature, because Δx rises to a maximum near
36 °C and falls above it, so a cap temperature would be ambiguous. κ_max
is the same for all layers, so it does not add groups when layers are
merged (5.4).

* **Capped** (κ_max = 1, default): at and above T_ref flies move at the
  given daily rate (calendar time); below T_ref movement declines with Δx
  and stops at the development threshold T_L = 14.5 °C (with
  `dev_rate_floor` > 0 it keeps a small positive value). κ falls below 1
  again only above about 41 °C, where Δx drops under Δx(T_ref).
* **Uncapped** (κ_max = `Inf`): movement per day follows Δx at all
  temperatures, up to about 1.57× the reference rate near 35–36 °C.
* **Calendar**: for wind-borne layers (wind transport does not scale
  with fly development) and for the MGDrivE2 cross-check.

Rates are given at T_ref, so field numbers (e.g. a daily leave fraction,
Matlock & Skoda 2009, measured at about 25–28 °C) go in directly. 27 °C is
a reference point near T_opt (27.2 °C), not a biological limit.

**What the clock controls.** Within a node, leaving competes with death
and ageing, so the clock mainly sets how **lifetime** dispersal changes
with temperature. Female leave rate 0.17 d⁻¹ at T_ref, default adult chain
(n_A = 10, X_A = 202.8 dd), temperature mortality μ_T of S2, constant
temperature (computed 2026-11-04):

| T (°C) | Adult life (d) | Δx / Δx(T_ref) | Moves per lifetime: calendar | capped | uncapped |
|---|---|---|---|---|---|
| 18 | 24.6 | 0.28 | 4.19 | 1.17 | 1.17 |
| 22 | 23.3 | 0.60 | 3.96 | 2.38 | 2.38 |
| 27 | 16.5 | 1.00 | 2.81 | 2.81 | 2.81 |
| 30 | 13.1 | 1.24 | 2.24 | 2.24 | 2.77 |
| 35 | 9.4 | 1.57 | 1.59 | 1.59 | 2.50 |

* Calendar: most lifetime movement in cool weather; flies at 18–22 °C
  leave at the full summer rate, also on days too cold to fly.
* Capped: lifetime movement falls in heat (life shortens while the daily
  rate is fixed).
* Uncapped: lifetime movement about constant above about 25 °C (the daily
  rate rises as calendar life shortens).
* Below T_ref capped and uncapped are identical. Both give less lifetime
  movement in the cold, because the calendar-time mortality μ_T ends life
  before many degree-days accumulate; scaling with Δx keeps dispersal only
  roughly consistent with ageing in degree-days (D-002).

**Evidence and assumptions.** No NWS study relates movement rate to
temperature, and no study supports 27 °C or the cap.

* In the cold, less movement has qualitative support (activity falls
  below about 21 °C in autumn and flies do poorly below about 7 °C;
  Oklahoma State and Drovers extension notes).
* In heat the two physiological variants differ. Capped claims nothing
  beyond the measured daily rate; at daily means above about 30 °C midday
  temperatures are high enough that flies may shelter rather than fly
  more. Uncapped keeps lifetime dispersal nearly independent of
  temperature, a reasonable null. The data do not separate them.
* The shape of the decline below T_ref is borrowed from the development
  curve as a convenience; a flight-activity curve has its own threshold.

**Sensitivity test (phase 3).** The first transect runs with all three
clocks (calendar, capped, uncapped). Capped vs uncapped differ only on days
above T_ref (south Texas summers); calendar vs physiological differ mainly
on cool days. If winter die-back and spring spread barely change, the
choice does not matter in practice.

### 5.4 Equations (ODE)

**(N1) Movement terms**, added to the v0.1 derivatives of every adult
compartment X of class s, summed over layers:

    dX_j/dt += Σ_ℓ [ Σ_i κ^ℓ_i(t) δ^ℓ_{s,i} π^ℓ_{ij} X_i − κ^ℓ_j(t) δ^ℓ_{s,j} X_j ]

where κ^ℓ_i(t) uses node i's temperature. Linear, so R₀ and the ODE keep
their form.

**(N2) Losses from the landscape.** The part 1 − Σ_j π_{ij} of the flies
leaving node i leaves the landscape:

    lost_i rate = Σ_ℓ κ^ℓ_i δ^ℓ_{s,i} (1 − Σ_j π^ℓ_{ij}) X_i

Rows summing to 1 are walls (closed edges); rows summing below 1 are open
edges. The network total changes by movement only through N2; with closed
edges it is unchanged.

**Layers are merged inside the model.** Rates of layers that share a clock
and a schedule add, so R merges them before the run into one group:
δ^G_{s,i} = Σ_{ℓ∈G} δ^ℓ_{s,i} and π^G_{ij} = Σ_{ℓ∈G} δ^ℓ_{s,i} π^ℓ_{ij} /
δ^G_{s,i} (per class, as rates differ by class). Layers are a concept of
the interface; the kernel sees only groups, at most one per (clock,
schedule) pair: in practice one in phases 1–3 and two or three in phase 4
(e.g. local on the physiological clock, long distance on the calendar
clock with a monsoon schedule), however many layers the user gives. The
ODE then needs one movement pass per group, and the chain-binomial step
draws once per compartment (5.6).

### 5.5 Source nodes (inflow from outside)

Inflow from outside the study area comes from **source nodes**: ordinary
nodes with their own weather, wound or host density and initial state,
marked `role = "source"`. They are simulated like any other node and feed
the study nodes through the layers. Outputs and spread metrics report study
nodes by default and source nodes on request. Example: for a south-Texas →
Kansas transect, a source node with northern-Mexico or Rio Grande valley
weather at the southern end.

### 5.6 Chain-binomial step

After the chain moves (S16) and before arrivals, for each node i in index
order and each adult compartment X of class s, with the total leave rate
over the groups of 5.4, λ_{s,i} = Σ_G κ^G_i m^G δ^G_{s,i}, and the
weights w_G = κ^G_i m^G δ^G_{s,i} / λ_{s,i}:

    E ~ Binomial(X, 1 − exp(−λ_{s,i} dt))                           emigrants
    p_ij = Σ_G w_G π^G_{ij},  p_lost = Σ_G w_G (1 − Σ_j π^G_{ij})
    (E_j)_j, E_lost ~ Multinomial(E; p_{i·}, p_lost)                 destinations, lost

drawn as sequential conditional binomials over destinations j in index
order (the union of the groups' edges, fixed at setup), then lost. One
emigration draw per compartment, whatever the number of layers, and the
result does not depend on the order of the layers; p is recomputed only
when κ or m changes (once per day). This is first-order operator
splitting (movement does not change age or genotype, as with mating
before ageing in D-022); emigrants face the destination's hazards from the
next step. Emigration, immigration
and losses are counted per node; with closed edges immigration and
emigration sum to the same network total exactly. The deterministic map
converges to the ODE as dt → 0.

### 5.7 Small populations in the ODE ("atto-flies")

In the ODE, movement puts tiny amounts (e.g. 10⁻¹² of a fly) into every
connected node almost at once, and a node dying out in winter decays
towards zero without reaching it. Without a mating Allee effect (the
default `male_density_half_mating = 0`) these fractions grow again in
spring. v0.1 has no population floors and uses the extinction threshold
only in summaries (MODEL_SPEC section 10, "Extinction"), so an ODE
landscape would show winter die-back and spring recolonisation made of
fractions of a fly.

**Node extinction rule** (ODE, optional): at the end of each day, a node
whose fly population (E, L, P, M, U, F, all genotypes) totals less than
`extinction_threshold` (default 1 individual) is set to zero. In the
host-wound model, its infested hosts end their infestation as when the
larvae leave and are counted in infestations_ended (MODEL_SPEC N5). The amount
removed is counted per node (a counter), so the conservation checks still
balance. The ODE is integrated one day at a time, so the rule costs
nothing and does not disturb the step-size control.

* Off by default, so one-node runs reproduce v0.1 exactly
  (`simulate_ode(..., node_extinction = TRUE)` to switch it on).
  `check_model()` notes multi-node ODE runs without it.
* Spread, die-back and recolonisation questions use either the ODE with
  the rule on or the chain-binomial solver (integer counts, no
  atto-flies). The phase 3 transect runs both.
* Off by default also for more than one node (decision 12, user,
  2026-11-04): the rule stays an optional feature, not a default that
  depends on the run.

### 5.8 Later options (not v0.2)

* Female leave rate rising as wounds get scarce (Mayer & Atzeni's
  favourability effect), e.g. a factor of Ψ (S8); nonlinear, so it would
  break the linear R₀ calculation. Fixed habitat differences (e.g.
  rangeland vs cropland, the "up to 8× farther in poor habitat" of
  section 2) are already covered by leave rates per node; this option is
  needed only when host numbers change within a run.
* Batch events of many flies carried together (MGDrivE2
  `batch_migration()`), if correlated arrival turns out to matter.
* Richer time variation (5.2): node × day multipliers, R-function inputs,
  destination matrices by time block, probability time series.
* Transport of infested hosts with their larvae, quarantine (v0.4).

## 6. Interface

### 6.1 Objects and constructors

* `sg_landscape(nodes, forcing, ...)`: per node id, area, coordinates
  (optional), forcing (one `sg_forcing` per node or one shared), wound or
  host density, role ("study" or "source"). Default: one study node.
* `sg_movement_layer(rate, dest, clock = "physiological")`: `rate` a named
  list or matrix of leave rates per day at T_ref (`movement_t_ref` in
  `sg_params()`; the cap of the physiological clock is
  `movement_kappa_max`, 5.3) for `male`, `virgin`, `mated` (scalar or one per node);
  `dest` a sparse destination matrix with rows ≤ 1, shared by the three
  classes. Validates dimensions, π_ii = 0, rows ≤ 1, rates ≥ 0.
* `sg_model(..., landscape = , movement = )`: `movement` is one layer or a
  list of layers; phases 1–3 accept one layer (a clear error for more,
  until phase 4). Single-node calls without these arguments are unchanged.
* `sg_releases(..., node = )`: releases per node.

Any network is given through `sg_movement_layer(rate, dest)`: a chain, a
grid, county neighbours, a trade network or long-distance links. The
constructors below only compute `rate` and `dest` for common cases.

Constructors and conversions (all return ordinary layers or numbers):

* `sg_movement_diffusion(area, border, dist, D, outer_border = NULL,
  outer_dist = NULL)`: local diffusion on **any** network of nodes by the
  finite-volume formula

      rate i → j = D · L_ij / (d_ij · A_i)        [d⁻¹ at T_ref]

  with `area` A_i (km²), `border` the matrix of shared border lengths L_ij
  (km, zero for nodes that do not touch), `dist` the centroid distances
  d_ij (km) and `D` (km² d⁻¹, one value per class: male, virgin, mated).
  Then δ_i = Σ_j rate_ij and π_ij = rate_ij / δ_i. Closed edges by default.
  Open edges: `outer_border` gives each node's border with the outside
  (km) and `outer_dist` the distance to a mirrored node outside (default
  twice the centroid-to-border distance, from `area` and `outer_border` if
  not given); the outward rate D · L_out / (d_out · A_i) is added to δ_i
  and becomes the lost share of π. Inputs are plain matrices; users compute
  them from shapes with their own tools (e.g. `sf::st_touches()`, shared
  border lengths, `sf::st_distance()` between centroids), so the package
  does not depend on `sf`.
  Special case, a transect of cells with spacing h: L = cell width,
  A = h × width, d = h, so each neighbour gets D/h² and an interior node
  leaves at δ = 2D/h²; a closed end leaves at D/h² to its one neighbour;
  an open end at 2D/h², half of it lost.
* `sg_movement_from_distance(dist, kernel, rate, cutoff = 1e-6, edges)`:
  `kernel` is any user function of distance giving relative weights;
  weights below `cutoff` dropped (sparsity); rows normalised to 1 (closed)
  or to the share of the kernel mass that falls on nodes (open).
* `sg_movement_from_prob(P)`: a daily transition-probability matrix to a
  layer, δ_i = −ln P_ii, π_ij = P_ij / (1 − P_ii) (rows of P below 1
  become losses), as MGDrivE2 `movement_prob2rate()`.
* Conversions: daily emigration fraction f → δ = −ln(1 − f) (e.g. Matlock
  & Skoda); lifetime leave probability → δ, computed exactly for the adult
  chain at T_ref (inverse of 6.2).

### 6.2 `check_movement()`

For each layer, node and class, at T_ref:

* daily leave probability 1 − e^{−δ};
* lifetime leave probability 1 − E[e^{−δτ}], with τ the adult lifespan of
  the chain at T_ref (exact from the Erlang chain with its hazard);
* expected number of moves per lifetime δ E[τ] and, with coordinates,
  expected path length per lifetime (δ E[τ] Σ_j π_ij d_ij) and implied
  D = δ Σ_j π_ij d_ij² / (2 · dimension);
* share of movers lost from the landscape;
* optionally, the lifetime numbers above over a range of constant
  temperatures, with the clock and κ_max in use (the table of 5.3), so the
  effect of the clock on lifetime dispersal is visible.

Also included in `check_model()`. Warnings for a lifetime leave probability
above 0.99 and for κδ dt above 0.5 in the chain-binomial step; the
largest-rate check includes κδ.

## 7. Outputs

### 7.1 Recorder: only what the user asks for

A run stores only the quantities requested, computed in C++ during the
run, not the full state (decision 9).

* **Linear quantities.** Almost every output is a weighted sum of states:
  adult females in a node, wild-type males per county, larvae over all
  nodes, copies of an allele (count × copies per genotype), mated females
  carrying sterile sperm, infested hosts. R turns each request, before the
  run, into a row of a sparse weight matrix W_x [quantities × states]; at
  each output time C++ computes y = W_x · x (one sparse matrix–vector
  product, about the cost of one rate evaluation). A new quantity is a new
  row built in R, not new C++ code.
* **Counters.** The daily counters (eggs laid, fertile batches or cases,
  matings, deaths by cause, emigrants, immigrants, losses), per node, are
  combined the same way with a second matrix W_c [quantities × counters].
  Counters can be summed over periods (e.g. weeks), so weekly cases per
  county compare directly with the Texas 1962 data.
* **Nonlinear quantities** (allele frequency, prevalence, wild-type
  share) are computed in R afterwards from recorded linear parts (e.g. a
  numerator and a denominator).
* **Final state** is always kept, whatever is recorded, so a run can be
  continued (e.g. after a burn-in).
* **Replicates.** In the chain-binomial solver the recorder applies to
  every replicate; memory is estimated from the recorder's size and
  checked against `max_memory` before the run (both solvers).

Interface (sketch):

```r
rec <- sg_record(
  sg_quantity("adults", stage = c("M", "U", "F"), by = "node"),
  sg_quantity("wild_F", stage = "F", genotype = "WW", node = study_nodes),
  sg_counter("cases", by = "node", every = "week"),
  sg_alleles("A", by = "node")
)
sim <- simulate_ode(model, record = rec, output_every = 1)
```

* `sg_quantity()`, `sg_counter()`, `sg_alleles()` select by stage, sex,
  genotype, mate genotype, node, host state, and group with `by`;
  `sg_quantity_weights(W)` takes any user matrix.
* Built-in sets: `"totals"` (the current block × genotype totals, per
  node; whether it becomes the default is decision 13) in phase 1; others (`"adults_by_node"`,
  `"cases_by_node"`, `"alleles"`) added when a scenario needs them.
* `output_every` (every k days) as now.
* `sg_totals()`, `sg_cases()`, `sg_hosts()` and allele frequencies read
  the recorded quantities; when something they need was not recorded they
  say so (e.g. "record `cases` to use `sg_cases()`").
* `record = "state"` keeps full states (all nodes or selected nodes), for
  tests and small runs.

Memory, 254 counties over 24 years:

| Request | Numbers stored |
|---|---|
| Full state every day | about 260 million (about 2 GB) |
| All counters every day | about 67 million (about 500 MB) |
| Weekly cases per county | about 320,000 (about 2.5 MB) |
| Adults per county every day | about 2.2 million (about 18 MB) |

**Single node, many genotypes.** The recorder does not depend on movement
or on the number of nodes; one node is a landscape of length one. It
matters at one node too when a cube has many genotypes, because mated
females carry the mate genotype and that block grows as the square of the
number of genotypes. States per node with the current defaults:

| Cube | Genotypes | States per node | Of which mated females |
|---|---|---|---|
| `cubeMendelian` | 3 | 207 | 90 |
| `cubeHomingDrive` | 10 | 1,390 | 1,000 |
| `cubeClvR2` | 60 | 38,340 | 36,000 |

Memory, one node over 50 years (about 18,260 days; estimates):

| Request | Homing drive (10) | ClvR2 (60) |
|---|---|---|
| Full state every day (v0.1 behaviour) | about 200 MB per run | about 5.6 GB per run |
| Same, 100 chain-binomial replicates | about 20 GB | about 560 GB |
| `"totals"` (block × genotype) every day | about 9 MB | about 53 MB |
| Allele counts by sex, weekly | well under 1 MB | well under 1 MB |

* Quantities over the mate genotype (e.g. females carrying sterile or
  drive sperm) are single rows of W_x, however large that block is.
* Long runs can be split into chunks (burn-in, then release years) through
  the kept final state, with nothing stored in between.
* The recorder saves memory, not run time: one output costs about one rate
  evaluation, but each evaluation still costs about 17–18 ns per state
  (about 0.7 ms per evaluation for ClvR2 at one node). Speed with many
  genotypes is limited by the female × mate block, which this design does
  not address.

When to build it (decided by the user, 2026-11-02; decision 10): in
phases 1–2 as planned (ODE in phase 1, chain-binomial in phase 2), not as
a stand-alone step before phase 1, although it needs nothing from
movement.

Open point (deferred by the user, 2026-11-04, decision 13: decided once
the recorder is fully integrated in the package; until then v0.1's full
`sim$state` stays the default output):

* **Change from v0.1.** v0.1 returns full states (`sim$state`), read in
  `R/outputs.R`, `R/checks.R`, `R/init-steady.R`, `R/simulate.R`,
  `R/simulate-cb.R` and most tests (18 files). With `"totals"` as the
  default, those uses ask for `record = "state"` explicitly and the
  summary functions read recorded quantities; one-node results are
  unchanged, the default output object is not. Recommendation: one
  default everywhere (`"totals"`), not a default that depends on the size
  of the run, which would surprise users.

### 7.2 Summaries and plots

* `sg_totals()`, `sg_cases()`, `sg_hosts()`, allele frequencies: a `node`
  column; per node or per km²; study nodes by default.
* Counters per node: emigrants, immigrants and losses by class.
* Spread metrics: first-arrival day per node (first day above a threshold),
  number of occupied nodes, front position along a transect.
* Plots: node × time heat map; front position over time.

## 8. Tests

* n_nodes = 1: all v0.1 regressions identical (bitwise for the ODE).
* Two identical nodes, symmetric movement, identical initial states: each
  equals the one-node run.
* Closed edges: the network total is unchanged by movement exactly (ODE to
  round-off; CB as an integer identity). Open edges: the change equals the
  loss counter.
* Genotype, mate genotype and sub-stage preserved by moves.
* Clocks: calendar rates constant; physiological rates scale with Δx below
  T_ref and are zero at and below the development threshold; with
  κ_max = 1 they equal the given rates from T_ref to about 41 °C; with
  κ_max = `Inf` they follow Δx / Δx(T_ref) at all temperatures; an
  intermediate κ_max caps κ at that value. The lifetime moves of the 5.3
  table are reproduced by `check_movement()`.
* R reference vs C++: rates (ODE) and draws (CB), as in
  `test-rates-cpp.R` and `test-backend-equivalence.R`.
* Deterministic CB map vs ODE as dt → 0.
* Transect (built with `sg_movement_diffusion()`): a marked cohort (no
  births, no deaths) spreads with variance 2Dt; front speed at constant
  temperature vs 2√(rD).
* Constructors: `sg_movement_diffusion()` on a chain gives D/h² to each
  neighbour and the closed and open end rates of 6.1; on a regular grid
  δ = 4D/h²; `sg_movement_from_prob()` inverts the matrix exponential to
  round-off; `check_movement()` returns the D used by
  `sg_movement_diffusion()` on a chain and a grid.
* General networks: an irregular network (a few polygons with known
  shared borders) gives rates that match a hand calculation of the
  finite-volume formula.
* Field emigration: with a leave rate from Matlock & Skoda, the share of
  females and males leaving per day is about 16% and 10%.
* Metapopulation (CB): weakly connected nodes go extinct and are
  recolonised; nodes behind a barrier (no edges) stay empty.
* Node extinction rule (ODE): off reproduces v0.1 bitwise; on, a node below
  the threshold at the end of a day is zero the next day, the removal
  counter equals the amount removed, and with closed edges the network
  total changes only by births, deaths and removals.
* Recorder: every recorded quantity equals the same sum computed in R
  from a full-state run (`record = "state"`), for both solvers and every
  replicate; period sums of counters equal the sums of the daily
  counters; `"totals"` matches the v0.1 totals for one node.
* Optional: two-node cross-check against MGDrivE2fast in the
  mosquito-like limit (M6), calendar clock, U rate 0 (MGDrivE2 does not
  move unmated females). The costliest test to set up; the symmetric
  two-node, conservation and diffusion-variance tests above already check
  the movement terms, so it is done only if time allows.

## 8a. Performance

Multi-node runs must finish in as little time as possible. Every work
package keeps to this section; each item names the phase that does it.

**Starting point.** The C++ kernel costs about 17–18 ns per state and
evaluation in every scenario (`bench/README.md`, M3). A 254-county SIT
run has about 30,000 states, so one evaluation takes about 0.5 ms; at
about 20 Dormand–Prince steps a day (6 evaluations each) the ODE takes
about 25 s per simulated year if cost scales linearly with nodes. The
aim is to keep it linear and avoid the traps below.

### 8a.1 Traps that could stop large runs

1. **Dense Jacobian in the stiff fallback** (phase 1).
   `simulate_ode(stiff_fallback = TRUE)` hands stiff days to `lsoda`,
   whose default Jacobian is dense: n² doubles, about 7 GB for 30,000
   states. For more than one node the fallback is off, and asking for it
   gives a clear error. `deSolve::lsodes` with the model's sparsity
   pattern (block diagonal by node plus the movement edges, built once in
   R) is deferred until a multi-node run actually turns out stiff.
2. **Full states every day** (phase 1 for the ODE, phase 2 for the
   chain-binomial). 30,000 states × 365 days × 24 years is about 260
   million doubles (2 GB) per run, more with replicates; all counters
   every day add about 500 MB. The recorder (section 7.1) stores only the
   requested quantities, computed in C++ during the run (e.g. `"totals"`
   per node); the default output is decision 13. The
   `max_memory` check applies to both solvers.

### 8a.2 Savings

3. **Skip empty nodes in the chain-binomial step** (phase 2). A node with
   every compartment at zero and no arrivals in the step is skipped. This
   is exact, and since `RandomDraws::binom()` returns 0 without a call to
   R's generator when n = 0, skipping consumes the same random numbers:
   the draw-for-draw match with the R reference step holds. Before the
   front arrives, and in northern nodes in winter, most nodes are empty.
4. **Forcing stored once per distinct series** (phase 1): e.g. per
   weather station, with a forcing index per node; 254 counties use 170
   stations. Cheap and saves memory. **Sharing the daily preparation**
   (day rates, the chain matrices of S16) across nodes with the same
   forcing is deferred until the phase 3 profile shows it matters: day
   rates are computed once per day, not per evaluation, so for the ODE the
   saving is probably small.
5. **Larval chains** (phase 2). With larval crowding the larval chain is
   rebuilt every step in every node (n_L = 6, small); profile it and
   optimise only if it shows up.
6. **Movement pass fitted to the memory layout** (phase 1). States are
   stored node by node, so a node's adult block of a class is contiguous;
   movement is one loop over the sparse edges, each adding a scaled
   contiguous block, which the compiler can vectorise. R's default
   optimisation flags; no `-march=native` (not portable, not CRAN).
7. **ODE tolerances** (phase 3 benchmark). Nearly empty nodes hold tiny
   values that the error control keeps refining: `atol` set in counts
   (e.g. 1e-6 individuals); test whether `rtol = 1e-6` is accurate enough
   for multi-node runs before changing the default.
8. **One step size for all nodes** (phase 1). Dormand–Prince uses one
   step size for the whole system, so the fastest-changing node sets it.
   Solver diagnostics report which node limits the step.

### 8a.3 Package-writing rules

* R does all setup once; C++ does everything per step. No calls back into
  R during a step (no `Rcpp::Function`, no R loops per day); this is why
  R-function inputs are deferred (decision 8).
* C++ receives flat numeric vectors and compressed sparse row arrays,
  validated once in the `ModelView` constructor; no lists or Rcpp sugar in
  per-step loops.
* The kernel object is reused across replicates (`Rcpp::XPtr`); output
  arrays are allocated once at their known size.
* Parallelism: replicates in separate R processes (`future.apply`,
  L'Ecuyer streams), as now. No C++ threads or OpenMP in v0.2: they would
  break the reproducible random-number rules (D-022) and OpenMP is not
  portable on macOS. Reconsider for the ODE only if benchmarks need it.
* Profile before optimising (Instruments on macOS, `bench::mark()` on
  parts); change only what the profile shows.

### 8a.4 Benchmarks and targets

`bench/m6_multinode_scaling.R` (phase 3; rerun in phase 4): run time of
both solvers against the number of nodes (1, 10, 100, 254) and edges, per
replicate for the chain-binomial solver; results in
`bench/results/m6_multinode_scaling.csv` and `bench/README.md`, like M3
and M4, so speed regressions are visible. Targets:

* linear scaling in the number of nodes (cost per state and evaluation
  stays about 18 ns);
* 254-county SIT (2 genotypes), one year: ODE under about 30 s;
  chain-binomial, one replicate: a target set from the first measurement;
* memory below `max_memory` with the default outputs.

## 9. Implementation plan (work packages)

### Phase 0. Specification (done 2026-11-04)
MODEL_SPEC section 8a (N1–N5) (N1, N2, 5.3, 5.6, the node extinction rule
of 5.7 as an exception to "no population floors"); D-032; parameter-table
entries (`movement_t_ref` = 27 °C, `movement_kappa_max` = 1); this note's
decisions.

### Phase 1. Landscape, one layer, ODE
* `sg_landscape()` with node roles; `sg_movement_layer()`;
  `sg_model(landscape = , movement = )` accepting a list with one layer.
* `ModelView`: node stride; forcing [n_days, n_nodes]; per-node scalars
  (H, N_H, h_M, area); per merged group (5.4; one in phases 1–3): rates
  [3, n_nodes], destinations in compressed sparse row form, loss shares,
  clock; T_ref once for the model.
* `RateKernel`: loop over nodes for day rates and aggregates (wound search,
  mating, larval crowding stay local); one movement pass over adult blocks
  along the edges, with losses; no allocation per evaluation. The rate
  multiplier of a layer for a node and day (the clock now, times the
  schedule in phase 4) is read in one place, so richer time variation can
  be added later without touching the movement pass.
* `.sg_rates_r()`: the same node loop and N1–N2; equivalence tests.
* `sg_releases(node = )`; `sg_init()` per node and burn-in with movement.
* Node extinction rule for the ODE (5.7): end-of-day check per node,
  `node_extinction` argument, removal counter, note in `check_model()`.
* Performance (8a): forcing stored per distinct series with an index per
  node; stiff fallback off for more than one node (an error if asked
  for); `max_memory` check for the ODE; diagnostics report the node that
  limits the step.
* Recorder (7.1) for the ODE: `sg_record()`, `sg_quantity()`,
  `sg_counter()`, `sg_quantity_weights()`; W_x and W_c in compressed
  sparse row form in `ModelView`; the products in the C++ output step;
  final state kept; the `"totals"` set; `record = "state"` (full states
  stay the default until decision 13 is made).

### Phase 2. Chain-binomial movement
Section 5.6 in `.sg_cb_step_r()` and `cb_step.cpp`; counters (emigrants,
immigrants, losses); draw-for-draw and map tests. Performance (8a): skip
empty nodes (exact, same random numbers); profile the per-step larval
chain. Recorder
(7.1) for the chain-binomial solver, per replicate, with the memory
estimate from the recorder's size.

### Phase 3. Constructors, diagnostics, outputs, first transect
* `bench/m6_multinode_scaling.R` and the targets of 8a.4 (on regular
  grids of 1–254 nodes, built with `sg_movement_diffusion()`); the
  tolerance study of 8a item 7; profiling for the deferred items of 8a.2
  (item 4).
* `sg_movement_diffusion()`, `sg_movement_from_prob()` (needed for the
  MGDrivE2 check), conversions. `sg_movement_from_distance()` moves to
  phase 4 with the long-distance layer it mainly serves.
* `check_movement()` and its place in `check_model()`.
* Outputs and plots of section 7: counters summed over periods
  (`every = "week"`), `sg_alleles()`, summary functions reading recorded
  quantities with clear errors, further built-in sets as scenarios need
  them.
* **First landscape (decision 5):** a transect from south Texas north,
  built with `sg_movement_diffusion()` (a few lines of user code, shown in
  the scenario, not a separate function), a source node at the southern
  end with its own weather (section 5.5), weather per node from the
  1961–62 station data (`data-raw/texas_nws/weather_1961_1962/`); run with
  the chain-binomial solver and with the ODE with the node extinction rule
  on (5.7), and once with the rule off to show the atto-fly effect. Checks:
  winter die-back in the north, spring recolonisation, front speed against
  80–160 km per generation with local movement only (expected to be too
  slow, section 2); the same run with all three clocks (calendar, capped
  κ_max = 1, uncapped κ_max = `Inf`), to see whether the clock choice
  matters for die-back and spring spread (section 5.3).
### Phase 4. Several layers
* Lists of layers in `sg_model()`; a schedule per layer (one daily
  multiplier, the same in every node; section 5.2). Mostly R work: layers
  are merged into groups by (clock, schedule) before the run (5.4), so the
  C++ kernel only loops over two or three groups instead of one.
* `sg_movement_from_distance()` (user kernel, cutoff for sparsity).
* A long-distance layer: sparse far edges, calendar clock, monsoon-season
  schedule, northward destinations; the transect rerun with it.
* **Documented example of a general network:** Texas counties as nodes.
  Areas, shared border lengths and centroid distances from county shapes
  (user code with `sf` and `tigris`, in a vignette or `data-raw/`, not a
  package dependency); weather per county from
  `county_station_map_1961_1962.csv` (170 stations, 254 counties, some
  share a station); about 1,500 edges. With D ≈ 0.5 km² d⁻¹ and counties
  of about 2,700 km² the rate to each neighbour is about 2 × 10⁻⁴ d⁻¹, so
  local diffusion hardly moves flies between counties and spread on this
  grid must come from the long-distance layer; `check_movement()` shows
  this as a very small lifetime leave probability. Moved here from phase
  3: without the long-distance layer the example would only show that
  nothing moves.
* `sg_release_barrier()`: barrier width and lane spacing to releases per
  node; barrier-width scenario.
* Texas 1962 county network: `county_station_map_1961_1962.csv` (170
  stations), weekly county cases from `texasScrewwormdat1962.xlsx`; fit
  weeks 1–26, test the sterile releases on weeks 27–52 (extends scenario 20
  from Uvalde to all counties).

### Phase 5 (v0.4). Host transport
Infested-host transport with larvae, quarantine; host-dependent dispersal;
batch events if needed.

Cost: 254 Texas counties with 2 genotypes and n_A = 10 have about 80 adult
states per node; movement costs O(edges × adult states per node), with
edges counted once per merged group (5.4), so long-distance layers are
kept sparse (cutoff). Realistic host numbers (about 1000 km² per node or
more) as for the host model. Performance rules and
targets: section 8a.

## 10. Decisions (made by the user, 2026-11-02; 12–13 and revision of 2 on 2026-11-04)

1. Unmated females (U) move.
2. Movement clock: degree-day scaling with an adjustable cap,
   κ = min(κ_max, Δx(T)/Δx(T_ref)), rates given per day at T_ref = 27 °C;
   `movement_kappa_max` = 1 by default (capped), `Inf` for uncapped;
   calendar clock available per layer. No study supports 27 °C, the cap or
   any of the clocks: the cap avoids claiming more daily movement in heat,
   uncapped keeps lifetime dispersal nearly constant in warm weather, and
   the decline in the cold has qualitative support (section 5.3). The
   clock is therefore a sensitivity dimension: the phase 3 transect runs
   calendar, capped and uncapped. Revised 2026-11-02 from uncapped to
   capped; revised 2026-11-04 to make the cap a parameter (user).
3. Parameterisation: raw rates and destination probabilities as the only
   model input (MGDrivE2 style), with constructors that take a
   user-supplied kernel, exact conversion helpers and `check_movement()`
   (section 3); no built-in kernel.
4. Long-distance movement: in the same rate framework, as additional layers;
   the interface takes a list of layers from the start, one layer is built
   and tested in phases 1–3, several in phase 4. Batch events deferred.
5. First landscape: a transect, built with the general diffusion
   constructor; any network can be given (`sg_movement_layer(rate,
   dest)`), with a Texas county network as the documented general example.
   No transect-only function (user decision, 2026-11-02).
6. Landscape edges: open or closed through π (rows ≤ 1; the remainder is
   lost and counted).
7. Inflow from outside: source nodes with their own weather.
8. Time variation of movement kept minimal: the clock only in phases 1–3;
   one daily multiplier per layer in phase 4; richer forms (node × day
   multipliers, R functions, time-varying destinations, probability time
   series) deferred until a scenario needs them.
9. Outputs: a recorder stores only the quantities requested (weighted sums
   of states and of counters, computed in C++ during the run; ratios in R
   afterwards), with full states available through `record = "state"`
   (section 7.1); the default output is decision 13.
10. The recorder is built in phases 1–2 as planned (ODE in phase 1,
    chain-binomial in phase 2), not earlier as a stand-alone step; it also
    serves one-node runs with many genotypes (section 7.1).
11. Simplifications and additions after review (2026-11-02):
    * an optional node extinction rule for the ODE against atto-flies,
      off by default; spread and die-back questions use it or the
      chain-binomial solver (5.7);
    * layers are merged inside the model by (clock, schedule); the
      chain-binomial step draws emigrants once at the total rate (5.4,
      5.6);
    * one destination matrix for all classes per layer (no
      `dest_by_class`) and one `movement_t_ref` for the model (5.2, 6.1);
    * `sg_movement_from_distance()` and the Texas county example move to
      phase 4;
    * stiff fallback off for more than one node, `lsodes` deferred;
      shared daily preparation deferred until profiled (8a);
    * fixed habitat differences through per-node rates; host-dependent
      dispersal stays a later option (5.8);
    * MGDrivE2fast cross-check optional (section 8).
12. The ODE node extinction rule stays optional and off by default, also
    for multi-node runs (user, 2026-11-04; section 5.7).
13. Whether `"totals"` replaces v0.1's full `sim$state` as the default
    output everywhere is decided once the recorder is fully integrated in
    the package (user, 2026-11-04; section 7.1). Until then full states
    remain the default and `"totals"` is available on request.

## References

* Mayer, D.G. & Atzeni, M.G. (1993) Estimation of dispersal distances for
  *Cochliomyia hominivorax* (Diptera: Calliphoridae). Environmental
  Entomology 22: 368–374. https://academic.oup.com/ee/article/22/2/368/2394476
* Hightower, B.G., Adams, A.L. & Alley, D.A. (1965) Dispersal of released
  irradiated laboratory-reared screw-worm flies. Journal of Economic
  Entomology 58: 373–374.
  https://scilit.com/publications/76df9d02ebf46fa0c36332432e43daeb
* Matlock, R.B. Jr & Skoda, S.R. (2009) Mark–recapture estimates of
  recruitment, survivorship and population growth rate for the screwworm
  fly, *Cochliomyia hominivorax*. Medical and Veterinary Entomology 23
  (Suppl. 1): 111–125.
* Robinson, A.S., Vreysen, M.J.B., Hendrichs, J. & Feldmann, U. (2009)
  Enabling technologies to improve area-wide integrated pest management
  programmes for the control of screwworms. Medical and Veterinary
  Entomology 23 (Suppl. 1): 1–7. (Cites Spradbery 1994 for 290 km,
  Rawlins 1985 for wind-borne reinfestation, Lindquist et al. 1992 for
  Libya.)
* Lyra, M.L., Klaczko, L.B. & Azeredo-Espin, A.M.L. (2009) Complex
  patterns of genetic variability in populations of the New World
  screwworm fly revealed by mitochondrial DNA markers. Medical and
  Veterinary Entomology 23 (Suppl. 1): 32–42. (Cites Torres et al. 2007
  on extinction and recolonisation in Uruguay.)
* McDonagh, L., García, R. & Stevens, J.R. (2009) Phylogenetic analysis
  of New World screwworm fly suggests genetic isolation of some Caribbean
  island populations. Medical and Veterinary Entomology 23 (Suppl. 1):
  14–22.
* FAO. Screwworm control and eradication in the southern United States of
  America. https://www.fao.org/4/u4220t/u4220T0a.htm
* Gutierrez, A.P., Ponti, L. & Arias, P.A. (2019) Deconstructing the
  eradication of new world screwworm in North America. Medical and
  Veterinary Entomology 33: 282–295.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6849717/
* Krafsur, E.S. (1978) Aggregations of male screwworm flies in south Texas.
  https://dr.lib.iastate.edu/entities/publication/35a60af5-83b9-4aec-9c4c-0a09c50c503c
* Behavior of sexually active male screwworm flies.
  https://scilit.net/publications/6209773e2c8e202ab726c23d0d6bfc6a
* Thomas, D.B. et al. (1992) Persistence of sterile screwworm flies at a
  release site. https://pubmed.ncbi.nlm.nih.gov/1593015
* A mathematical model for the efficient control of the New World
  screwworm (2026). PLOS ONE.
  https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0355963
* USDA (2025) Mexico confirms case of New World screwworm in Nuevo León.
  https://www.usda.gov/about-usda/news/press-releases/2025/09/21/mexico-confirms-case-new-world-screwworm-nuevo-leon
* Wu, S.L. et al. (2021) MGDrivE 2. PLOS Computational Biology 17:
  e1009030 (SI section 1.1.5, movement).
