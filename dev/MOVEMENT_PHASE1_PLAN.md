# Phase 1 implementation plan: landscape, one movement layer, ODE

Status: **accepted**, 2026-11-04 (choices of section 8 confirmed by the
user the same day). Implements phase 1 of
dev/MOVEMENT_DESIGN.md (section 9) against MODEL_SPEC section 8a (N1–N3,
N5) and D-032. No code has been written yet. The chain-binomial step (N4)
is phase 2; constructors, `check_movement()`, period sums of counters and
the first transect are phase 3; several layers and schedules are phase 4.

## 1. Goal and ground rules

At the end of phase 1 a user can build a landscape of N nodes, each with
its own weather, area, wound or host density and releases. They can give
one movement layer (leave rates by class and node, a destination matrix,
a physiological or calendar clock), run `simulate_ode()` with the C++ or
R backend, optionally switch on the node extinction rule, store full
states, per-node totals or chosen quantities, and read per-node results.

Ground rules, from CONVENTIONS and the design note:

* **One node reproduces v0.1 bitwise** (ODE): same states, counters,
  state and counter names, model hash, print output and snapshots, apart
  from the two new parameters in the parameter table.
* Every model change goes into the R reference (`.sg_rates_r()`) and the
  C++ kernel together, with an equivalence test.
* R validates and flattens once; C++ does the per-step work, allocates
  nothing per evaluation and never calls back into R.
* No new package dependencies (Imports stay Rcpp, stats, utils; no
  `Matrix`, no `sf`).
* No slowdown of single-node runs: `bench/m3_cpp_vs_reference.R` within
  5% of the current numbers.

Features not in phase 1 stop with a clear error that names the phase:
`simulate_cb()` on a model with more than one node or with movement
(phase 2), more than one movement layer or a schedule (phase 4),
`calibrate_population()` and `sg_model(adult_females = )` with a landscape
(later; both assume one node), `check_movement()` (phase 3).

## 2. Starting point (what the code does now)

* **State layout** (`.sg_layout()`, R/model.R): one vector of blocks E, L,
  P, M, U, F (then HS, HW, HI with hosts). `dims$n_node = 1` is multiplied
  into each block's length, so the node is the slowest index *inside* each
  block (block-major). MODEL_SPEC 8a asks for node-major order (all blocks
  of node 1, then node 2). With one node the two are identical.
* **Model view** (R/model-view.R, src/model_view.h/.cpp): scalars and
  vectors for one node: `H`, `h_M`, `N_H`, daily forcing vectors of length
  `n_days`, block offsets `offsets[6]`, `host_offsets[3]`.
* **Rate kernel** (src/rates.cpp): `day_rates(day)`, `aggregate()` and
  `eval()` work on one node through the block offsets; counters follow
  `.sg_counter_names()`. `max_rate()` covers one node.
* **ODE driver** (src/simulate.cpp, R/simulate.R): day loop with releases
  at the start of each day, Dormand–Prince with counters integrated
  alongside, `lsoda` stiff fallback, full states stored on output days.
* **Outputs** read `sim$state` (ODE) or `sim$state` / `sim$totals`
  (`simulate_cb(keep = "totals")`) through `.sg_sim_array()`, and parse
  block and genotype from state names ("M[AA,1]", "F[AA|aa,3]").
* **Parameters** are validated as finite numbers (`.sg_check_value()`),
  so `movement_kappa_max = Inf` needs a change there.
* **Snapshots** (tests/testthat/_snaps) cover cube, forcing, init, params
  and releases; the params snapshot will change with the two new
  parameters.

## 3. Design choices for phase 1

Choices marked **(confirmed)** were recommendations confirmed by the user
on 2026-11-04 (section 8).

### 3.1 Node-major state layout

The state vector becomes the concatenation of per-node vectors. Each
per-node vector has today's layout (E, L, P, M, U, F, then HS, HW, HI), so
every current offset is an offset within a node and node i starts at
(i − 1) × S, with S the number of states per node.

* `model$layout` keeps one row per block with offsets **within a node**;
  new fields `model$n_states_node` (S) and `model$n_states` (N × S).
* `sg_state_index(model, block, genotype, mate, substage, node = NULL)`:
  new `node` argument (ids or numbers; default all), result ordered with
  the node slowest.
* The model view passes `node_stride = S`; C++ adds `node * node_stride`
  to every offset. With one node nothing changes.

Why node-major and not block-major: the per-node kernel code can be reused
unchanged on a pointer `y + node * S`, a node's adult block of each class
is contiguous (the movement pass adds scaled contiguous blocks), and
skipping empty nodes in phase 2 is a range check.

### 3.2 State and counter names (confirmed)

* One node: names unchanged ("M[AA,1]", "eggs_laid[AA]").
* Several nodes: the node id is appended after "@": "M[AA,1]@n3",
  "HW[hosts]@n3", "eggs_laid[AA]@n3". Every current parser reads the block
  from the first character and the genotype from the first bracket, so they
  keep working; the node is `sub(".*@", "", name)`.

### 3.3 Landscape object

`sg_landscape(nodes, forcing)`:

* `nodes`: data frame with `id` (character, unique, no "@"), `area`
  (km², > 0) and optional `x`, `y` (km, for later diagnostics), `role`
  (`"study"` default or `"source"`), `wound_density`, `host_density` (per
  km²; missing means the value in `sg_params()`).
* `forcing`: one `sg_forcing` for all nodes, or a list with one per node
  (names = node ids). All must have the same dates and the same forcing
  parameters (as `.sg_check_forcing_params()` now checks for one).
* Forcing is stored once per distinct series (identical objects are found
  by a hash) with an index per node, so 254 counties on 170 stations store
  170 series (design note 8a.2 item 4).
* Default when no landscape is given: one study node with `area` from
  `sg_model(area = )`; that path is today's code.

### 3.4 Movement layer object

`sg_movement_layer(rate, dest, clock = "physiological")`:

* `rate`: list with `male`, `virgin`, `mated`, each a scalar or one value
  per node (d⁻¹ at T_ref), all ≥ 0 and finite.
* `dest` **(confirmed)**: either a numeric N × N matrix (row = from, column
  = to) or a data frame of edges `from`, `to`, `prob` (node ids). Stored as
  a sorted edge list; zero entries dropped. Checks: π_ii = 0, entries in
  [0, 1], row sums ≤ 1 + 1e-12, ids known. No `Matrix` dependency.
* `clock`: `"physiological"` or `"calendar"`.
* `sg_model(movement = )` takes one layer or a list of layers; phase 1
  accepts exactly one (an error for more names phase 4). The merge into
  groups (MODEL_SPEC 8a) is written in general form now, as it is short and
  makes phase 4 an R-only change.

### 3.5 Per-node quantities in the model

`sg_model(params, cube, forcing, releases, landscape = NULL, movement =
NULL, ...)`. With a landscape, `forcing` and `area` come from it (giving
`area` too is an error). New model fields:

* `model$landscape`: nodes table, forcing series and node → series index.
* `model$counts`: `H`, `h_M`, `N_H` as vectors of length N (H_i = wound
  density_i × A_i, h_M A_i, N_H,i = host density_i × A_i).
* `model$movement`: merged groups with per-class rates [3, N], edges
  (CSR: `row_ptr`, `col`, per-class probabilities [3, edges]), loss shares
  [3, N], clock; plus `dx_ref` = `sg_degree_days(movement_t_ref, params)`
  (S1 with the floor) and `kappa_max`.
* `model$meta$hash` gains the landscape and movement **only when they are
  not the default**, so single-node hashes, `sg_state` checks and
  snapshots are unchanged.

### 3.6 Counters

* Per-node counters: today's list for each node, in node order (names as
  in 3.2).
* Movement counters, present when the model has movement, per node and
  class: `emigrants_M`, `emigrants_U`, `emigrants_F`, `immigrants_*`,
  `lost_*` (rates in the ODE, N2).
* Node extinction counters `removed_E` … `removed_F` per node, present only
  when `simulate_ode(node_extinction = TRUE)`. They are jumps added by the
  day loop to that day's row, not integrated rates. Infested hosts returned
  to S are added to that node's `infestations_ended` (D-032).
* `sg_flow_balance()` adds immigrants − emigrants − lost and − removed to
  each node's adult (or stage) balance.

### 3.7 Outputs and the recorder

Decision 13 keeps full states as the default output, so
`simulate_ode(record = "state")` is the default and nothing changes for
existing code.

* `record = "totals"`: block × genotype totals per node (and host blocks),
  in `sim$totals` with the same column naming as
  `simulate_cb(keep = "totals")`, so `.sg_sim_array()` and the summary
  functions already handle it.
* `record = sg_record(...)`: chosen quantities in `sim$recorded`
  (matrix [output day, quantity]) and, if `sg_counter()` terms are given,
  chosen counter combinations in `sim$counters`; otherwise all counters
  daily as now.
* `sim$final_state` is always kept (decision 9), so a run can be continued
  whatever was recorded.
* Summary functions that need something not recorded stop with a message
  naming what to record. Reading arbitrary recorded quantities in
  `sg_totals()` etc. is phase 3.

### 3.8 Releases and initial states

* `sg_releases(..., node = NULL)`: node ids; **(confirmed)** required when
  the model has more than one node (an error otherwise names the ids), and
  ignored for one node. Release sizes per km² use A_i.
* `sg_init(..., nodes = NULL)`: the nodes that receive the initial
  population for `"constant"` and `"hardy_weinberg"` (default all nodes;
  others start empty), with `per = "node"` or `"km2"` as now (A_i per
  node). `"state"` takes a full vector of length N × S. `"burn_in"` runs
  the multi-node ODE with movement. `"equilibrium"` and `"threshold"` are
  computed per node from that node's forcing (a constant temperature per
  series) without movement (D-032); `temperature` may be one value or one
  per node.

## 4. Work packages

Sizes are rough: S under 150 lines of code and tests, M 150–500, L over
500.

### WP1. Parameters (S)

Files: R/param-table.R, R/params.R, tests/testthat/test-params.R,
test-param-table.R, the params snapshot.

* Rows `movement_t_ref` (27, degC, `(dev_temp_min, 45]`, spec N3, symbol
  T_ref) and `movement_kappa_max` (1, "-", `(0, Inf]`, spec N3, symbol
  kappa_max) in a new group `"movement"`.
* `.prow(allow_inf = FALSE)` and `.sg_check_value()` accepting `Inf` only
  where `allow_inf = TRUE` (only `movement_kappa_max`).
* Cross-check: `movement_t_ref` above `dev_temp_min`, so Δx(T_ref) > 0
  with the default floor.
* Not forcing parameters (they do not change `sg_forcing()`).
* Done when: parameter tests pass; the params snapshot is updated in a
  separate commit so the change is visible.

### WP2. Landscape and movement layer objects (M)

Files: new R/landscape.R, R/movement.R; tests test-landscape.R,
test-movement-layer.R.

* `sg_landscape()` and `sg_movement_layer()` with validation, `print()`
  and `as.data.frame()` methods (3.3, 3.4).
* Internal `.sg_movement_groups(layers, node_ids)`: merge by (clock,
  schedule) into per-class rates, CSR edges, per-class probabilities and
  loss shares (MODEL_SPEC 8a "Layers and groups"); exact for one layer.
* Internal `.sg_forcing_series(forcings)`: distinct series and index.
* Tests: every validation error (dimensions, π_ii ≠ 0, rows > 1, negative
  or non-finite rates, unknown ids, mismatched forcing dates or
  parameters); merge of two layers by hand (tested now, used in phase 4);
  CSR layout; series deduplication.

### WP3. Multi-node model and node-major layout (M)

Files: R/model.R, R/releases.R, tests test-model.R, test-releases.R, new
test-model-nodes.R.

* `sg_model(landscape = , movement = )` (3.5); per-node `counts`; host
  setup per node (`n` a vector); `dims$n_node = N`.
* `.sg_layout()` per node; `n_states_node`; `sg_state_index(node = )`;
  `sg_state_names()` with the "@id" suffix for N > 1.
* Releases per node (3.8): `sg_releases(node = )`, `c()`, `print()`,
  `.sg_resolve_releases()` with per-node area, `.sg_release_events()`
  with the node offset.
* `print.sg_model()`: node count, study and source nodes, forcing series,
  movement summary (rates by class, edges, open-edge share, clock,
  T_ref, κ_max).
* Hash rule of 3.5.
* Errors for `adult_females` with a landscape and for more than one layer.
* Tests: one-node model identical to today's (`identical()` on layout,
  names, counts, releases, hash); index arithmetic for 3 nodes with hosts
  and an F block; names; releases routed to the right node.

### WP4. Model view and C++ kernel per node (L)

Files: R/model-view.R, src/model_view.h/.cpp, src/rates.h/.cpp,
tests test-model-view.R, test-rates-cpp.R, new test-rates-nodes.R.

Model view additions (validated in the `ModelView` constructor):

* `n_node`, `node_stride`, `n_series`, `forcing_index` [N];
  forcing fields `dx`, `dx_egg`, `mu_T`, `mu_E`, `phi_ov`, `moisture` as
  [n_days, n_series] column-major (one column for one node: same values,
  same memory order as today).
* Per-node `H`, `h_M`, `N_H` [N] (`h_M_param` stays the switch).
* Movement: `n_groups`; per group the clock, rates [3, N], CSR
  `row_ptr` [N + 1], `col`, probabilities [3, edges], loss shares [3, N];
  `dx_ref`, `kappa_max` (Inf passed as is). Empty when there is no
  movement.
* Counter count = N × per-node counters + movement counters.

Kernel changes:

* `DayRates` computed per node from the node's forcing series.
  *Implemented (2026-11-05):* recomputed at each evaluation, like the
  movement clock κ, instead of cached per day: a few operations per node
  against the node's aggregates, so a cache would save little; revisit if
  the phase 3 profile shows it.
* `eval()` loops over nodes and runs today's body on `y + i·S`,
  `dydt + i·S`, `counters + i·C` with node i's `H`, `h_M`, `N_H` and day
  rates (the body becomes `eval_node()`, with the same floating-point
  operations in the same order). `Aggregates` stays one work struct
  reused per node; phase 2 will keep one per node.
* **Movement pass (N1–N3)** after all nodes: for each group and node, κ
  for the day (min(κ_max, dx_i / dx_ref), or 1 for the calendar clock),
  λ_s = κ δ_s,i; for each class block (M: n_g·n_A, U: n_g·n_A, F:
  n_g²·n_A states, contiguous) subtract λ x_i from node i and add λ π x_i
  to each destination along the CSR row; accumulate emigrant, immigrant
  and lost counters. κ and λ are computed once per day with the day rates.
  No pass when there is no movement.
* `max_rate()`: maximum over nodes, plus λ for adults when movement is on.
* R reference: `.sg_reference_setup()` gains per-node quantities;
  `.sg_rates_r()` loops nodes over today's body (moved into
  `.sg_rates_node()`) and then applies `.sg_movement_r()`, written
  directly from N1–N2 with matrices so it reads like the specification.
* Tests: R = C++ (derivatives and counters to `sg_tol_exact` relative) for
  3 nodes × {static, hosts} × {Mendelian, homing} × {no movement,
  physiological, calendar, open edges}; one-node results `identical()` to
  the current kernel output saved as a fixture before the change; movement
  terms sum to zero over nodes with closed edges (per compartment, to
  round-off); with open edges the network change equals the loss counter;
  genotype, mate genotype and sub-stage preserved (moving a single F
  compartment); clock values (calendar 1; physiological = Δx ratio below
  T_ref, 0 at and below T_L with floor 0, small positive with floor 0.01;
  capped at 1 with κ_max = 1, Δx ratio with `Inf`, capped at 1.3 with 1.3).

### WP5. ODE driver for several nodes (M)

Files: src/simulate.cpp, src/dopri5.h, R/simulate.R, tests
test-simulate-cpp.R, test-simulate.R, new test-movement-ode.R.

* `simulate_ode(..., stiff_fallback = NULL, max_memory = 1e9)`:
  `stiff_fallback` defaults to `TRUE` for one node and `FALSE` for more;
  `TRUE` with more than one node is an error that explains why (dense
  `lsoda` Jacobian). Other deSolve methods stay available; for `lsoda`,
  `bdf` and `radau` with more than one node, an error when the dense
  Jacobian (n² × 8 bytes) would exceed `max_memory`.
* `max_memory` check before the run: output days × recorded columns plus
  days × counters, × 8 bytes (as `simulate_cb()`).
* Diagnostics: new column `limiting_node` (NA for one node). Dopri5 records,
  at each accepted step, the component with the largest scaled error; the
  day loop reports the node that held it most often that day. Cost: one
  argmax per step, done inside the existing error-norm loop.
* Releases per node are already state positions (WP3), so the day loop is
  unchanged.
* Tests: two identical nodes with symmetric movement and identical states
  each equal the one-node run (to `sg_tol_ode`); a marked cohort with no
  births and no deaths on a closed chain keeps its total and spreads; flies
  appear in an empty node only through movement; source node feeds a study
  node; stiff fallback error; memory error; C++ vs R backend (`ode45`) on a
  3-node run to `sg_tol_ode`.

### WP6. Node extinction rule (S–M)

Files: src/simulate.cpp, R/simulate.R, tests new test-node-extinction.R.

* `simulate_ode(node_extinction = FALSE)`. When `TRUE`, at the end of each
  day (after the day is integrated, before the next day's releases):
  for each node whose fly total is below `extinction_threshold`, set its
  fly states to zero, add the removed amounts by stage to `removed_*`, and
  with hosts move Σ HI to HS and add it to `infestations_ended` (N5).
  Implemented in the C++ day loop and in the R deSolve loop.
* Tests: `FALSE` is bitwise the current run; `TRUE` zeroes a dying node the
  next day and the counters equal the amounts removed; host conservation
  (S + W + I constant) holds; with closed edges the network total changes
  only by births, deaths, removals; works with one node too; the "atto-fly"
  demonstration (a node far from the source fills with fractions of a fly
  without the rule and stays empty with it).

### WP7. Initial states per node (M)

Files: R/model.R (helper), R/init.R, R/init-steady.R, R/checks.R
(`sg_r0()` help); tests test-init.R, test-init-steady.R, new
test-init-nodes.R. Only R code.

Approach: every initial-state method is already correct for one node, so a
helper gives node i as a one-node model and each method runs per node; the
results are placed into the full state vector (node slowest).

1. `.sg_node_model(model, i)`: the model reduced to node i (its forcing
   series, area, H, h_M, host count; no landscape, movement or releases),
   made by replacing those fields so the numbers are exactly the
   multi-node model's; `.sg_place_node(x, x_node, i)`. Tested against a
   one-node `sg_model()` built from scratch for each node.
2. `"constant"` and `"hardy_weinberg"`: new `nodes` argument (ids or
   numbers; default all nodes). Selected nodes get the one-node result
   (`per = "km2"` scales by each node's area); other nodes get no flies;
   every node gets its hosts at their fly-free steady state.
3. `"equilibrium"` and `"threshold"`: per node from its one-node model,
   without movement (D-032). `temperature` one value or one per node;
   default each node's constant forcing temperature (an error naming the
   node whose forcing is not constant). Nodes without an equilibrium start
   empty, with one message listing them. The attribute is a per-node table
   (node, temperature, R0, eggs per day, adult females).
4. `sg_equilibria()`: guard removed; with several nodes one data frame with
   a `node` column (attribute `R0` named by node).
5. `"burn_in"`: the multi-node ODE with movement from the `"constant"`
   start (so `nodes` sets where flies start). Without `weather` each node
   repeats the first year of its own forcing series; `weather` is one data
   frame for all nodes or a list per node (as in `sg_landscape()`). The
   start-effect check stays.
6. `"state"`: a full-length vector, or a named list of per-node vectors.
7. `print.sg_state()`: with several nodes, stage totals by node, hosts per
   node and per-node equilibrium information.
8. Clean-up: guards in `sg_init()` and `sg_equilibria()` removed; help of
   `sg_r0()` notes that R0 at low density does not depend on a node's area
   or wounds; `calibrate_population()` stays one-node only.
9. Tests: one-node results unchanged (existing tests and the baseline);
   helper against `sg_model()`; constant and Hardy-Weinberg per node, `nodes`
   and `per = "km2"`, hosts in every node, unknown nodes; equilibrium and
   threshold per node against one-node models (different temperatures,
   temperature vector, hosts, a cold node starting empty); `sg_equilibria()`
   rows per node; burn-in with movement fills an empty connected node, two
   identical symmetric nodes equal the one-node burn-in, each node repeats
   its own series; `simulate_ode(m)` runs with the default start; print.
10. Checks: full suite, bitwise baseline, lint with the package installed,
    sanitizers.

*Implemented (2026-11-05).* Found on the way: with the host-wound model, an
equilibrium below the oviposition threshold (e.g. 13 degC) stopped with an
error in the one-node code too (fertile batches per egg became 0 / 0); it now
gives an empty state with the usual message (regression test in
test-init-steady.R). Only R code changed, so the sanitizers were not rerun.

Choices confirmed by the user (2026-11-05):

* `sg_init(node_extinction = FALSE)`, passed only to the burn-in run, so a
  burn-in with movement can start without fractions of flies in every node.
* `sg_equilibria()` with several nodes returns one data frame with a `node`
  column.
* Burn-in `weather`: one data frame for all nodes or a list per node.
* `"state"` also accepts a named list of per-node vectors.
* Default start for landscapes (option A of four discussed): `sg_init(m)`
  and the default start of `simulate_ode(m)` put the Gutierrez & Ponti
  starting state (counts per node) into every node, study and source alike,
  as for one node; `nodes` and `per = "km2"` restrict or scale it. Rejected:
  per km^2 by default, source nodes only, and an error when
  `simulate_ode()` has no `init` on a landscape.

### WP8. Summaries, checks and flow balance per node (M)

Files: R/outputs.R, R/checks.R, R/plots.R, R/simulate.R (print), tests
test-outputs.R, test-checks.R.

* `sg_totals()` and `as.data.frame()`: `"node"` allowed in `by`; a
  `nodes` argument (default: study nodes) for which nodes to include;
  `per = "km2"` divides by the summed area of those nodes.
* `sg_females()`, `sg_cases()`, `sg_hosts()`, `sg_allele_freq()`: sum
  over the chosen nodes (default study nodes) with a `nodes` argument;
  per-node columns are phase 3 (design note 7.2).
* `sg_flow_balance()`: per node, with movement and removal counters.
* `check_model()`: per-node table (area, wounds or hosts, forcing series),
  movement summary (largest λ at T_ref per class, open-edge share), and a
  note for more than one node that ODE runs without `node_extinction` can
  show recolonisation by fractions of a fly.
* `print.sg_sim()`: totals over study nodes and the number of occupied
  nodes at the end.
* Plots: existing plots use the study-node sums; heat maps and front
  plots are phase 3.
* Tests: per-node sums against direct sums of state columns; study vs
  source selection; flow balance closes to `sg_tol_exact` per node with
  movement, open edges and extinction removals.
* *Implemented (2026-11-05):* `nodes` argument (default study nodes) in
  `sg_totals()`, `sg_females()`, `sg_allele_freq()`, `sg_cases()`,
  `sg_hosts()`, the plots and `sg_plot_data()`; `as.data.frame()` gives all
  nodes by node. For multi-node runs `sim$releases` has one column per stage
  and node (`"M@north"`), which `sg_flow_balance()` needs for per-node
  balances. `check_model()` checks each node as a one-node model and runs
  the population-scale check on the whole landscape with movement. The flow
  balance closes to about 1e-12 relative.

### WP9. Recorder for the ODE (L)

Files: new R/record.R; src/simulate.cpp, R/simulate.R (day loops), R/outputs.R
(readers), R/simulate-cb.R (shared totals matrix); tests new
test-recorder.R. The chain-binomial recorder is phase 2; counter sums by
period, `sg_alleles()`, further built-in sets and summaries that read
arbitrary recorded quantities are phase 3.

Goal: a run stores only what is asked for, computed in C++ during the run;
full states stay the default (decision 13), so existing code is unchanged.

```r
rec <- sg_record(
  "totals",
  sg_quantity("adults", stage = c("M", "U", "F"), by = "node"),
  sg_quantity("wild_F", stage = "F", genotype = "AA", node = c("a", "b")),
  sg_counter("cases", counter = "fertile_batches", by = "node")
)
sim <- simulate_ode(model, record = rec)
```

1. `sg_quantity(name, stage, sex, genotype, mate, substage, node, host,
   by)`: selects states (`stage` E, L, P, M, U, F or `"adults"`,
   `"juveniles"`; `host` HS, HW, HI); `by` (`"node"`, `"genotype"`,
   `"stage"`, `"mate"`) splits it into one quantity per group, named in the
   style of state names, e.g. `adults[AA]@north`.
2. `sg_counter(name, counter, node, by)`: sums of daily counters by base
   name (e.g. `"fertile_batches"`, `c("deaths_M", "senescence_M")`); unknown
   names give an error listing those available (including `removed_*` with
   the node extinction rule).
3. `sg_quantity_weights(W, names)`: any weights, as a dense matrix
   [quantities x states] or a table (`quantity`, `state`, `weight`); covers
   allele copies and other weighted sums.
4. `sg_record(...)`: built-in sets `"state"`, `"totals"` and `"counters"`
   with quantities and counter sums. `record = "state"` (default) keeps
   `sim$state`; `"totals"` gives `sim$totals` (block x genotype per node and
   host blocks, the column names of `simulate_cb(keep = "totals")` with
   "@id" for several nodes); quantities go to `sim$recorded` [output day
   (with day 0), quantity]; counter sums replace the daily counters in
   `sim$counters` unless `"counters"` is also requested.
   `sim$final_state` is always kept.
5. Internals: R builds W_x [quantities x states] and W_c [counter sums x
   counters] in CSR form (no `Matrix` dependency).
6. C++ day loop: on output days y = W_x x (and the state copied when
   `"state"` is recorded); every day W_c c when counter sums are requested.
   The recorder is an argument of the run, not part of the model view, and
   its arrays are validated in C++.
7. R paths (deSolve loop, single-node fallback day): the same products from
   the same arrays.
8. `max_memory` from the recorder's size (quantities x output days, counter
   columns x days, final state).
9. Readers: `.sg_sim_array()` reads ODE totals; summaries work with
   `"state"` or `"totals"` (what needs the full state, e.g. grouping by mate,
   gives the existing message or `NA`, as for chain-binomial totals);
   with a custom record without either they stop with a hint such as
   `record "totals" (or the state) to use sg_totals()`. `sg_cases()`,
   `sg_flow_balance()` and the counter lines of `print()` need all counters
   and stop with a hint to add `"counters"` otherwise.
10. `print.sg_sim()` with only recorded quantities: final values of the
    quantities and the solver lines.
11. Tests: every recorded quantity equals the same weighted sum of a
    `record = "state"` run (both backends, several genotypes, hosts, three
    nodes with movement); counter sums equal sums of the full counters
    (with `removed_*`); `"totals"` equals totals from the state and uses the
    chain-binomial column names; summaries agree between `"totals"` and
    `"state"`; the single-node baseline is unchanged; two runs continued
    from `final_state` equal one long run; `by` splitting and names;
    `sg_quantity_weights()` in both forms; errors for unknown stages,
    genotypes, nodes and counters; memory limit; hint messages.
12. Checks: full suite, bitwise baseline, lint with the package installed,
    sanitizers (the C++ day loop changes), a quick A/B timing (the default
    `record = "state"` is not slower; a small record is at least as fast).

Choices confirmed by the user (2026-11-05):

1. The recorder is an argument of the run, not part of the model view.
2. Quantity names follow the state names: `name[group]@node`.
3. With `sg_counter()` terms only those counter sums are kept; adding
   `"counters"` keeps all counters too; without them all counters are kept.
4. Summaries with a custom record (no `"state"` or `"totals"`) stop with a
   "record ... to use ..." message; reading recorded quantities is phase 3.
5. `sg_quantity_weights()` takes a dense matrix or a table (`quantity`,
   `state`, `weight`), without a `Matrix` dependency.
6. Counter sums by period stay in phase 3 (weekly totals can be computed
   from daily counter sums in R).

*Implemented (2026-11-05).* Two details settled while coding: counter sums
are stored in their own `sim$counter_sums` (and `sim$counters` is `NULL` when
only sums are kept), so both can be kept side by side; and `node = NULL` in
`sg_quantity()` and `sg_counter()` means the study nodes, as in the summary
functions. Totals are W_x rows too, so the C++ loop has one product for all
recorded sums.

### WP10. Documentation and housekeeping (S–M)

* roxygen for every new exported function (purpose, arguments with units,
  value, runnable example on a 3-node landscape); `NAMESPACE` via
  roxygen; `_pkgdown.yml` reference section "Landscapes and movement".
* NEWS entry; README line.
* D-032 "Findings during implementation" (as D-030 does) for anything that
  changes in the specification while coding; MODEL_SPEC updated in the
  same PR if so.
* MOVEMENT_DESIGN.md: phase 1 marked done.

## 5. Order, dependencies and pull requests

```
WP1 ─┐
WP2 ─┼─> WP3 ─> WP4 ─> WP5 ─┬─> WP6
     │                      ├─> WP7
     │                      ├─> WP8
     │                      └─> WP9 ─> WP10
```

WP4 is the core and the largest risk; WP6–WP9 depend on it but not on each
other. **(confirmed)** Three pull requests, each green on CI:

1. **PR A** (WP1–WP5): parameters, objects, multi-node model, kernel with
   movement, ODE driver. Single-node bitwise test and R = C++ tests.
2. **PR B** (WP6–WP8): extinction rule, initial states, summaries and
   checks per node.
3. **PR C** (WP9–WP10): recorder, documentation.

## 6. Verification gates (every PR)

* `devtools::test()`: all current tests pass unchanged, except the params
  snapshot (WP1), plus the new tests.
* Single-node bitwise check: a fixture of states and counters from
  `main` (gp2014 preset, McAllen, Mendelian and homing cubes, static and
  host mode) saved before WP3; `identical()` after every PR.
  *Found on CI (2026-11-05):* bitwise equality holds only on the machine and
  R version that made the fixture (CI differs by about 1e-15 relative, and
  the hash differs on Windows), so the test is bitwise with
  `SCREWDRIVE_BITWISE_BASELINE=true` (`dev/check.R`) and uses tolerances of
  1e-9 or tighter elsewhere.
* `lintr::lint_package()` clean (with the package installed so
  `object_usage_linter` resolves), `Rscript dev/style.R`, clang-format.
* `bash dev/sanitize.sh` after src/ changes (WP4, WP5, WP6, WP9).
* `devtools::check()` 0 errors, 0 warnings, 0 notes; coverage of new code
  ≥ 85% (`dev/coverage.R`).
* `bench/m3_cpp_vs_reference.R`: single-node C++ time within 5% of the
  stored result; a quick 1, 10, 50-node timing of the ODE (closed chain)
  to confirm cost grows linearly (the full M6 benchmark is phase 3).

## 7. Risks

* **Bitwise identity for one node.** Moving the kernel body into a node
  loop must not reorder floating-point operations. Mitigation: the
  fixture test from `main` from the first commit of WP4.
* **Forcing as matrices.** Every reader of `forcing$dx` etc. in R (init,
  calibration, checks) must use the node's series. Mitigation: one helper
  `.sg_node_forcing(model, i)` returning a one-node forcing; grep for all
  `forcing$` uses in WP3.
* **Name parsing.** Regexes in outputs, checks and plots that assume no
  suffix. Mitigation: tests with "@id" names in WP8.
* **ODE step size.** One step size for all nodes; a node with fast
  dynamics slows all. The `limiting_node` diagnostic shows which.
* **Scope.** The recorder (WP9) is the largest single piece and is
  independent of movement; it can slip to its own PR without blocking
  phase 2.

## 8. Choices confirmed by the user (2026-11-04)

1. State and counter names carry an "@id" suffix for more than one node;
   one-node names are unchanged (3.2).
2. The destination matrix is given as a dense matrix or an edge table
   (`from`, `to`, `prob`), without a `Matrix` dependency (3.4).
3. `node` is required in `sg_releases()` when the model has more than one
   node (3.8).
4. Three pull requests: A (WP1–WP5), B (WP6–WP8), C (WP9–WP10)
   (section 5).
