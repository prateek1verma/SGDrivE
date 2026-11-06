# Phase 2 implementation plan: movement in the chain-binomial solver

Status: **accepted**, 2026-11-05 (choices of section 4 confirmed by the user
the same day); not implemented. Implements phase 2 of dev/MOVEMENT_DESIGN.md
(section 9) against MODEL_SPEC section 8a (N4) and D-032. Phase 1 (ODE,
recorder) is in dev/MOVEMENT_PHASE1_PLAN.md.

## 1. Goal and ground rules

`simulate_cb()` runs landscapes with adult movement in whole numbers of
flies (MODEL_SPEC N4): local extinction and recolonisation without
fractions of flies, and the recorder of phase 1 for stochastic runs.

* One-node runs stay draw-for-draw identical to the code before phase 2
  (same seed, same random numbers, same results).
* Every change goes into the R reference step (`R/cb-reference.R`) and the
  C++ step (`src/cb_step.h`, `src/cb_step.cpp`) together, with the
  draw-for-draw test (CONVENTIONS, random numbers).
* No calls into R during a step; no new dependencies; replicates stay in
  separate R processes for parallel runs.
* The node extinction rule (N5) stays ODE-only: counts are whole numbers.

## 2. Starting point

* The C++ stepper (`CbStepper`) works on one node: `prepare_day()` builds
  the chains of the day (egg, pupa, adult and infestation chains; the larval
  chain every step with crowding), `step()` draws in the order aggregates
  (S15), eggs (S17), hosts (H3), mating (S18), chain moves (S16), emergence
  (S19), arrivals.
* `.sg_cb_step_r()` follows the same order; with the same seed the two give
  identical results (`test-cb-cpp.R`).
* `simulate_cb()` refuses landscapes (`.sg_one_node_only()`), stores
  releases per block for one node, and offers `keep = "totals"` instead of a
  recorder; totals are computed in R after each replicate.

## 3. Work packages

### WP1. Specification and baseline (S)

* MODEL_SPEC section 10 ("Order of draws in a step") and N4: with several
  nodes a step has three parts:
  1. for each node in index order: aggregates (S15), eggs (S17), hosts
     (H3), mating (S18), chain moves (S16);
  2. movement (N4) of all nodes: emigrants per adult compartment in state
     order (node slowest), destinations and lost; immigrants added after
     all draws;
  3. for each node in index order: emergence (S19) and arrivals.
  With one node and no movement this is the current order. D-032 records
  it.
* Single-node chain-binomial baseline: seeded stochastic runs (`keep =
  "state"`) and deterministic-map runs of three configurations (default
  Mendelian with releases, gp2014 homing, host model with crowding and an
  Allee effect), saved from `main` before any change
  (`data-raw/make_cb_baseline.R`, fixture `cb_single_node_baseline.rds`,
  test `test-cb-single-node-baseline.R`): bitwise with
  `SCREWDRIVE_BITWISE_BASELINE=true`, tolerances elsewhere.
* *Done (2026-11-05):* MODEL_SPEC section 10 and N4 updated, D-032 note.
  The baseline compares stochastic runs and steps draw for draw only in the
  bitwise mode: on another platform a last-bit difference in a probability
  can change one draw and then the whole path; elsewhere it compares the
  deterministic map and the deterministic step.

### WP2. Stepper per node, without movement (L)

* C++: `step()` runs parts 1 and 3 per node on the node's part of the state
  (`node_stride` offsets) with the node's forcing series, wounds, hosts and
  mate finding (`RateKernel::set_node()`). A small work area per node holds
  what part 1 hands to part 3 (eggs by genotype, chain outflows, matings,
  deaths, senescence, host flows, fertile batches).
* Chains are prepared once per day per distinct forcing series and shared
  by the nodes on that series (exact); with larval crowding the larval
  chain is rebuilt per node and step.
* Skip empty nodes: a node without flies skips the fly parts and still does
  its host draws. Draws of zero individuals use no random numbers, so the
  result is the same; only time is saved.
* Flows per node in the counter layout of the ODE (per-node blocks, "@id").
* R reference: the same three parts per node (`.sg_cb_step_r()` over the
  node setups of `.sg_reference_setup()`).
* Tests: one node bitwise equal to the WP1 baseline; C++ = R draw for draw
  on 3 nodes without movement (static and hosts, Mendelian and homing); in
  the deterministic map each node equals its own one-node model.
* *Done (2026-11-05):* the R step is split into `.sg_cb_part1()` and
  `.sg_cb_part3()` per node; the C++ stepper has `SeriesChains` per forcing
  series and a `NodeWork` area per node; the day loop's extinction checks
  span all nodes. One-node runs are identical draw for draw to the WP1
  baseline; one-node `simulate_cb()` is not slower (about 15-20% faster in
  an A/B run, both versions compiled alike).

### WP3. The movement step N4 (M)

* R and C++: for each adult compartment, E ~ Binomial(X, 1 - exp(-lambda dt));
  destinations over the union of the groups' edges in increasing node order,
  then lost, as sequential conditional binomials; destination shares
  computed once a day per node and class from the groups' rates and the
  clock kappa; immigrants added after all draws (no fly moves twice in a
  step). Emigrant, immigrant and lost counters per node and class (whole
  numbers). Deterministic map: expectations.
* Tests: C++ = R draw for draw for 3 nodes x {static, hosts} x {Mendelian,
  homing} x {no movement, physiological, calendar, open edges}; R and C++
  deterministic maps agree; closed edges: emigrants = immigrants and the
  network total unchanged exactly; open edges: change = lost; movers keep
  genotype, mate and sub-stage; one step from one node on a ring 1 -> 2 -> 3
  leaves node 3 empty; statistics of emigrant counts and destination splits
  (binomial and multinomial means and variances).
* *Done (2026-11-05):* `.sg_cb_move_r()` and `.sg_cb_move_shares()` in R,
  `CbStepper::prepare_movement()` (once a day, with `prepare_day()`) and
  `CbStepper::movement()` in C++; the union of the groups' edges per node
  is built once in the constructor. C++ = R draw for draw in all 16
  configurations and for two groups put together by hand (one layer per
  model is accepted for now). Models without movement draw no extra random
  numbers, so the WP1 baseline and the WP2 tests are unchanged.

### WP4. `simulate_cb()` for landscapes and multi-node outputs (M)

* Guard removed; releases per stage and node ("M@north", as in the ODE);
  `extinct_day` and `wild_extinct_day` for the whole landscape; early stop
  when every node is empty and no release is left; `print()` for several
  nodes; memory estimate.
* Outputs on replicate arrays: `sg_cases()` and `sg_flow_balance()` per node
  on [time, column, replicate] arrays; `sg_totals()` and the other summaries
  through the existing node selection.
* Tests: summaries against direct sums of states and counters per
  replicate; flow balance per node with residual exactly 0; parallel
  replicates; print.

### WP5. Recorder for `simulate_cb()` (M)

* The C++ recorder (`Csr`, `Recorder`) moves from `simulate.cpp` into a
  shared header and is used by the chain-binomial day loop: weighted state
  sums on output days, counter sums every day, final state, per replicate.
* `simulate_cb(record = )` with the meaning of `simulate_ode(record = )`;
  `keep` keeps working as a shortcut (`keep = "totals"` is `record =
  "totals"`); giving both is an error. `sim$recorded` [time, quantity,
  replicate], `sim$counter_sums` [day, quantity, replicate],
  `sim$final_state` [state, replicate]; `max_memory` from the recorder's
  size times the replicates.
* Tests: with the same seed, recorded sums equal the sums of a full-state
  run per replicate (recording uses no random numbers); totals equal the
  former `keep = "totals"`; counter sums; final states; memory limit.

### WP6. Performance, validation and documentation (M)

* Speed: A/B timing of one-node `simulate_cb()` against the code before
  phase 2 (within 5%); time against 1, 10 and 50 nodes (linear growth,
  saving from skipped empty nodes); profile of the per-step larval chain
  with crowding (design note 8a item 5), optimised only if it shows up.
* Validation: the deterministic map converges to the multi-node ODE as
  dt -> 0; the replicate mean converges to the deterministic map;
  metapopulation: weakly connected nodes go extinct and are recolonised,
  nodes behind a barrier (no edges) stay empty.
* Documentation: help of `simulate_cb()`, NEWS, D-032 findings, design note
  and roadmap status; `R CMD check` 0/0/0; sanitizers.

### Pull requests

* PR D: WP1–WP3 (specification and baseline, stepper per node, movement
  step).
* PR E: WP4–WP6 (`simulate_cb()` for landscapes and outputs, recorder,
  performance, validation, documentation).

Size: large, about 900–1,200 lines with tests, mostly the C++ restructuring
of WP2.

## 4. Choices confirmed by the user (2026-11-05)

1. Order of draws across nodes: the three parts of WP1 (as N4 reads); not a
   complete step per node with arrivals set aside until after movement.
2. A single-node chain-binomial baseline is saved before any change.
3. `simulate_cb()` gains `record` with the meaning of the ODE; `keep` keeps
   working as a shortcut; both together are an error.
4. Extinction days stay landscape-wide; per-node first-arrival and
   extinction days are the phase 3 spread metrics.
5. Prepared chains are shared by the nodes on the same forcing series.
6. Two pull requests: D (WP1–WP3) and E (WP4–WP6).
