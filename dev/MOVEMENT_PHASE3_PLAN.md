# Phase 3 implementation plan: constructors, diagnostics, outputs, first transect

Status: **accepted**, 2026-11-05 (choices of section 4 confirmed by the user
the same day); not implemented. Implements phase 3 of dev/MOVEMENT_DESIGN.md (section
9) with the interface of sections 6 and 7 and the tests of section 8.
Phases 1 (ODE, recorder) and 2 (chain-binomial movement, recorder) are in
dev/MOVEMENT_PHASE1_PLAN.md and dev/MOVEMENT_PHASE2_PLAN.md; findings in
D-032.

## 1. Goal and ground rules

Make landscapes practical to build, check and read, and run the first
landscape: a transect from south Texas north.

* Users build movement from biology (a diffusion coefficient, a daily
  leave fraction, a transition-probability matrix) instead of raw rates on
  a particular grid, and `check_movement()` translates rates back into
  biology, so a rate copied to the wrong grid shows up (design note,
  section 3).
* Long runs of many nodes store only what is needed: counters summed over
  weeks, and summaries that read recorded quantities.
* Spread is measured directly: first-arrival days, occupied nodes, front
  position.
* The transect answers the open questions of phases 1–2: does the clock
  choice matter for winter die-back and spring spread (decision 2), how
  large is the atto-fly effect without the node extinction rule, and does
  the chain-binomial step order (emerged adults do not move in their
  emergence step; D-032 phase 2 findings) change the front at dt = 0.1.
* No change to the model equations (MODEL_SPEC N1–N5) unless WP8 shows the
  step-order effect matters (then a separate, specified change).
* One-node results stay identical (bitwise baselines); no new
  dependencies (`sf` and `tigris` stay out; ggplot2 stays in Suggests).

## 2. Starting point

* Movement is given only through `sg_movement_layer(rate, dest, clock)`
  (rates per class at T_ref, destinations as a matrix or an edge table);
  one layer per model. `sg_landscape()` stores optional coordinates `x`,
  `y` (km), not yet used.
* `check_model()` on a landscape reports per-node rates, R0 and wounds and,
  for movement, the clock, number of edges, largest rate per class and the
  nodes with open edges; nothing about lifetime dispersal.
* Recorder (both solvers): `sg_record()` with `"state"`, `"totals"`,
  `"counters"`, `sg_quantity()`, `sg_quantity_weights()` and
  `sg_counter()` (daily sums). Summaries need full states or totals
  (`sg_totals()`, `sg_females()`, `sg_allele_freq()`, `sg_hosts()`) or all
  counters (`sg_cases()`, `sg_flow_balance()`) and say so when they are
  missing. Full states remain the default output (decision 13 open).
* No spread metrics, no node × time plots.
* Weather for the transect: `data-raw/texas_nws/weather_1961_1962/` (170
  GHCN stations, 253 county centroids with their nearest station, gap
  filled, 1961–1962), already used by scenario 20 (Uvalde). Not shipped
  with the package; scenarios read it from the repository.

## 3. Work packages

### WP1. `sg_movement_diffusion()` (M)

* Finite-volume local diffusion on any network (design note 6.1):
  rate i → j = D · L_ij / (d_ij · A_i) [d⁻¹ at T_ref], with `area` A_i
  (km²), shared border lengths L_ij and centroid distances d_ij (km), and
  `D` per class (male, virgin, mated; one value for all or one per class).
  Returns an ordinary `sg_movement_layer` (δ_i = Σ_j rate_ij,
  π_ij = rate_ij / δ_i), so users can inspect and change it.
* Closed edges by default; open edges with `outer_border` (each node's
  border with the outside) and `outer_dist` (default twice the
  centroid-to-border distance); the outward rate becomes the lost share.
* Inputs as matrices or as an edge table (choice 1); the package computes
  nothing from shapes.
* Tests: on a chain of spacing h each neighbour gets D/h², an interior
  node leaves at 2D/h², a closed end at D/h², an open end at 2D/h² with
  half lost; on a regular grid δ = 4D/h²; an irregular network of a few
  polygons with known borders matches a hand calculation; a marked cohort
  on a long chain (no births or deaths, generator exponentiated in the
  test) spreads with variance 2Dt.
* *Done (2026-11-05):* `sg_movement_diffusion(area, border, dist, D,
  outer_border, outer_dist, clock)` in `R/movement.R`; returns an
  ordinary layer (edge table with node ids, `n_nodes` set). A pair in an
  edge table may be listed once or both ways (equal values); the default
  outer distance is `area / outer_border`. Tests in
  `test-movement-diffusion.R`.

### WP2. `sg_movement_from_prob()` and conversions (S)

* `sg_movement_from_prob(P)`: a daily transition-probability matrix to a
  layer, δ_i = −ln P_ii, π_ij = P_ij / (1 − P_ii) (rows of P summing below
  1 become losses), as MGDrivE2 `movement_prob2rate()`. Exact for the
  daily leave probability; the destination split ignores several moves in
  a day, as in MGDrivE2 (documented).
* Conversions (choice 2): daily emigration fraction f → δ = −ln(1 − f)
  (Matlock & Skoda 2009: 16% of females, 10% of males per day); lifetime
  leave probability → δ, the exact inverse of WP3's calculation for the
  adult chain at T_ref.
* Tests: round trips (rate → fraction → rate; rate → lifetime probability
  → rate) to round-off; `sg_movement_from_prob()` reproduces the daily
  leave probabilities and destination shares of P; Matlock & Skoda
  fractions give the expected rates.
* *Done (2026-11-05):* `R/movement-tools.R`: `sg_movement_from_prob()`,
  `sg_rate_from_fraction()`, `sg_rate_from_lifetime(p, model, class, clock,
  genotype)`; internal `.sg_adult_life()`, `.sg_lifetime_leave()` and
  `.sg_expected_life()` (exact recursions over the adult chain, for WP3).
  The lifespans reproduce the design note table of 5.3 to its printed
  digits; `sg_movement_from_prob()` equals MGDrivE2 `movement_prob2rate()`
  on closed rows. Tests in `test-movement-tools.R`.

### WP3. `check_movement()` and `check_model()` (M)

* For each node and class at T_ref (design note 6.2): daily leave
  probability 1 − e^(−δ); lifetime leave probability 1 − E[e^(−δτ)], with
  τ the adult lifespan of the wild type's adult chain and hazard at T_ref,
  computed exactly from the chain (no simulation); expected moves per
  lifetime δ E[τ]; with coordinates, expected path length per lifetime and
  implied D = δ Σ_j π_ij d_ij² / (2 · dimension) (dimension 1 for nodes
  with collinear neighbours, else 2); share of movers lost.
* Optional table over constant temperatures with the clock and κ_max in
  use, reproducing the table of design note 5.3 (calendar, capped,
  uncapped).
* In `check_model()` (summary and warnings): lifetime leave probability
  above 0.99; κδ dt above 0.5 for the chain-binomial step (with
  `check_model(dt = )`); the largest-rate check includes κδ.
* Tests: the 5.3 table to its printed precision; `check_movement()`
  returns the D given to `sg_movement_diffusion()` on a chain and a grid;
  the Matlock & Skoda rates give 16% and 10% daily leave probability;
  warnings appear for an implausible rate (a 1 km-scale rate on a 50 km
  grid) and for a large κδ dt.
* *Done (2026-11-06):* `R/movement-check.R`: `check_movement(model,
  temperatures, genotype)` with `.sg_movement_shares()` (the shares of the
  chain-binomial step, for any kappa per group and node; tested against
  `.sg_cb_move_shares()` with two groups) and `.sg_clock_kappa()`; the
  implied D uses dimension 1 when all nodes lie on a line. In
  `check_model()`: `movement$check`, lifetime leave and D ranges by class,
  `nodes$max_move_rate` (kappa at the node's warmest day), the largest
  step probability, and the two warnings. Tests in `test-movement-check.R`.

### WP4. Counter sums over periods and a built-in cases set (M)

* `sg_counter(every = "day")` with `"week"`, `"month"`, `"year"`
  (choice 3): the recorder adds each day's W_c c into the row of its
  period, in C++ for both solvers, and stores one row per period (the
  memory saving that makes weekly county cases over decades cheap: about
  2.5 MB for 254 counties over 24 years instead of about 500 MB of daily
  counters).
* A built-in set `"cases"`: the counter sums `sg_cases()` needs (fertile
  batches or infestations, egg batches, superinfestations), per node,
  daily or by period; `sg_cases()` reads it when full counters were not
  kept.
* Tests: period sums equal the sums of the daily counters of a full run,
  for both solvers and every replicate (ODE to round-off, chain-binomial
  exactly); partial last periods; `sg_cases()` from the set equals
  `sg_cases()` from full counters; memory estimates count period rows.
* *Done (2026-11-06):* `every` is an argument of `sg_record()`, not of
  `sg_counter()`: a run has one period for all its counter sums (and the
  `"cases"` set), so `sim$counter_sums` stays one matrix (array for
  `simulate_cb()`) with `sim$counter_periods`. `Recorder::add_sums()` in
  `src/recorder.h` adds a day into its period's row (`c_period` from R);
  the ODE driver adds C++ chunks and the R paths (deSolve backend, stiff
  fallback days) into period rows, so a period may span them. The
  `"cases"` set names its sums like the counters (`fertile_batches@a`) and
  covers every node; `sg_cases()` reads it through `.sg_cases_source()`.
  Tests in `test-recorder-periods.R`.

### WP5. Allele quantities and summaries reading recorded quantities (M)

* `sg_alleles(name, alleles = NULL, stage = "adults", sex = NULL, by =)`:
  allele copies as recorded quantities (one W_x row per allele and group;
  default alleles split from the genotype names as in `sg_allele_freq()`);
  `sg_allele_freq()` reads them when neither states nor totals were kept.
* Summaries look for what they need in the order states, totals, recorded
  quantities, and name the missing piece otherwise (e.g. "record
  sg_alleles() or "totals" to use sg_allele_freq()").
* Decision 13 (choice 4, decided): the default output of both solvers
  changes here from full states to `"totals"` plus all counters
  (`sim$final_state` always kept); code and tests that need states ask for
  `record = "state"` (`keep = "state"` stays for `simulate_cb()`). Totals
  keep their current columns (mated females summed over mates), so by
  default `sg_females()$fertile_equivalent` is `NA` and
  `sg_totals(by = "mate")` needs `record = "state"`; help pages and NEWS
  say so.
* Tests: every summary gives the same result from states, from totals and
  from the recorded quantities it accepts (both solvers, 1 and 3 nodes);
  clear errors otherwise.
* *Done (2026-11-06):* `sg_alleles(name, alleles, stage, sex, node, by)`
  (rows weighted by allele copies, `.sg_allele_rows()`; described in
  `meta$record$alleles`) and `.sg_allele_freq_recorded()`. Default
  `record = "totals"` in both solvers; `check_model()` reads totals and
  burn-ins read `final_state`. The tests that read `sim$state` ask for
  `record = "state"` (added with a parser-based rewrite to 150 calls;
  `sim_ode()` in `helper-backends.R` defaults to it). README prose updated
  in `README.Rmd` and `README.md` (not re-knitted: no pandoc locally).
  Tests in `test-summaries-recorded.R`.

### WP6. Spread metrics and plots (M)

* `sg_spread(sim, threshold = , nodes = , position = )` (choice 5): per
  node the first-arrival day (first day with at least `threshold` flies),
  the last occupied day and the days occupied; over time the number of
  occupied nodes and the front position (largest `position` of an
  occupied node; default the nodes' `x` coordinate). Per replicate for the
  chain-binomial solver. Works from states, totals or a recorded
  per-node quantity.
* Plots: `plot(sim, what = "nodes")` (node × time heat map of flies or
  adults per km², study nodes by default) and `what = "front"` (front
  position over time, replicates as thin lines).
* Tests: metrics against direct calculations on small runs (both
  solvers); a front moving through a chain at a known order of arrival;
  plot data shapes (`sg_plot_data()`), as the existing plot tests do.
* *Done (2026-11-06):* `R/spread.R`: `sg_spread(sim, threshold, nodes,
  position, quantity)` with a print method; the simulation's node
  information now carries `x`, `y` and `extinction_threshold`. Plots
  `"nodes"` (viridis heat map, a panel per replicate, the first replicate
  by default) and `"front"` in `R/plots.R`. Tests in `test-spread.R`.

### WP7. Benchmarks, tolerance study, profiles (M)

* `bench/m6_multinode_scaling.R` (design note 8a.4): both solvers on
  regular grids of 1, 10, 100 and 254 nodes built with
  `sg_movement_diffusion()`, SIT (2 genotypes), one year; per replicate
  for the chain-binomial solver; results in
  `bench/results/m6_multinode_scaling.csv` and `bench/README.md`.
  Targets: linear growth in nodes (about 18 ns per state and evaluation
  for the ODE); 254-node SIT ODE year under about 30 s; a chain-binomial
  target set from this first measurement.
* Tolerance study (8a item 7): `rtol` 1e-6 vs 1e-8 and `atol` in counts on
  multi-node runs, accuracy against speed; the defaults change only if the
  study supports it (a separate decision).
* Profiles (Instruments, `bench::mark()`): sharing the daily preparation
  of the ODE across nodes on one forcing series (8a.2 item 4); the
  chain-binomial larval chain under crowding (D-032 phase 2 findings).
  Optimise only what the profile shows, keeping the draw-for-draw match and
  the baselines.

### WP8. First landscape: the south Texas transect (L)

* `inst/scenarios/21_texas_transect.R` (user code, no transect function;
  decision 5), layout of choice 6: cells in a line from the Rio Grande
  valley north, built with `sg_movement_diffusion()`; weather per cell
  from the 1961–1962 stations; a source cell at the southern end with Rio
  Grande valley weather (design note 5.5); D from Matlock & Skoda (about
  0.4 km² d⁻¹ for females, less for males), with a ×10 sensitivity run.
* Runs (choice 7): 1961 as a burn-in year, 1962 reported; the
  chain-binomial solver (replicates at dt = 0.1 and 0.05), the ODE with
  the node extinction rule, and the ODE without it (atto-fly effect); all
  three clocks (calendar, capped κ_max = 1, uncapped κ_max = `Inf`).
* Checks, written up in the scenario and D-032: winter die-back in the
  north after the 1961–62 winter; spring recolonisation; front speed
  against 80–160 km per generation with local movement only (expected to
  be far too slow; design note section 2); whether the clock choice
  changes die-back and spring spread; the front at dt = 0.1 vs 0.05 (if
  the difference matters, a specified fix for movement of emerged adults
  follows as its own change; D-032 phase 2 findings).
* Figures in `inst/scenarios/scenario21_figures/` (heat maps, fronts).

### WP9. Documentation and release checks (S)

* Help pages and examples for the new functions; `_pkgdown.yml`
  reference; NEWS; D-032 phase 3 findings; design note and roadmap status;
  `R CMD check` 0/0/0; sanitizers; full matrix on the pull requests.

### Optional WP10. Cross-check against MGDrivE2 (M)

Left optional by the user (choice 9); not scheduled. Two forms:
* movement only: two nodes of adults with death and movement but no births
  or ageing, calendar clock, mating at emergence, so both packages solve
  the same linear system; also `sg_movement_from_prob()` against MGDrivE2
  `movement_prob2rate()`. Small and exact;
* full mosquito-like limit (milestone M6): a one-node match of the life
  cycle first (constant temperature, instant mating, exponential adult
  life, no wound limit, MGDrivE2's larval density dependence; some terms
  may not map exactly), then two nodes. MGDrivE2 is installed;
  MGDrivE2fast is not.
The conservation, symmetric-node and diffusion-variance tests already
check the movement terms (design note section 8).

### Pull requests

* PR F: WP1–WP3 (constructors, conversions, `check_movement()`).
* PR G: WP4–WP6 (period counter sums, allele quantities, summaries from
  recorded quantities, spread metrics, plots).
* PR H: WP7–WP9 (benchmarks, transect, documentation), and WP10 if done.

Size: large, about 1,500–2,000 lines with tests and the scenario, mostly R;
C++ only for period sums (WP4) and whatever the profiles of WP7 justify.

## 4. Choices confirmed by the user (2026-11-05)

Choices 1–8 as recommended; for choice 4 the plain `"totals"` (not a
variant keeping the mate genotype of mated females, which was discussed);
choice 9: WP10 stays optional.

1. **Inputs of `sg_movement_diffusion()`.** Dense matrices (`area`,
   `border`, `dist`) as in the design note, or also an edge table
   (`from`, `to`, `border`, `dist`), as `sg_movement_layer()` accepts
   both. *Recommended: both* (an edge table is what users get from
   `sf::st_touches()` and border lengths, and it stays sparse for large
   networks).
2. **Conversion helpers.** Two small functions,
   `sg_rate_from_fraction(f)` and `sg_rate_from_lifetime(p, model)`, or
   one `sg_movement_rate(fraction = , lifetime = )`. *Recommended: two
   functions* (each does one thing, and the lifetime one needs a model).
3. **Period sums of counters.** Computed in C++ with only period rows
   stored (saves memory for long runs), or computed in R from daily sums
   after the run (simpler, no saving). And the week: 7-day blocks from the
   first simulated day (as `sg_cases(period = "week")` now), or calendar
   weeks. *Recommended: C++ with period rows; 7-day blocks from the first
   day* (runs starting on 1 January match the 1962 weekly data).
4. **Decision 13, the default output.** Switch the default of both solvers
   from full states to `"totals"` plus all counters in WP5 (the design
   note's recommendation: one default everywhere, small outputs for big
   runs; full states with `record = "state"`), or keep full states as the
   default and decide at the v0.2 release. *Recommended: switch in WP5*
   (the recorder is now in both solvers and the summaries read totals; it
   changes the default output object, not any result).
5. **Spread metrics.** A node counts as occupied when all its flies (E, L,
   P, M, U, F) total at least `threshold`, default `extinction_threshold`
   (1 individual), as the node extinction rule; or adults only. Front
   position from the nodes' `x` coordinate unless `position` is given.
   *Recommended: all flies, threshold 1, `x` by default.*
6. **Transect layout.** Regular cells (for example 18 cells of 50 km × 50
   km from about 26°N to 34°N along about 98°W, each with the weather of
   the station nearest its centre), or a chain of real counties (irregular
   spacing; county areas would have to be added to the data). *Recommended:
   regular cells* (clean D/h², known areas of 2,500 km², realistic host
   numbers; the county network is phase 4).
7. **Transect runs.** Chain-binomial with 20 replicates at dt = 0.1 and
   0.05, ODE with and without the node extinction rule, three clocks, D
   and 10 × D. *Recommended as listed* (about 30 runs, a few minutes each
   at most).
8. **Pull requests.** Three (F, G, H as above) or two (constructors and
   diagnostics; outputs, benchmarks and transect). *Recommended: three.*
9. **MGDrivE2 cross-check (WP10).** Include it in phase 3 or leave it
   optional. *Recommended: optional*, done only if time allows.
