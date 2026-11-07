# Phase 4 implementation plan: several layers and long-distance movement

Status: **draft**, 2026-11-06, for the user's review (choices in section
4); not implemented. Implements phase 4 of dev/MOVEMENT_DESIGN.md (section
9). Phases 1-3 are in dev/MOVEMENT_PHASE1_PLAN.md to
dev/MOVEMENT_PHASE3_PLAN.md; findings in D-032.

## 1. Goal and ground rules

Give the model the second scale of movement that the transect showed is
missing: local diffusion moves the front about 30 km per generation, the
1960s spread was 80-160 km per generation (D-032 phase 3 findings; design
note section 2). Long-distance movement (wind, long flights of gravid
females) comes as additional movement layers with their own clock and a
seasonal schedule.

* Several layers per model, merged into groups by (clock, schedule) before
  the run (design note 5.4); the kernels loop over two or three groups.
* Time variation stays minimal (decision 8): one daily multiplier per
  layer, the same in every node.
* One-node results and every existing one-layer landscape stay identical
  (bitwise baselines; a one-layer model without a schedule draws the same
  random numbers and gives the same ODE results).
* No new package dependencies: county shapes are prepared by user code in
  `data-raw/` (with `sf`, `tigris`), and the package ships at most a small
  derived table.

## 2. Starting point

* `sg_movement_layer(rate, dest, clock)` builds one layer; `sg_model()`
  stops for a list of several. `.sg_movement_groups()` already merges
  layers by clock (rates add, destinations weighted by rate), and both
  kernels (ODE in C++ and R, chain-binomial in C++ and R) already loop over
  groups and union their edges; tests put two groups together by hand.
* No schedules: a group's rate varies in time only through its clock.
* Constructors: `sg_movement_diffusion()` (local), `sg_movement_from_prob()`;
  no distance-kernel constructor.
* Landscapes: the transect (scenarios 21-22) and a 50-county demonstration
  with cells approximated from centroids (branch `texas-50-counties-video`,
  not merged). County centroids and stations for 253 counties are in
  `data-raw/texas_nws/weather_1961_1962/county_station_map_1961_1962.csv`
  (the row named "Denton" is DeWitt County); weekly 1962 cases per county
  are in `texasScrewwormdat1962.xlsx` (Wk_1..Wk_52).

## 3. Work packages

### WP1. Specification (S)

* MODEL_SPEC 8a: the layer schedule m^l(t) >= 0 (daily, one value per
  layer and day, the same in every node) in N1-N4: lambda = sum over groups
  of kappa^G m^G delta^G; groups by (clock, schedule). D-032 note.

### WP2. Several layers and schedules (M)

* `sg_movement_layer(..., schedule = )` (choice 1) and `sg_model(movement =
  list(...))`; groups by (clock, schedule), the schedule aligned with the
  forcing's dates (an error when they do not cover it).
* Model view: `group_schedule` [n_days, n_groups]; the ODE kernel, the
  chain-binomial day preparation and both R references multiply kappa by
  the day's multiplier; `check_model()`, `check_movement()` and `print()`
  report each group (lifetime numbers on days when the layer is active).
* Tests: a one-layer model without a schedule is unchanged (baselines); a
  schedule of all ones equals none; zero days move nothing; two layers equal
  the hand-merged groups of phase 2; C++ = R (rates, draws); the network
  total is conserved with closed edges.

### WP3. `sg_movement_from_distance()` (M)

* `sg_movement_from_distance(dist, kernel, rate, cutoff = 1e-6, clock,
  schedule, direction = NULL)`: a user function of distance gives relative
  weights; weights below `cutoff` times the row's largest are dropped
  (sparsity); rows normalised as choice 2 says; an optional `direction`
  (a bearing in degrees with a spread, or a user function of the bearing)
  for directed transport such as monsoon winds (choice 4).
* Tests: a kernel that is 1 for nearest neighbours reproduces a chain;
  cutoffs drop the right edges; directed weights only point the right way;
  row sums and losses as specified.

### WP4. Long-distance layer on the transect (M)

* `inst/scenarios/23_texas_transect_long_distance.R`: scenario 21 plus a
  long-distance layer (choice 4: northward fat-tailed kernel, calendar
  clock, monsoon-season schedule); a sweep over its rate and reach to find
  what gives 80-160 km per generation; both solvers; the front, winter
  die-back and spring recolonisation as in scenario 21.
* Findings in D-032: which long-distance rates are consistent with the
  observed spread, and how they compare with the 290 km recoveries of
  marked flies (Hightower et al. 1965).

### WP5. Texas county network (L)

* `data-raw/texas_counties/` (choice 3): user code that builds, for the 254
  Texas counties, areas, shared border lengths, centroid distances and
  coordinates, saved as small CSV files; county weather from the 170
  stations of 1961-1962 (as scenario 21; the "Denton"/DeWitt row
  corrected in the derived table).
* A documented example (scenario 24 or a vignette, choice 6): the 254
  counties with local diffusion and the long-distance layer of WP4;
  `check_movement()` shows how little local diffusion moves between
  counties (design note phase 4).
* Tests: the derived table is consistent (symmetric borders, areas summing
  to Texas's area within a tolerance, every county with at least one
  neighbour except where expected).

### WP6. `sg_release_barrier()` (M)

* Releases per node for a sterile-fly barrier: a band (coordinates, width)
  and lane spacing to release events per node, in proportion to the part
  of each node inside the band (choice 5); a barrier-width scenario on the
  transect with the long-distance layer.
* Tests: total released equals band area times density; nodes outside the
  band get none; per-node shares as specified.

### WP7. Performance and documentation (S)

* Benchmarks with several layers (long-distance edges are the cost:
  design note phase 5 note); help, NEWS, D-032 phase 4 findings, design
  note and roadmap status, `R CMD check`, sanitizers.

### Pull requests

* PR I: WP1-WP4 (layers, schedules, the distance constructor, the
  long-distance transect).
* PR J: WP5-WP7 (county network, barrier helper, documentation).

Size: large, about 1,500 lines with tests and scenarios; C++ only for the
schedule multiplier in both kernels.

## 4. Choices for the user

1. **Schedule input.** A numeric vector of daily multipliers aligned with
   the forcing's dates, or also a helper `sg_schedule(months = 7:9)` that
   builds one from months (and optionally values per month).
   *Recommended: both* (the vector is exact; the helper covers the common
   seasonal case).
2. **Losses in `sg_movement_from_distance()`.** Rows always normalised to 1
   (closed; flies only land on nodes), or open rows with the kernel mass
   that falls outside the nodes lost, which needs the kernel as a 2-D
   density and a numerical integral per node. *Recommended: closed by
   default, with an optional `lost` share per node given by the user*
   (simple and exact; open rows from kernel mass can follow if needed).
3. **County geometry.** Download the Census TIGER/Line county shapes once
   (installing `sf` and `tigris` on this computer; needs your OK for the
   downloads) and commit only the derived CSV tables; or keep the
   centroid-cell approximation of the 50-county script. *Recommended:
   download once* (real areas and borders; the approximation misplaces
   borders by tens of km).
4. **The long-distance layer.** Northward fat-tailed (Cauchy-like) kernel
   on the calendar clock with a monsoon-season schedule (July-September,
   Gutierrez, Ponti & Arias 2019), or an undirected fat-tailed kernel all
   year, or both as variants. *Recommended: both*, the directed seasonal
   one as the main case (it is the hypothesis in the literature), the
   undirected year-round one as the null.
5. **Barrier releases.** Shares of each node inside the band from node
   coordinates and areas only (approximate: a node counts in full or by the
   overlap of its equivalent square), or from the county shapes of WP5
   where available. *Recommended: the shapes where given, the
   approximation otherwise.*
6. **The county example.** A scenario script (runs in a few minutes, like
   scenario 22) or a vignette built with the package (must run on CRAN-like
   checks, so needs shipped data and a short run). *Recommended: a scenario
   script*; a vignette can follow once the data question is settled.
7. **The 1962 county calibration** (fit weeks 1-26, test the sterile
   releases on weeks 27-52; extends scenario 20). In phase 4, or as its
   own phase after it, once the county network exists; it also needs the
   1962 sterile release schedule by county, which is not in the repository.
   *Recommended: its own phase*, starting with the release data question.
8. **Pull requests.** Two (I, J as above) or three. *Recommended: two.*
