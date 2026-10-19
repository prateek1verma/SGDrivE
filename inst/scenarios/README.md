# Scenario scripts

Short, self-contained R scripts that run simple scenarios with SGDrivE, written
as a user would write them. They serve to test the model scientifically and as
material for tutorials and vignettes.

After installing the package, the scripts are in

```r
system.file("scenarios", package = "SGDrivE")
```

A script can be run line by line (for example after opening it in RStudio
with `file.edit()`), or all at once with `echo = TRUE`, which prints results
and plots as an interactive session would:

```r
file.edit(system.file("scenarios", "01_wild_type_dynamics.R", package = "SGDrivE"))
source(system.file("scenarios", "01_wild_type_dynamics.R", package = "SGDrivE"), echo = TRUE)
```

## Scripts

| Script | Scenario |
|---|---|
| `01_wild_type_dynamics.R` | Wild-type screwworm population at McAllen, Texas (1971-1972): model check, burn-in, deterministic and stochastic runs, plots of stages, females and cases, R0 by temperature |
| `02_constant_temperature.R` | Wild-type population at a constant 20 degC in a region sized for 1000 adult females (`sg_model(adult_females = )`), started at equilibrium: deterministic and stochastic runs |
| `03_allee_threshold.R` | Allee threshold at 20 degC for several half-saturation male densities h_M: populations starting just below the threshold (`sg_init(method = "threshold")`) die out, those just above it grow to the equilibrium |
| `04_host_limited_cases.R` | Host-wound model at McAllen (1971) in a 1000 km^2 region with 30,000 animals: animals newly infested compared with egg batches holding viable eggs, and infested animals in deterministic and stochastic runs |
| `05_larval_crowding.R` | Constant 25 degC, 1000 km^2: growth from a small population with and without larval crowding (`larval_crowding_mode = "per_wound"`); both grow alike at first, crowding levels the population off lower (`sg_equilibria()`) |
| `06_rainfall_pupae.R` | Constant 25 degC, 1000 km^2: wet (600 mm per year, no stress) and dry (300 mm per year) climates with drought mortality of pupae (`pupa_drought_mortality`); from the wet equilibrium the dry population falls to a lower level |
| `07_classic_sit.R` | Classic SIT at constant 25 degC as in mgdrive-nws analysis 01: 26 weekly releases of sterile males (eta = 0.4) and sterile females into 20,000 wild adults, two-genotype W/S cube, stochastic replicates of wild-type adult females |
| `08_classic_sit_mcallen_cycle.R` | The classic SIT releases of scenario 07 under the mean annual temperature cycle of McAllen (each calendar day averaged over 1960-1983), two years from 1 January: region sized for a yearly mean of 20,000 adults, burn-in start, temperature cycle and wild-type females with and without releases |

## Conventions

* One scenario per script, numbered in the order of a tutorial
  (`01_...`, `02_...`).
* Only exported functions of SGDrivE (no `:::`), so each script works as a
  user's script and as vignette code.
* A header saying what the scenario shows, then numbered steps with comments
  that explain the modelling choices.
* Packages beyond SGDrivE (MGDrivE for inheritance cubes, ggplot2 for plots)
  are loaded at the top.
* Stochastic runs set `seed`, so results are reproducible.
* Each script runs in about a minute or less on a laptop.
