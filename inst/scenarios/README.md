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
