# data-raw

Scripts and inputs that produce data shipped with, or used to test, ScrewDrivE.
This folder is excluded from the built package (`.Rbuildignore`).

## Contents

| Path | What it is |
|---|---|
| `reference_impl/scripts/` | **Unchanged** copies of the legacy R translations of Gutierrez & Ponti (2014): `gutierrez_ponti_2014_pbdm.R` (line-by-line port of the Delphi Pascal `ScrewWormPop` unit; discrete daily distributed delay) and `gp2014_ode.R` (Erlang-chain ODE counterpart). Do not edit: they are the regression baseline. |
| `reference_impl/mcallen_1960_1983_filled.csv` | Weather input used for the reference runs (same file as `inst/extdata/mcallen_tx_1960_1983.csv`). |
| `make_reference.R` | Regenerates `reference/`: `Rscript data-raw/make_reference.R` (about 1 minute). |
| `reference/*.rds` | Reference outputs on McAllen weather: `discrete`, `ode`, `ode_no_floors`, `ode_match_larvae`, plus the parameter list. |
| `make_test_fixtures.R` | Builds `tests/testthat/fixtures/`: legacy temperature functions on a grid, legacy McAllen forcing terms, and short raw NOAA extracts around backup-filled and interpolated days (needs the raw NOAA files; see the script). |
| `make_reference_tight.R` | Tight-tolerance legacy run without floors (rtol 1e-10, atol 1e-12) for the M2 regression (D-019): writes `reference/gp2014_mcallen_ode_no_floors_tight.rds` (1960-1983), `reference/PROVENANCE_tight.txt` and the test fixture `tests/testthat/fixtures/gp2014_mcallen_ode_tight_3yr.rds` (daily totals 1960-1962). About 1 minute. A 3-year rerun at rtol 1e-11 changes the legacy solution by at most 2e-9 (recorded in the provenance). |
| `reference/PROVENANCE.txt` | Date, R and deSolve versions, md5 of inputs and scripts, run times. |

## Weather data (`inst/extdata/mcallen_tx_1960_1983.csv`)

Daily weather for McAllen, Texas, 1960-01-01 to 1983-12-31.

* Source: NOAA NCEI GHCN-Daily, station `USC00415701` (McAllen), with gaps
  filled from `USW00012959` (McAllen Miller International Airport, < 1 km
  away), then by the mean of the neighbouring days, then by the mean of the
  preceding and following 5 days. This follows Gutierrez et al. (2019),
  Appendix S1. The gap-filling function is `gp_prepare_ghcn()` in
  `reference_impl/scripts/gutierrez_ponti_2014_pbdm.R`.
* Raw files: <https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/USC00415701.csv>
  and `.../USW00012959.csv` (values in tenths of °C and tenths of mm).
* Columns: `date`, `tmax`, `tmin` (°C), `rain` (mm/day), `source`
  (`primary`, `backup` or `interpolated`), `tmean` = (tmax + tmin)/2.
* Days by source: 8612 primary, 152 backup, 2 interpolated.
* Check: the cold-stress index computed from this file matches the values
  reported by Gutierrez & Ponti / Gutierrez et al. (e.g. 1971 = 6.67,
  1979 = 9.80).
