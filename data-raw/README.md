# data-raw

Scripts and inputs that produce data shipped with, or used to test, SGDrivE.
This folder is excluded from the built package (`.Rbuildignore`).

## Contents

| Path | What it is |
|---|---|
| `reference_impl/scripts/` | **Unchanged** copies of the legacy R translations of Gutierrez & Ponti (2014): `gutierrez_ponti_2014_pbdm.R` (line-by-line port of the Delphi Pascal `ScrewWormPop` unit; discrete daily distributed delay) and `gp2014_ode.R` (Erlang-chain ODE counterpart). Do not edit: they are the regression baseline. |
| `reference_impl/mcallen_1960_1983_filled.csv` | Weather input used for the reference runs (same file as `inst/extdata/mcallen_tx_1960_1983.csv`). |
| `make_reference.R` | Regenerates `reference/`: `Rscript data-raw/make_reference.R` (about 1 minute). |
| `reference/*.rds` | Reference outputs on McAllen weather: `discrete`, `ode`, `ode_no_floors`, `ode_match_larvae`, plus the parameter list. |
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
