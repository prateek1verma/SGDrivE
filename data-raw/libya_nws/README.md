# Libya NWS data (1989-1991)

Inputs for the Libyan validation (see `libya_release_schedule.csv`,
`libya_monthly_cases.csv`, `libya_anchor_data.csv`).

## Weather: `tripoli_airport_daily_1989_1991.csv`

Daily weather, 1989-07-01 to 1991-08-31 (792 days).

* Source: NOAA NCEI GHCN-Daily, station `LYE00147744` (Tripoli Airport,
  32.67 N, 13.15 E, 81 m), downloaded from
  `https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/LYE00147744.csv`.
* Columns: `date`, `tmax`, `tmin`, `tmean` (°C), `rain` (mm), `source`.
* `tmean` is the observed TAVG where available, otherwise (TMAX + TMIN)/2.
* Coverage in the window: TAVG is nearly complete, but TMAX and TMIN are
  observed on only about half of the days, so many days needed filling.
* `source` flags:
  * `primary`: TMAX and TMIN observed (249 days).
  * `dtr_fill`: TMAX and TMIN estimated from TAVG and the monthly mean
    diurnal range (TMAX − TMIN) from all complete days in the full record
    (1943-2014). Tmax = TAVG + DTR/2, Tmin = TAVG − DTR/2 (541 days).
  * `no_record`: no GHCN record for that date (2 days). Values are blank.
* `rain` is blank on 90 days (88 missing in GHCN, plus the 2 no-record days).
* Nearest alternatives were checked and are less suitable:
  `LYE00147743` (Tripoli City) and `LYE00147742` (Sidi El Mesri) have no
  data in this window. `LYM00062012` (El Khoms, ~60 km east) starts
  June 1990. `LYE00147745` (Misrata, ~190 km east) covers the window but
  is far from the region.
* Al-Aziziyah (~41-43 km southwest of Tripoli) has no station of its own in
  GHCN; Tripoli Airport is the nearest available record.

The DTR-based fill is a judgement call, not a verified observation. A
cleaner option, if needed, is to fill from the neighbouring stations with
more complete records or from ERA5 reanalysis.
