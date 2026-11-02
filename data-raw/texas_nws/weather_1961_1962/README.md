# Texas weather, 1961-1962 (step 1 for the 1962 validation)

Daily weather for the Texas counties of the 1962 myiasis data, 1961-01-01 to 1962-12-31 (730 days).

## Files

* `county_station_map_1961_1962.csv`: one row per county (253). Gives the
  county centroid (from `texasScrewwormdat1962.xlsx`), the assigned GHCN-Daily
  station, the distance in km, and that station's count of valid TMAX, TMIN and
  PRCP days in 1961-1962.
* `station_daily_1961_1962.csv`: observed daily TMAX, TMIN (°C) and PRCP (mm)
  for the 170 assigned stations (124,100 rows). Blank values are missing in GHCN.
  Not gap-filled.
* `station_daily_1961_1962_filled.csv`: the same rows, gap-filled. Columns:
  `tmax`, `tmin`, `tmean` = (tmax + tmin)/2, `prcp`, and `tmax_src`, `tmin_src`,
  `prcp_src` flagging each value's origin. **Use this file for modelling.**

## Method

* Source: NOAA NCEI GHCN-Daily, via the NCEI access service (daily-summaries, metric units).
* Station pool: all Texas (`TX`) GHCN stations with TMAX, TMIN and PRCP records
  spanning 1961-1962 in the inventory (357 stations). Stations with fewer than 650
  of the 730 days with both TMAX and TMIN were dropped, leaving 170.
* Each county is assigned its nearest station in that pool (great-circle distance).
  Median distance 7.5 km, maximum 118.7 km (West Texas counties).
* Gap-filling, following the McAllen procedure in
  `data-raw/reference_impl/scripts/gutierrez_ponti_2014_pbdm.R`, in order:
  1. **backup**: the nearest other pool station with a value on the same day,
     within 100 km (used as-is, no adjustment);
  2. **interpolated**: the mean of the days immediately before and after;
  3. **interpolated_window**: the mean of the 5 days before and after.

## Gap-fill counts

| Variable | primary | backup | interpolated | interpolated_window | missing |
|---|---|---|---|---|---|
| TMAX | 122,880 | 1,205 | 4 | 11 | 0 |
| TMIN | 122,713 | 1,373 | 7 | 7 | 0 |
| PRCP | 123,719 | 381 | 0 | 0 | 0 |

## Gap-fill error (estimated)

The backup values are not adjusted for distance or elevation. To estimate the
error, I compared pairs of stations (40 sampled station-nearest-neighbour pairs,
median separation 44 km) on days when both were observed:

| Variable | Mean difference | SD of difference |
|---|---|---|
| TMAX | −0.14 °C | 3.1 °C |
| TMIN | 0.00 °C | 2.6 °C |
| PRCP | 0.13 mm | 8.4 mm |

These are rough, and the backup gaps are only about 1% of days, so the effect on
a monthly total should be small. But the substitution error is much larger than the
McAllen case, where the backup station was under 1 km away. Check the results
for counties with a backup distance over 50 km.

## Known cautions

* Several counties have a station more than 50 km away.
* Stations with fewer than 650 valid days were excluded even if they cover the
  county, so some counties are assigned a farther station than necessary.
* Precipitation is locally variable; the backup PRCP values are the least reliable.
