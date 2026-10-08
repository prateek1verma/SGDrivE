
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# SGDrivE <small>Screwworm Genetic-Control Dynamics Explorer</small>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/prateek1verma/SGDrivE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/prateek1verma/SGDrivE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

SGDrivE will simulate genotype-resolved population dynamics of the New
World screwworm (*Cochliomyia hominivorax*) under genetic control: the
sterile insect technique (SIT), precision-guided SIT (pgSIT),
female-specific RIDL (fsRIDL) and gene drives.

> **Status: under development (milestone M1, inputs, complete).** The
> package builds model inputs (parameters, inheritance cubes, weather
> forcing, releases) and assembles models, but has no solvers yet
> (milestones M2-M4). Do not use it for analysis.

## Planned model (v0.1)

| Component | Basis |
|----|----|
| Life history: degree-day development, temperature-dependent mortality, adult ageing chains, pre-oviposition period, wound-limited oviposition | Gutierrez & Ponti (2014) |
| Genotype inheritance, read from unmodified inheritance cubes | MGDrivE (Sanchez et al. 2020) |
| Mating: unmated → mated females, male availability | MGDrivE2 (Wu et al. 2021), adapted to NWS |
| Solvers: deterministic ODE and stochastic chain-binomial, one shared C++ rate kernel | SGDrivE |

Populations are tracked as counts per node. The full specification is in
[`dev/MODEL_SPEC.md`](dev/MODEL_SPEC.md) and the design decisions in
[`dev/DECISIONS.md`](dev/DECISIONS.md).

## Installation

The repository is private during development. With access:

``` r
# install.packages("remotes")
remotes::install_github("prateek1verma/SGDrivE")
```

A C++17 compiler is required (Xcode command line tools on macOS, Rtools
on Windows).

## Check the installation

``` r
library(SGDrivE)
sg_backend_info()
#> $package_version
#> [1] "0.0.0.9000"
#> 
#> $cplusplus
#> [1] 201703
#> 
#> $compiler
#> [1] "clang 21.0.0 (clang-2100.3.34.2)"
```

## Building a model (milestone M1)

``` r
# parameters (definitions: sg_param_table())
p <- sg_params()
# any MGDrivE inheritance cube
cube <- sg_cube(MGDrivE::cubeMendelian())
# daily forcing from the bundled McAllen, Texas weather
wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1972-12-31")
# weekly releases of 500 adult males of genotype aa
rel <- sg_releases("aa", "adult", "male", times = seq(60, by = 7, length.out = 20), number = 500)
mod <- sg_model(p, cube, wx, rel, area = 1)
mod
#> <sg_model> 207 states; 3 genotypes; 731 days (1971-01-01 to 1972-12-31)
#>   sub-stages: eggs 3, larvae 6, pupae 10, adults 10 per sex
#>   female chain: [0, 61) dd: 3; [61, 202.8) dd: 7
#>   male chain: [0, 202.8) dd: 10
#>   node area 1 km^2; wounds H = 5 per node; mating "rate"; search competitors "searching_females"
#>   releases: 20 events on 20 days (days 60 to 193), 10,000 individuals
#>   hash: a5844ca5f27e0c07362d7651d7b9b72a
```

The McAllen cold-stress index reproduces the values published by
Gutierrez & Ponti (2014) and Gutierrez et al. (2019):

``` r
ci <- sg_cold_index(sg_example_weather())
ci[ci$year %in% c(1971, 1979), c("year", "mu_cold")]
#>    year  mu_cold
#> 12 1971 6.667763
#> 20 1979 9.798071
```

## Roadmap

| Milestone | Content |
|----|----|
| M0 ✔ | Foundations: package skeleton, tests, linting, CI, design documents |
| M1 ✔ | Inputs: parameters, inheritance cubes, weather forcing, releases, model assembly |
| M2 | Reference genotype-resolved ODE in pure R |
| M3 | C++ rate kernel and ODE solver |
| M4 | C++ chain-binomial stochastic solver |
| M5 | Checks, summaries, calibration and plots |
| M6 | Scientific validation |
| M7 | Documentation and v0.1.0 release |

See [`dev/ROADMAP.md`](dev/ROADMAP.md) for completion criteria.

## References

Gutierrez, A.P. & Ponti, L. (2014). The new world screwworm: prospective
distribution and role of weather in eradication. *Agricultural and
Forest Entomology* 16:158–173. <https://doi.org/10.1111/afe.12046>

Sanchez C., H.M., Wu, S.L., Bennett, J.B. & Marshall, J.M. (2020).
MGDrivE: A modular simulation framework for the spread of gene drives
through spatially explicit mosquito populations. *Methods in Ecology and
Evolution* 11:229–239. <https://doi.org/10.1111/2041-210X.13318>

Wu, S.L. et al. (2021). MGDrivE 2: A simulation framework for gene drive
systems incorporating seasonality and epidemiological dynamics. *PLoS
Computational Biology* 17(5):e1009030.
<https://doi.org/10.1371/journal.pcbi.1009030>

## Licence

GPL (\>= 3).
