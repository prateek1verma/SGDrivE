
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

> **Status: under development (milestone M0, foundations).** The package
> builds and checks, but contains no model code yet. Do not use it for
> analysis.

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

## Roadmap

| Milestone | Content |
|----|----|
| **M0** | Foundations: package skeleton, tests, linting, CI, design documents |
| M1 | Inputs: parameters, inheritance cubes, weather forcing, releases, model assembly |
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
