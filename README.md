# adnuts

main: [![R-CMD-check](https://github.com/Cole-Monnahan-NOAA/adnuts/workflows/R-CMD-check/badge.svg?branch=main)](https://github.com/Cole-Monnahan-NOAA/adnuts/actions?query=workflow%3AR-CMD-check) dev: [![R-CMD-check](https://github.com/Cole-Monnahan-NOAA/adnuts/workflows/R-CMD-check/badge.svg?branch=dev)](https://github.com/Cole-Monnahan-NOAA/adnuts/actions?query=workflow%3AR-CMD-check) [![codecov](https://codecov.io/gh/Cole-Monnahan-NOAA/adnuts/branch/dev/graph/badge.svg)](https://codecov.io/gh/Cole-Monnahan-NOAA/adnuts) [![](https://www.r-pkg.org/badges/version/adnuts)](https://www.r-pkg.org/pkg/adnuts) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/adnuts)](https://www.r-pkg.org/pkg/adnuts)

The aim of 'adnuts' (pronounced A-D nuts) is to provide advanced MCMC sampling for 'TMB' and 'ADMB' models. For TMB models it uses the sparse NUTS (SNUTS; Monnahan et al. in prep) algorithm to decorrelate the posterior using the joint precision matrix. The R package 'tmbstan' (available on CRAN) provides an alternative for TMB which more closely links to Stan. Until the development of SNUTS, 'adnuts' was primarily used for ADMB models. For the foreseeable future SNUTS via 'adnuts' is likely to be the best general option for TMB users.

For ADMB it mimics 'Stan' in functionality and feel, specifically providing no-U-turn (NUTS) sampling with adaptive mass matrix and parallel execution. Development of ADMB features is winding down, but functionality expected to be maintained in the coming years.

See the following papers for an introduction to the package capabilities, and contrast with tmbstan:

Monnahan CC, Kristensen K (2018) No-U-turn sampling for fast Bayesian inference in ADMB and TMB: Introducing the adnuts and tmbstan R packages. PLoS ONE 13(5):e0197954. <https://doi.org/10.1371/journal.pone.0197954>

Monnahan CC, Thorson, J.T., Kristensen, K, and Carpenter, B (in prep). Leveraging sparsity to improve no-u-turn sampling efficiency of hierarchical Bayesian models.

## Installation


To use SNUTS with TMB first install the [StanEstimators](https://github.com/andrjohns/StanEstimators) package which is not on CRAN but can be installed as:

```         
# we recommend running this is a fresh R session or restarting your current session
install.packages('StanEstimators', repos = c('<https://andrjohns.r-universe.dev>', '<https://cloud.r-project.org>'))
```

Then install the GitHub version of 'adnuts':

````
devtools::install_github('Cole-Monnahan-NOAA/adnuts')

````

A quick test of functionality is:

```         
library(adnuts)
TMB::runExample('simple')
mcmc <- sample_snuts(obj, num_samples=500)
```

A brief [demonstration file](https://github.com/Cole-Monnahan-NOAA/adnuts/blob/dev/inst/demo_SNUTS.R) is the best place to help get you started, and there is also a user guide: see `vignette('adnuts')` for the basics and this [online article](https://cole-monnahan-noaa.github.io/adnuts/docs/articles/SNUTS-for-TMB-models.html) for more detailed information.



## Disclaimer

“The United States Department of Commerce (DOC) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. DOC has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce stemming from the use of its GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.”

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" alt="NOAA Fisheries" height="75"/>

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National Oceanographic and Atmospheric Administration](https://www.noaa.gov) \| [NOAA Fisheries](https://www.fisheries.noaa.gov/)
