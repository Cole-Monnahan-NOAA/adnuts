# adnuts

main: [![R-CMD-check](https://github.com/Cole-Monnahan-NOAA/adnuts/workflows/R-CMD-check/badge.svg?branch=main)](https://github.com/Cole-Monnahan-NOAA/adnuts/actions?query=workflow%3AR-CMD-check)
dev: [![R-CMD-check](https://github.com/Cole-Monnahan-NOAA/adnuts/workflows/R-CMD-check/badge.svg?branch=dev)](https://github.com/Cole-Monnahan-NOAA/adnuts/actions?query=workflow%3AR-CMD-check) [![codecov](https://codecov.io/gh/Cole-Monnahan-NOAA/adnuts/branch/dev/graph/badge.svg)](https://codecov.io/gh/Cole-Monnahan-NOAA/adnuts)
[![](https://www.r-pkg.org/badges/version/adnuts)](https://www.r-pkg.org/pkg/adnuts)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/adnuts)](https://www.r-pkg.org/pkg/adnuts)

The aim of 'adnuts' (pronounced A-D NUTS like A-D MB) is to provide
advanced MCMC sampling for 'ADMB' and 'TMB' models. It mimics 'Stan' in
functionality and feel, specifically providing no-U-turn (NUTS) sampling
with adaptive mass matrix and parallel execution.

The R package 'tmbstan' (available on CRAN) has largely replaced the TMB
capabilities since original development. As such, adnuts is primarily used
for ADMB models. See the following paper for an introduction to the package
capabilities, and contrast with tmbstan:

Monnahan CC, Kristensen K (2018) No-U-turn sampling for fast Bayesian
inference in ADMB and TMB: Introducing the adnuts and tmbstan R
packages. PLoS ONE 13(5):e0197954. 
https://doi.org/10.1371/journal.pone.0197954

'adnuts' was designed specifically for use in fisheries stock assessments,
and interested authors are referred to:

Monnahan, C.C., T.A. Branch, J.T. Thorson, I.J. Stewart, C.S. Szuwalksi
(2020) Overcoming long Bayesian run times in integrated fisheries stock
assessments. ICES Journal of Marine
Science. https://dx.doi.org/10.1093/icesjms/fsz059


## Usage
The 'sample_rwm' and 'sample_nuts' functions draw posterior samples from an
ADMB model using an MCMC algorithm (random walk Metropolis or no-U-turn
sampler). The returned fitted object contains samples and other
information. The function 'extract_samples' can be used to get posterior
samples (post warmup and thinning) into a data frame for inference, while
'launch_shinyadmb' can be used for interactive diagnostics based on
'ShinyStan'.

A brief [demonstration
file](https://github.com/Cole-Monnahan-NOAA/adnuts/blob/master/inst/demo.R)
is the best place to help get you started, and there is also a user guide:
`vignette('adnuts')` for more detailed information.

## Installation

To use the ADMB functionality you need to build your model with version
12.0 (released December 2017) or later, otherwise this functionality is not
available. See [the ADMB installation
instructions](https://www.admb-project.org/docs/install/) for more
information. ADMB 12.2 is highly recommended because it provides better
console output, fixes bugs, and adds improved adaptation capabilities as
compared to 12.0. You can check the ADMB version of a compiled model from
the command line with a command `model.exe -version` which prints the
version among other things.

The adnuts R package version 1.1.2 can be installed from CRAN:
`install.packages('adnuts')`. Future minor releases [listed
here](https://github.com/Cole-Monnahan-NOAA/adnuts/releases) may not be
released on CRAN so the latest stable version can be installed as:

`devtools::install_github('Cole-Monnahan-NOAA/adnuts')`

The development version can be installed as:
`devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='dev')`

## Known issues
Windows users may experience issues if their model name is too long. In
some cases the OS will rename the output files using a "short"
version. You'll see files like "MODEL~1.par". The package tries to handle
this but it is **highly recommended** to simply shorten your filename. So
instead of 'model_filename_2021.tpl' use e.g. 'model_21'.

Analyses are reproducible by setting the same initial values and a seed in
`sample_rwm` or `sample_nuts` (passed to ADMB as '-mcseed'). However, they
may not be entirely consistent across OS platforms. The chains will start
the same but may eventually diverge. This is likely due to minuscule
differences in the gradient and log-posterior calculations between systems
and compilers.

## Disclaimer

“The United States Department of Commerce (DOC) GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its
use. DOC has relinquished control of the information and no longer has
responsibility to protect the integrity, confidentiality, or availability
of the information. Any claims against the Department of Commerce stemming
from the use of its GitHub project will be governed by all applicable
Federal law. Any reference to specific commercial products, processes, or
services by service mark, trademark, manufacturer, or otherwise, does not
constitute or imply their endorsement, recommendation or favoring by the
Department of Commerce. The Department of Commerce seal and logo, or the
seal and logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.”

