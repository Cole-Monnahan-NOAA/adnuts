# adnuts
The aim of 'adnuts' (pronounced A-D NUTS like A-D MB) is to provide
advanced MCMC sampling for 'ADMB' and 'TMB' models. It mimics 'Stan' in
functionality and feel, specifically providing no-U-turn (NUTS) sampling
with adaptive mass matrix and parallel execution.

The R package 'tmbstan' (available on CRAN) has largely replaced the TMB
capabilities since original development. As such, adnuts is primiarly used
for ADMB models. See the following paper for an introduction to the package
capabilities, and contrast with tmbstan:

Monnahan CC, Kristensen K (2018) No-U-turn sampling for fast Bayesian
inference in ADMB and TMB: Introducing the adnuts and tmbstan R
packages. PLoS ONE 13(5):e0197954. https://doi.org/10.1371/journal.pone.0197954

## Usage
The 'sample_admb' function draws posterior samples from an ADMB model using
an MCMC algorithm (NUTS by default). The returned fitted object contains
samples and other information. The function 'extract_samples' can be used
to get posterior samples (post warmup and thinning) into a data frame for
inference, while 'launch_shinyadmb' can be used for interactive diagnostics
based on 'ShinyStan'.

A brief [demonstration file](https://github.com/colemonnahan/adnuts/blob/master/inst/demo.R) is
provided to help get you started, and there is also a user guide:
`vignette('adnuts')` for more detailed information.

## Installation

To use the ADMB functionality you need to build your model with version
12.0 (released December 2017) or later, otherwise this functionality is not
available. See [the ADMB installation
instructions](http://www.admb-project.org/docs/install/) for more
information.

The adnuts R package version 1.0.1 can be installed from CRAN:
`install.packages('adnuts')`.

The development version of 'adnuts' can be installed with
`devtools::install_github('colemonnahan/adnuts', ref='dev')`
