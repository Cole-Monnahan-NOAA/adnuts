# adnuts
The aim of 'adnuts' is to provide advanced MCMC sampling for 'ADMB' and 'TMB' models. It  mimics 'Stan' in functionality and feel, specifically providing no-U-turn (NUTS) sampling with adaptive mass matrix and parallel execution.

## Usage
The 'sample_admb' and 'sample_tmb' functions draw posterior samples using an MCMC algorithm (NUTS by default). The returned fitted object contains samples and other information. The function 'extract_samples' can be used to get posterior samples (post warmup and thinning) into a data frame for inference, while 'launch_shinyadmb' can be used for interactive diagnostics based on 'ShinyStan'. 

A brief [demonstration file](https://github.com/colemonnahan/adnuts/blob/master/inst/demo.R) is provided to help get you started, and there is also a user guide: `vignette('adnuts')` for more detailed information

## Installation

The adnuts R package can be installed from CRAN: `install.packages('adnuts')`. To use the ADMB functionality you need to build your model with version 12.0 (released December 2017) or later, otherwise this functionality is not available.
