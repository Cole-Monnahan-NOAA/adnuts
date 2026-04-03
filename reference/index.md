# Package index

## All functions

- [`adfit()`](https://cole-monnahan-noaa.github.io/adnuts/reference/adfit.md)
  : Constructor for the "adfit" (A-D fit) class

- [`adnuts`](https://cole-monnahan-noaa.github.io/adnuts/reference/adnuts.md)
  : adnuts: No-U-turn sampling for AD Model Builder (ADMB)

- [`as.data.frame(`*`<adfit>`*`)`](https://cole-monnahan-noaa.github.io/adnuts/reference/as.data.frame.adfit.md)
  :

  Convert object of class adfit to data.frame. Calls `extract_samples`

- [`check_identifiable()`](https://cole-monnahan-noaa.github.io/adnuts/reference/check_identifiable.md)
  : Check identifiability from model Hessian

- [`.check_ADMB_version()`](https://cole-monnahan-noaa.github.io/adnuts/reference/dot-check_ADMB_version.md)
  : Check that the model is compiled with the right version of ADMB
  which is 12.0 or later

- [`.check_console_printing()`](https://cole-monnahan-noaa.github.io/adnuts/reference/dot-check_console_printing.md)
  : Check if the session is interactive or Rstudio which has
  implications for parallel output

- [`.check_model_path()`](https://cole-monnahan-noaa.github.io/adnuts/reference/dot-check_model_path.md)
  : Check that the file can be found

- [`.getADMBHessian()`](https://cole-monnahan-noaa.github.io/adnuts/reference/dot-getADMBHessian.md)
  : Read in admodel.hes file

- [`.sample_admb()`](https://cole-monnahan-noaa.github.io/adnuts/reference/dot-sample_admb.md)
  : Hidden wrapper function for sampling from ADMB models

- [`.update_model()`](https://cole-monnahan-noaa.github.io/adnuts/reference/dot-update_model.md)
  : Convert model name depending on system

- [`extract_sampler_params()`](https://cole-monnahan-noaa.github.io/adnuts/reference/extract_sampler_params.md)
  : Extract sampler parameters from a fit.

- [`extract_samples()`](https://cole-monnahan-noaa.github.io/adnuts/reference/extract_samples.md)
  : Extract posterior samples from a model fit.

- [`is.adfit()`](https://cole-monnahan-noaa.github.io/adnuts/reference/is.adfit.md)
  : Check object of class adfit

- [`launch_shinyadmb()`](https://cole-monnahan-noaa.github.io/adnuts/reference/launch_shinyadmb.md)
  : Launch shinystan for an ADMB fit.

- [`launch_shinytmb()`](https://cole-monnahan-noaa.github.io/adnuts/reference/launch_shinytmb.md)
  : Launch shinystan for a TMB fit.

- [`pairs_admb()`](https://cole-monnahan-noaa.github.io/adnuts/reference/pairs_admb.md)
  : Plot pairwise parameter posteriors and optionally the MLE points and
  confidence ellipses.

- [`plot(`*`<adfit>`*`)`](https://cole-monnahan-noaa.github.io/adnuts/reference/plot.adfit.md)
  : Plot object of class adfit

- [`plot_marginals()`](https://cole-monnahan-noaa.github.io/adnuts/reference/plot_marginals.md)
  : Plot marginal distributions for a fitted model

- [`plot_sampler_params()`](https://cole-monnahan-noaa.github.io/adnuts/reference/plot_sampler_params.md)
  : Plot adaptation metrics for a fitted model.

- [`plot_uncertainties()`](https://cole-monnahan-noaa.github.io/adnuts/reference/plot_uncertainties.md)
  : Plot MLE vs MCMC marginal standard deviations for each parameter

- [`print(`*`<adfit>`*`)`](https://cole-monnahan-noaa.github.io/adnuts/reference/print.adfit.md)
  : Print summary of adfit object

- [`sample_admb()`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_admb.md)
  : Deprecated version of wrapper function. Use sample_nuts or
  sample_rwm instead.

- [`sample_inits()`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_inits.md)
  : Function to generate random initial values from a previous fit using
  adnuts

- [`sample_tmb()`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_tmb.md)
  : Bayesian inference of a TMB model using the no-U-turn sampler.

- [`sample_tmb_hmc()`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_tmb_hmc.md)
  : Draw MCMC samples from a model posterior using a static HMC sampler.

- [`sample_tmb_nuts()`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_tmb_nuts.md)
  : Draw MCMC samples from a model posterior using the No-U-Turn (NUTS)
  sampler with dual averaging.

- [`sample_tmb_rwm()`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_tmb_rwm.md)
  : \[Deprecated\] Draw MCMC samples from a model posterior using a
  Random Walk Metropolis (RWM) sampler.

- [`summary(`*`<adfit>`*`)`](https://cole-monnahan-noaa.github.io/adnuts/reference/summary.adfit.md)
  : Print summary of object of class adfit

- [`sample_nuts()`](https://cole-monnahan-noaa.github.io/adnuts/reference/wrappers.md)
  [`sample_rwm()`](https://cole-monnahan-noaa.github.io/adnuts/reference/wrappers.md)
  : Bayesian inference of an ADMB model using the no-U-turn sampler
  (NUTS) or random walk Metropolis (RWM) algorithms.
