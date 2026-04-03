# Package index

## All functions

- [`adfit()`](adfit.md) : Constructor for the "adfit" (A-D fit) class

- [`adnuts`](adnuts.md) : adnuts: No-U-turn sampling for AD Model
  Builder (ADMB)

- [`as.data.frame(`*`<adfit>`*`)`](as.data.frame.adfit.md) :

  Convert object of class adfit to data.frame. Calls `extract_samples`

- [`check_identifiable()`](check_identifiable.md) : Check
  identifiability from model Hessian

- [`.check_ADMB_version()`](dot-check_ADMB_version.md) : Check that the
  model is compiled with the right version of ADMB which is 12.0 or
  later

- [`.check_console_printing()`](dot-check_console_printing.md) : Check
  if the session is interactive or Rstudio which has implications for
  parallel output

- [`.check_model_path()`](dot-check_model_path.md) : Check that the file
  can be found

- [`.getADMBHessian()`](dot-getADMBHessian.md) : Read in admodel.hes
  file

- [`.sample_admb()`](dot-sample_admb.md) : Hidden wrapper function for
  sampling from ADMB models

- [`.update_model()`](dot-update_model.md) : Convert model name
  depending on system

- [`extract_sampler_params()`](extract_sampler_params.md) : Extract
  sampler parameters from a fit.

- [`extract_samples()`](extract_samples.md) : Extract posterior samples
  from a model fit.

- [`is.adfit()`](is.adfit.md) : Check object of class adfit

- [`launch_shinyadmb()`](launch_shinyadmb.md) : Launch shinystan for an
  ADMB fit.

- [`launch_shinytmb()`](launch_shinytmb.md) : Launch shinystan for a TMB
  fit.

- [`pairs_admb()`](pairs_admb.md) : Plot pairwise parameter posteriors
  and optionally the MLE points and confidence ellipses.

- [`plot(`*`<adfit>`*`)`](plot.adfit.md) : Plot object of class adfit

- [`plot_marginals()`](plot_marginals.md) : Plot marginal distributions
  for a fitted model

- [`plot_sampler_params()`](plot_sampler_params.md) : Plot adaptation
  metrics for a fitted model.

- [`plot_uncertainties()`](plot_uncertainties.md) : Plot MLE vs MCMC
  marginal standard deviations for each parameter

- [`print(`*`<adfit>`*`)`](print.adfit.md) : Print summary of adfit
  object

- [`sample_admb()`](sample_admb.md) : Deprecated version of wrapper
  function. Use sample_nuts or sample_rwm instead.

- [`sample_inits()`](sample_inits.md) : Function to generate random
  initial values from a previous fit using adnuts

- [`sample_tmb()`](sample_tmb.md) : Bayesian inference of a TMB model
  using the no-U-turn sampler.

- [`sample_tmb_hmc()`](sample_tmb_hmc.md) : Draw MCMC samples from a
  model posterior using a static HMC sampler.

- [`sample_tmb_nuts()`](sample_tmb_nuts.md) : Draw MCMC samples from a
  model posterior using the No-U-Turn (NUTS) sampler with dual
  averaging.

- [`sample_tmb_rwm()`](sample_tmb_rwm.md) : \[Deprecated\] Draw MCMC
  samples from a model posterior using a Random Walk Metropolis (RWM)
  sampler.

- [`summary(`*`<adfit>`*`)`](summary.adfit.md) : Print summary of object
  of class adfit

- [`sample_nuts()`](wrappers.md) [`sample_rwm()`](wrappers.md) :
  Bayesian inference of an ADMB model using the no-U-turn sampler (NUTS)
  or random walk Metropolis (RWM) algorithms.
