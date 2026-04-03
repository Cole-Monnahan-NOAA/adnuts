# Hidden wrapper function for sampling from ADMB models

Hidden wrapper function for sampling from ADMB models

## Usage

``` r
.sample_admb(
  model,
  path = getwd(),
  iter = 2000,
  init = NULL,
  chains = 3,
  warmup = NULL,
  seeds = NULL,
  thin = 1,
  mceval = FALSE,
  duration = NULL,
  cores = NULL,
  control = NULL,
  verbose = TRUE,
  algorithm = "NUTS",
  skip_optimization = TRUE,
  skip_monitor = FALSE,
  skip_unbounded = TRUE,
  admb_args = NULL
)
```

## Arguments

- model:

  Name of model (i.e., 'model' for model.tpl). For non-Windows systems
  this will automatically be converted to './model' internally. For
  Windows, long file names are sometimes shortened from e.g.,
  'long_model_filename' to 'LONG\_~1'. This should work, but will throw
  warnings. Please shorten the model name. See
  https://en.wikipedia.org/wiki/8.3_filename.

- path:

  Path to model executable. Defaults to working directory. Often best to
  have model files in a separate subdirectory, particularly for
  parallel.

- iter:

  The number of samples to draw.

- init:

  A list of lists containing the initial parameter vectors, one for each
  chain or a function. It is strongly recommended to initialize multiple
  chains from dispersed points. A of NULL signifies to use the starting
  values present in the model (i.e., `obj$par`) for all chains.

- chains:

  The number of chains to run.

- warmup:

  The number of warmup iterations.

- seeds:

  A vector of seeds, one for each chain.

- thin:

  The thinning rate to apply to samples. Typically not used with NUTS.

- mceval:

  Whether to run the model with `-mceval` on samples from merged chains.

- duration:

  The number of minutes after which the model will quit running.

- cores:

  The number of cores to use for parallel execution. Default is number
  available in the system minus 1. If `cores=1`, serial execution occurs
  (even if `chains>1`), otherwise parallel execution via package
  snowfall is used. For slow analyses it is recommended to set
  `chains`\<=`cores` so each core needs to run only a single chain.

- control:

  A list to control the sampler. See details for further use.

- verbose:

  Flag whether to show console output (default) or suppress it
  completely except for warnings and errors. Works for serial or
  parallel execution.

- algorithm:

  The algorithm to use, one of "NUTS" or "RWM"

- skip_optimization:

  Whether to run the optimizer before running MCMC. This is rarely need
  as it is better to run it once before to get the covariance matrix, or
  the estimates are not needed with adaptive NUTS.

- skip_monitor:

  Whether to skip calculating diagnostics (effective sample size, Rhat)
  via the
  [`rstan::monitor`](https://mc-stan.org/rstan/reference/monitor.html)
  function. This can be slow for models with high dimension or many
  iterations. The result is used in plots and summaries so it is
  recommended to turn on. If model run with `skip_monitor=FALSE` you can
  recreate it post-hoc by setting
  `fit$monitor=rstan::monitor(fit$samples, fit$warmup, print=FALSE)`.

- skip_unbounded:

  Whether to skip returning the unbounded version of the posterior
  samples in addition to the bounded ones. It may be advisable to set to
  FALSE for very large models to save space.

- admb_args:

  A character string which gets passed to the command line, allowing
  finer control
