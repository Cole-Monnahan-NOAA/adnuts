# Bayesian inference of an ADMB model using the no-U-turn sampler (NUTS) or random walk Metropolis (RWM) algorithms.

Draw Bayesian posterior samples from an AD Model Builder (ADMB) model
using an MCMC algorithm. \`sample_nuts\` and \`sample_rwm\` generates
posterior samples from which inference can be made.

## Usage

``` r
sample_nuts(
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
  parallel = FALSE,
  cores = NULL,
  control = NULL,
  skip_optimization = TRUE,
  verbose = TRUE,
  skip_monitor = FALSE,
  skip_unbounded = TRUE,
  admb_args = NULL,
  extra.args = NULL
)

sample_rwm(
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
  parallel = FALSE,
  cores = NULL,
  control = NULL,
  skip_optimization = TRUE,
  verbose = TRUE,
  skip_monitor = FALSE,
  skip_unbounded = TRUE,
  admb_args = NULL,
  extra.args = NULL
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

- parallel:

  A deprecated argument, use cores=1 for serial execution or cores\>1
  for parallel (default is to parallel with cores equal to the
  available-1)

- cores:

  The number of cores to use for parallel execution. Default is number
  available in the system minus 1. If `cores=1`, serial execution occurs
  (even if `chains>1`), otherwise parallel execution via package
  snowfall is used. For slow analyses it is recommended to set
  `chains`\<=`cores` so each core needs to run only a single chain.

- control:

  A list to control the sampler. See details for further use.

- skip_optimization:

  Whether to run the optimizer before running MCMC. This is rarely need
  as it is better to run it once before to get the covariance matrix, or
  the estimates are not needed with adaptive NUTS.

- verbose:

  Flag whether to show console output (default) or suppress it
  completely except for warnings and errors. Works for serial or
  parallel execution.

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

- extra.args:

  Deprecated, use a `admb_args` instead.

## Details

Adaptation schemes are used with NUTS so specifying tuning parameters is
not necessary. See vignette for options for adaptation of step size and
mass matrix. The RWM algorithm provides no new functionality not
available from previous versions of ADMB. However, \`sample_rwm\` has an
improved console output, is setup for parallel execution, and a smooth
workflow for diagnostics.

Parallel chains will be run if argument \`cores\` is greater than one.
This entails copying the folder, and starting a new R session to run
that chain, which are then merged back together. Note that console
output is inconsistent when using parallel, and may not show. On Windows
the R terminal shows output live, but the GUI does not. RStudio is a
special case and will not show live, and instead is captured and
returned at the end. It is strongly recommended to start with serial
execution as debugging parallel chains is very difficult.

Note that the algorithm code is in the ADMB source code, and 'adnuts'
provides a wrapper for it. The command line arguments are returned and
can be examined by the user. See vignette for more information.

This function implements algorithm 6 of Hoffman and Gelman (2014), and
loosely follows package `rstan`. The step size can be adapted or
specified manually. The metric (i.e., mass matrix) can be unit diagonal,
adapted diagonal (default and recommended), a dense matrix specified by
the user, or an adapted dense matrix. Further control of algorithms can
be specified with the `control` argument. Elements are:

- adapt_delta:

  The target acceptance rate. D

- metric:

  The mass metric to use. Options are: "unit" for a unit diagonal
  matrix; `NULL` to estimate a diagonal matrix during warmup; a matrix
  to be used directly (in untransformed space).

- adapt_delta:

  Whether adaptation of step size is turned on.

- adapt_mass:

  Whether adaptation of mass matrix is turned on. Currently only allowed
  for diagonal metric.

- adapt_mass_dense:

  Whether dense adaptation of mass matrix is turned on.

- max_treedepth:

  Maximum treedepth for the NUTS algorithm.

- stepsize:

  The stepsize for the NUTS algorithm. If `NULL` it will be adapted
  during warmup.

- adapt_init_buffer:

  The initial buffer size during mass matrix adaptation where sample
  information is not used (default 50)

- adapt_term_buffer:

  The terminal buffer size (default 75) during mass matrix adaptation
  (final fast phase)

- adapt_window:

  The initial size of the mass matrix adaptation window, which gets
  doubled each time thereafter.

- refresh:

  The rate at which to refresh progress to the console. Defaults to even
  10 progress updates.

The adaptation scheme (step size and mass matrix) is based heavily on
those by the software Stan, and more details can be found in that
documentation and this vignette.

## Warning

The user is responsible for specifying the model properly (priors,
starting values, desired parameters fixed, etc.), as well as assessing
the convergence and validity of the resulting samples (e.g., through the
`coda` package), or with function
[`launch_shinytmb`](https://cole-monnahan-noaa.github.io/adnuts/reference/launch_shinytmb.md)
before making inference. Specifically, priors must be specified in the
template file for each parameter. Unspecified priors will be implicitly
uniform.

## Author

Cole Monnahan

## Examples

``` r
if (FALSE) { # \dontrun{
## This is the packaged simple regression model
path.simple <- system.file('examples', 'simple', package='adnuts')
## It is best to have your ADMB files in a separate folder and provide that
## path, so make a copy of the model folder locally.
path <- 'simple'
dir.create(path)
trash <- file.copy(from=list.files(path.simple, full.names=TRUE), to=path)
## Compile and run model
oldwd <- getwd()
setwd(path)
system('admb simple.tpl')
system('simple')
setwd('..')
init <- function() rnorm(2)
## Run NUTS with defaults
fit <- sample_nuts(model='simple', init=init, path=path)
unlink(path, TRUE) # cleanup folder
setwd(oldwd)
} # }
```
