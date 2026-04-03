# Bayesian inference of a TMB model using the no-U-turn sampler.

Draw Bayesian posterior samples from a Template Model Builder (TMB)
model using an MCMC algorithm. This function generates posterior samples
from which inference can be made. Adaptation schemes are used so
specification tuning parameters are not necessary, and parallel
execution reduces overall run time.

## Usage

``` r
sample_tmb(
  obj,
  iter = 2000,
  init,
  chains = 3,
  seeds = NULL,
  warmup = floor(iter/2),
  lower = NULL,
  upper = NULL,
  thin = 1,
  parallel = FALSE,
  cores = NULL,
  path = NULL,
  algorithm = "NUTS",
  laplace = FALSE,
  control = NULL,
  ...
)
```

## Arguments

- obj:

  A TMB model object.

- iter:

  The number of samples to draw.

- init:

  A list of lists containing the initial parameter vectors, one for each
  chain or a function. It is strongly recommended to initialize multiple
  chains from dispersed points. A of NULL signifies to use the starting
  values present in the model (i.e., `obj$par`) for all chains.

- chains:

  The number of chains to run.

- seeds:

  A vector of seeds, one for each chain.

- warmup:

  The number of warmup iterations.

- lower:

  A vector of lower bounds for parameters. Allowed values are -Inf and
  numeric.

- upper:

  A vector of upper bounds for parameters. Allowed values are Inf and
  numeric.

- thin:

  The thinning rate to apply to samples. Typically not used with NUTS.

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

- path:

  Path to model executable. Defaults to working directory. Often best to
  have model files in a separate subdirectory, particularly for
  parallel.

- algorithm:

  The algorithm to use. NUTS is the default and recommended one, but
  "RWM" for the random walk Metropolis sampler and "HMC" for the static
  HMC sampler are available. These last two are deprecated but may be of
  use in some situations. These algorithms require different arguments;
  see their help files for more information.

- laplace:

  Whether to use the Laplace approximation if some parameters are
  declared as random. Default is to turn off this functionality and
  integrate across all parameters with MCMC.

- control:

  A list to control the sampler. See details for further use.

- ...:

  Further arguments to be passed to samplers

## Value

A list containing the samples, and properties of the sampler useful for
diagnosing behavior and efficiency.

## Details

This function implements algorithm 6 of Hoffman and Gelman (2014), and
loosely follows package `rstan`. The step size can be adapted or
specified manually. The metric (i.e., mass matrix) can be unit diagonal,
adapted diagonal (default and recommended), or a dense matrix specified
by the user. Further control of algorithms can be specified with the
`control` argument. Elements are:

- adapt_delta:

  The target acceptance rate.

- metric:

  The mass metric to use. Options are: "unit" for a unit diagonal
  matrix; "diag" to estimate a diagonal matrix during warmup; a matrix
  to be used directly (in untransformed space).

- adapt_engaged:

  Whether adaptation of step size and metric is turned on.

- max_treedepth:

  Maximum treedepth for the NUTS algorithm.

- stepsize:

  The stepsize for the NUTS algorithm. If `NULL` it will be adapted
  during warmup.

## Warning

This is deprecated and will cease to exist in future releases

## See also

[`extract_samples`](https://cole-monnahan-noaa.github.io/adnuts/reference/extract_samples.md)
to extract samples and
[`launch_shinytmb`](https://cole-monnahan-noaa.github.io/adnuts/reference/launch_shinytmb.md)
to explore the results graphically which is a wrapper for the
[`launch_shinystan`](https://mc-stan.org/shinystan/reference/launch_shinystan.html)
function.

## Author

Cole Monnahan

## Examples

``` r
## Build a fake TMB object with objective & gradient functions and some
## other flags
if (FALSE) { # \dontrun{
f <- function(x, order=0){
  if(order != 1) # negative log density
    -sum(dnorm(x=x, mean=0, sd=1, log=TRUE))
  else x # gradient of negative log density
}
init <- function() rnorm(2)
obj <- list(env=list(DLL='demo', last.par.best=c(x=init()), f=f,
  beSilent=function() NULL))
## Run NUTS for this object
fit <- sample_tmb(obj, iter=1000, chains=3, init=init)
## Check basic diagnostics
mon <- rstan::monitor(fit$samples, print=FALSE)
Rhat <- mon[,"Rhat"]
max(Rhat)
ess <- mon[, 'n_eff']
min(ess)
## Or do it interactively with ShinyStan
launch_shinytmb(fit)
} # }
```
