# \[Deprecated\] Draw MCMC samples from a model posterior using a Random Walk Metropolis (RWM) sampler.

\[Deprecated\] Draw MCMC samples from a model posterior using a Random
Walk Metropolis (RWM) sampler.

## Usage

``` r
sample_tmb_rwm(
  iter,
  fn,
  init,
  alpha = 1,
  chain = 1,
  warmup = floor(iter/2),
  thin = 1,
  seed = NULL,
  control = NULL
)
```

## Arguments

- iter:

  The number of samples to draw.

- fn:

  A function that returns the log of the posterior density.

- init:

  A list of lists containing the initial parameter vectors, one for each
  chain or a function. It is strongly recommended to initialize multiple
  chains from dispersed points. A of NULL signifies to use the starting
  values present in the model (i.e., `obj$par`) for all chains.

- alpha:

  The amount to scale the proposal, i.e, Xnew=Xcur+alpha\*Xproposed
  where Xproposed is generated from a mean-zero multivariate normal.
  Varying `alpha` varies the acceptance rate.

- chain:

  The chain number, for printing only.

- warmup:

  The number of warmup iterations.

- thin:

  The thinning rate to apply to samples. Typically not used with NUTS.

- seed:

  The random seed to use.

- control:

  A list to control the sampler. See details for further use.

## Value

A list containing samples and other metadata.

## Details

This algorithm does not yet contain adaptation of `alpha` so some trial
and error may be required for efficient sampling.

## References

Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., Teller, A.H.,
Teller, E., 1953. Equation of state calculations by fast computing
machines. J Chem Phys. 21:1087-1092.

## See also

[`sample_tmb`](sample_tmb.md)
