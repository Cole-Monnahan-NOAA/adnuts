# Draw MCMC samples from a model posterior using the No-U-Turn (NUTS) sampler with dual averaging.

Draw MCMC samples from a model posterior using the No-U-Turn (NUTS)
sampler with dual averaging.

## Usage

``` r
sample_tmb_nuts(
  iter,
  fn,
  gr,
  init,
  warmup = floor(iter/2),
  chain = 1,
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

- gr:

  A function that returns a vector of gradients of the log of the
  posterior density (same as `fn`).

- init:

  A list of lists containing the initial parameter vectors, one for each
  chain or a function. It is strongly recommended to initialize multiple
  chains from dispersed points. A of NULL signifies to use the starting
  values present in the model (i.e., `obj$par`) for all chains.

- warmup:

  The number of warmup iterations.

- chain:

  The chain number, for printing only.

- thin:

  The thinning rate to apply to samples. Typically not used with NUTS.

- seed:

  The random seed to use.

- control:

  A list to control the sampler. See details for further use.

## Details

This function implements algorithm 6 of Hoffman and Gelman (2014), which
includes adaptive step sizes (`eps`) via an algorithm called dual
averaging. It also includes an adaptation scheme to tune a diagonal mass
matrix (metric) during warmup.

These `fn` and `gr` functions must have Jacobians already applied if
there are transformations used.

## References

Hoffman and Gelman (2014). The No-U-Turn sampler: Adaptively setting
path lengths in Hamiltonian Monte Carlo. J. Mach. Learn. Res.
15:1593-1623.

## See also

`sample_tmb`
