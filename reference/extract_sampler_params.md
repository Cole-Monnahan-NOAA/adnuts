# Extract sampler parameters from a fit.

Extract information about NUTS trajectories, such as acceptance ratio
and treedepth, from a fitted object.

## Usage

``` r
extract_sampler_params(fit, inc_warmup = FALSE)
```

## Arguments

- fit:

  A list returned by `sample_admb`.

- inc_warmup:

  Whether to extract the warmup samples or not (default). Warmup samples
  should never be used for inference, but may be useful for diagnostics.

## Value

An invisible data.frame containing samples (rows) of each parameter
(columns). If multiple chains exist they will be rbinded together.

## Details

Each trajectory (iteration) in NUTS has associated information about the
trajectory: stepsize, acceptance ratio, treedepth, and number of
leapfrog steps. This function extracts these into a data.frame, which
may be useful for diagnosing issues in certain cases. In general, the
user should not need to examine them, or preferably should via
[`plot_sampler_params`](plot_sampler_params.md) or
[`launch_shinyadmb`](launch_shinyadmb.md).

## See also

[`launch_shinyadmb`](launch_shinyadmb.md).

## Examples

``` r
fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
sp <- extract_sampler_params(fit, inc_warmup=TRUE)
str(sp)
#> 'data.frame':    3000 obs. of  8 variables:
#>  $ chain        : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ iteration    : num  1 2 3 4 5 6 7 8 9 10 ...
#>  $ accept_stat__: num  1.33e-06 9.97e-01 9.62e-01 9.79e-01 7.32e-01 ...
#>  $ stepsize__   : num  0.0292 0.0288 0.0322 0.0404 0.0347 ...
#>  $ treedepth__  : num  4 4 2 4 1 3 1 3 4 2 ...
#>  $ n_leapfrog__ : num  15 15 3 11 1 7 1 7 11 3 ...
#>  $ divergent__  : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ energy__     : num  28.3 17.5 12 12.6 12.6 ...
```
