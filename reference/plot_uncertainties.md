# Plot MLE vs MCMC marginal standard deviations for each parameter

Plot MLE vs MCMC marginal standard deviations for each parameter

## Usage

``` r
plot_uncertainties(fit, log = TRUE, plot = TRUE)
```

## Arguments

- fit:

  A fitted object returned by
  [`sample_admb`](https://cole-monnahan-noaa.github.io/adnuts/reference/sample_admb.md)

- log:

  Whether to plot the logarithm or not.

- plot:

  Whether to plot it or not.

## Value

Invisibly returns data.frame with parameter name and estimated
uncertainties.

## Details

It can be helpful to compare uncertainty estimates between the two
paradigms. This plots the marginal posterior standard deviation vs the
frequentist standard error estimated from the .cor file. Large
differences often indicate issues with one estimation method.

## Examples

``` r
fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
x <- plot_uncertainties(fit, plot=FALSE)

head(x)
#>   par   sd.post  sd.mle
#> a   a 0.2028936 0.15547
#> b   b 0.9112519 0.70394
```
