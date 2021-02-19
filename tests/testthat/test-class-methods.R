## Note these were tested with ADMB 12.2 on 8/3/2020
test_that("class methods working", {
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, iter=500,
                     control=list(refresh=-1))
  x <- as.data.frame(fit)
  summary(fit)
  print(fit)
  plot_marginals(fit)
})
