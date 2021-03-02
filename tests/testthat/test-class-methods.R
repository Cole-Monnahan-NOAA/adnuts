## Note these were tested with ADMB 12.2 on 8/3/2020
test_that("class methods working", {
  fit <- readRDS('fit.RDS')
  x <- as.data.frame(fit)
  expect_is(x, 'data.frame')
  summary(fit)
  print(fit)
  plot_marginals(fit)
})
