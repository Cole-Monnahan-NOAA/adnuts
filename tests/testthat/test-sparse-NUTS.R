test_that("sparse metric works", {
 skip_if(skip_TMB)
    TMB::runExample('simple')
  Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  M <- as.matrix(solve(Q))
  fits <- list()
  library(StanEstimators)
  for(m in c('dense', 'sparse', 'diag')){
    fits[[m]] <- sample_sparse_tmb(obj, iter=1000,
                                   warmup=200, cores=1, chains=1, seed=1,
                                   metric=m, skip_optimization = TRUE,
                                   Q=Q, Qinv=M)
  }
 expect_equal(length(fits),3)
})
