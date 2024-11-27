test_that("all three metrics work", {
 skip_if(skip_TMB)
  TMB::runExample('simple')
  Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  M <- as.matrix(solve(Q))
  fits <- list()
  for(m in c('dense', 'sparse', 'diag', 'unit')){
    fits[[m]] <- sample_sparse_tmb(obj, iter=1000,
                                   warmup=200, cores=1, chains=1, seed=1,
                                   metric=m)
  }
  expect_equal(length(fits),4)
  out <- lapply(fits, function(x) as.numeric(tail(as.data.frame(x), n=1)[1]))
  expect_equal(out$dense,-1.802417, tolerance=1e-5)
  expect_equal(out$sparse,-1.802417, tolerance=1e-5)
  expect_equal(out$diag,-0.3480269, tolerance=1e-5)
  expect_equal(out$unit,-0.0162497, tolerance=1e-5)
})

test_that("Laplace works", {
  skip_if(skip_TMB)
  TMB::runExample('simple')
  fit1 <- sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense')
   # has RE and uses Laplace
  fit2 <- sample_sparse_tmb(obj, iter=1000, laplace=TRUE,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense')
  expect_equal(ncol(as.data.frame(fit1)),118)
  expect_equal(ncol(as.data.frame(fit2)),4)
})

test_that("metrics are robust to model type",{
  skip_if(skip_TMB)
  TMB::runExample('simple')
  ## rebuild without RE so it fails
  obj <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
  expect_error(sample_sparse_tmb(obj, iter=1000, laplace=TRUE,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense'),
               regexp = 'No random effects found')
  expect_error(sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='sparse'),
               regexp='sparse metric only allowed with random effects')
  ## should fail since M not available
  expect_error(sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense'),
               regexp = 'Some standard errors estimated to be NaN'
              )
  ## should work
  fit4 <- sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='unit')
  expect_equal(ncol(as.data.frame(fit4)),118)
  ## should fail since M not available
  expect_error(sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='diag'),
               "Some standard errors estimated to be NaN")
  ## should work if variance term is turned off (penalized ML)
  obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                         map=list(logsdu=factor(NA)),
                         random=NULL, silent=TRUE,
                         DLL=obj$env$DLL)
  fit6 <- sample_sparse_tmb(obj2, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense')
 expect_equal(ncol(as.data.frame(fit6)), 117)
})

test_that("parallel works", {
  skip_if(skip_TMB)
  TMB::runExample('simple')
  fit <- sample_sparse_tmb(obj, iter=1000, warmup=200, cores=4, chains=4, seed=1, metric='sparse')
})


test_that("auto metric selection is robust to model type",{
  skip_if(skip_TMB)
  TMB::runExample('simple')
  ## normal case of RE, with and without laplace
  fit1 <- sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='auto')
  expect_equal(as.numeric(tail(as.data.frame(fit1),1)[1]),-1.802417, tolerance =1e-6)
  fit2 <- sample_sparse_tmb(obj, iter=1000,  laplace=TRUE,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='auto')
  expect_equal(as.numeric(tail(as.data.frame(fit2),1)[1]), 51.80767, tolerance =1e-6)

  ## rebuild without RE so it fails: behavior on model w/o mode
  TMB::runExample('simple')
  obj <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
  expect_error(sample_sparse_tmb(obj, iter=1000, laplace=TRUE,
                                 warmup=200, cores=1, chains=1, seed=1,
                                 metric='auto'),
               regexp = 'No random effects found')
  # this breaks if init='last.par.best' b/c the inits are so bad it can't recover
  fit3 <- sample_sparse_tmb(obj, iter=1000, laplace=FALSE, init='random',
                                 warmup=200, cores=1, chains=1, seed=1,
                                 metric='auto', skip_optimization = TRUE)
  expect_equal(as.numeric(tail(as.data.frame(fit3),1)[1]), -1.24197, tolerance =1e-6)
  ## rebuild as penalized ML
  TMB::runExample('simple')
  obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                         map=list(logsdu=factor(NA)),
                         random=NULL, silent=TRUE,
                         DLL=obj$env$DLL)
  expect_error(sample_sparse_tmb(obj2, iter=1000, laplace=TRUE,
                                 warmup=200, cores=1, chains=1, seed=1,
                                 metric='auto'),
               regexp = 'No random effects found')
  # this breaks if init='last.par.best' b/c the inits are so bad it can't recover
  fit5 <- sample_sparse_tmb(obj2, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='auto')
  #expect_equal(as.numeric(tail(as.data.frame(fit5),1)[1]), -0.8898186, tolerance =1e-6)
})


test_that("RTMB works", {
  skip_if(skip_RTMB) # not sure why this fails when testing but not locally?
  ## from RTMB beginner help page
  ## suppressWarnings(detach("package:TMB", unload = TRUE))
  library(RTMB)
  data(ChickWeight)
  parameters <- list(
    mua=0,          ## Mean slope
    sda=1,          ## Std of slopes
    mub=0,          ## Mean intercept
    sdb=1,          ## Std of intercepts
    sdeps=1,        ## Residual Std
    a=rep(0, 50),   ## Random slope by chick
    b=rep(0, 50)    ## Random intercept by chick
  )
  f <- function(parms) {
    getAll(ChickWeight, parms, warn=FALSE)
    ## Optional (enables extra RTMB features)
    weight <- OBS(weight)
    ## Initialize joint negative log likelihood
    nll <- 0
    ## Random slopes
    nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
    ## Random intercepts
    nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
    ## Data
    predWeight <- a[Chick] * Time + b[Chick]
    nll <- nll - sum(dnorm(weight, predWeight, sd=sdeps, log=TRUE))
    ## Get predicted weight uncertainties
    ADREPORT(predWeight)
    ## Return
    nll
  }
  f(parameters)
  obj <- RTMB::MakeADFun(f, parameters, random=c("a", "b"))
  #obj$fn()
  #fit <- sample_sparse_tmb(obj, iter=2000, warmup=1000, chains=4, cores=4, seed=1)
  #expect_equal(ncol(as.data.frame(fit)),105)
})

