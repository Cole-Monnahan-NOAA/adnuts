test_that("all three metrics work", {
 skip_if(skip_TMB)
  TMB::runExample('simple')
  Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  M <- as.matrix(solve(Q))
  fits <- list()
  library(StanEstimators)
  for(m in c('dense', 'sparse', 'diag')){
    fits[[m]] <- sample_sparse_tmb(obj, iter=1000,
                                   warmup=200, cores=1, chains=1, seed=1,
                                   metric=m)
  }
 expect_equal(length(fits),3)
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
  fit <- sample_sparse_tmb(obj, iter=1000, warmup=200, cores=4, chains=4, seed=1)
})

test_that("RTMB works", {
  skip_if(TRUE) # not sure why this fails when testing but not locally?
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

