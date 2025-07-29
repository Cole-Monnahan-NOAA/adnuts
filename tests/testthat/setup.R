## Skip consistency and reproducibility tests? Only need to run
## these locally when ADMB changes.
skip_consistency <- skip_reproducibility <- skip_ADMB <- TRUE

# Only TMB is tested as of July 2025
skip_TMB <- FALSE
skip_RTMB <- FALSE

suppressWarnings(library(TMB))
# setup simple object once, have to be careful b/c of the dynamic
# links and I intentionally break the model downstream in tests
# so that messes up obj. Hence the quick function to rebuild it
# back to the MLE
TMB::runExample('simple')
obj0 <- obj
par0 <- obj$env$parList()

get_simple_obj <- function() {
  TMB::MakeADFun(data=obj0$env$data, parameters=par0, random=obj0$env$random,
                 dll=obj0$env$DLL, silent=TRUE)}

obj <- get_simple_obj()

### Skip all this if on CRAN. Otherwise locally or on CI, need to
### build the executables and run them so they're available for
### the tests. Then cleanup. On CRAN only a .RDS file is read in
### and really simple tests are performed.
if(!skip_ADMB){
if(Sys.getenv("NOT_CRAN")=='true'){
  oldwd <- getwd()
  setwd('../simple')
  system("admb simple", ignore.stdout = TRUE)
  system('./simple', ignore.stdout = TRUE)
  expect_equal(readLines('simple.par')[2], '# a:') # hack to test something
  dir.create('../simple_long_filename')
  trash <- file.copy('../simple/simple.tpl',
                     to='../simple_long_filename/simple_long_filename.tpl')
  trash <- file.copy('../simple/simple.dat',
                     to='../simple_long_filename/simple_long_filename.dat')
  setwd('../simple_long_filename')
  system("admb simple_long_filename", ignore.stdout = TRUE)
  system('./simple_long_filename', ignore.stdout = TRUE)
  setwd(oldwd)

  ## Clean up files to pass checks locally
  if(requireNamespace('withr')){
    withr::defer({
      files <- list.files('../simple', full.names = TRUE)
      ignore <- file.remove(files[-grep('.dat|.tpl', x=files)])
      ignore <- file.remove('../simple/mceval.dat')
      unlink('../simple_long_filename', TRUE)
      unlink("../simple_chain_1", TRUE)
      unlink("../simple_chain_2", TRUE)
      unlink("../simple_chain_3", TRUE)
      ## dev.off()
      ## plotout <- 'Rplots.pdf'
      ## trash <- if(file.exists(plotout)) file.remove(plotout)
    }, teardown_env())
  }

}
}
