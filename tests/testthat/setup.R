## Skip consistency and reproducibility tests? Only need to run
## these locally when ADMB changes.
skip_consistency <- TRUE
skip_reproducibility <- TRUE

### Skip all this if on CRAN. Otherwise locally or on CI, need to
### build the executables and run them so they're available for
### the tests. Then cleanup. On CRAN only a .RDS file is read in
### and really simple tests are performed.
if(Sys.getenv("NOT_CRAN")=='true'){
  oldwd <- getwd()
  setwd('../simple')
  system("admb simple")
  system('./simple')
  expect_equal(readLines('simple.par')[2], '# a:') # hack to test something
  dir.create('../simple_long_filename')
  trash <- file.copy('../simple/simple.tpl',
                     to='../simple_long_filename/simple_long_filename.tpl')
  trash <- file.copy('../simple/simple.dat',
                     to='../simple_long_filename/simple_long_filename.dat')
  setwd('../simple_long_filename')
  system("admb simple_long_filename")
  system('./simple_long_filename')
  setwd(oldwd)

  ## Clean up files to pass checks locally
  if(requireNamespace('withr')){
    withr::defer({
      files <- list.files('../simple', full.names = TRUE)
      ignore <- file.remove(files[-grep('.dat|.tpl', x=files)])
      unlink('../simple_long_filename', TRUE)
      unlink("../simple_chain_1", TRUE)
      unlink("../simple_chain_2", TRUE)
      unlink("../simple_chain_3", TRUE)
      plotout <- 'Rplots.pdf'
      trash <- if(file.exists(plotout)) file.remove(plotout)
    }, teardown_env())
  }

}
