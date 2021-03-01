## Skip consistency and reproducibility tests? Only need to run
## these locally when ADMB changes.
skip_consistency <- FALSE
skip_reproducibility <- FALSE

## Setup executable whether testing locally (Windows) or through CI (linux 18.04).
oldwd <- getwd()
setwd('../simple')
# system("admb simple")
system('./simple')
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
  }, teardown_env())
}
