
## Setup executable whether testing locally (Windows) or through CI (linux 18.04).
oldwd <- getwd()
setwd('../simple')
system("admb simple")
system('simple')
setwd(oldwd)
## Clean up files to pass checks locally
if(requireNamespace('withr')){
  withr::defer({
    files <- list.files('../simple', full.names = TRUE)
    ignore <- file.remove(files[-grep('.dat|.tpl', x=files)])
    unlink("../simple_chain_1", TRUE)
    unlink("../simple_chain_2", TRUE)
    unlink("../simple_chain_3", TRUE)
  }, teardown_env())
}
