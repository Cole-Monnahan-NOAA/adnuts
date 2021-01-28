library(testthat)
library(adnuts)

if(!.Platform$OS.type=='windows'){
  ## Swap out the Win executable for unix one
  if(!file.exists('../simple/simple_unix'))
    stop("Failed to copy Unix executable b/c not found")
  file.rename('../simple/simple.exe', '../simple/simple_windows.exe')
  file.copy('../simple/simple_unix', to='../simple/simple')
}
test_check("adnuts")
