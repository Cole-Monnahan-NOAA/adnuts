library(testthat)
library(adnuts)

if(!.Platform$OS.type=='windows'){
  ## Swap out the Win executable for unix one
  file.rename('simple/simple.exe', 'simple/simple_windows.exe')
  file.copy('simple/simple_unix', to='simple/simple')
}
test_check("adnuts")
