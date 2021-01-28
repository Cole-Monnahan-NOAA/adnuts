
## For now I'm providing copies of executables and just swapping
## out. Only temporary and really bad form.
oldwd <- getwd()
setwd('../simple')
if(!.Platform$OS.type=='windows'){
  ## Swap out the Win executable for unix one
  if(!file.exists('../simple/simple_unix'))
    stop("Failed to copy Unix executable b/c not found")
  file.copy('../simple/simple_unix', to='../simple/simple')
  system('sudo chmod a+x simple') # give permission
  system('./simple -nox')
} else {
  if(!file.exists('../simple/simple_windows.exe'))
    stop("Failed to copy Windows executable b/c not found")
  file.copy('../simple/simple_windows.exe',
            '../simple/simple.exe')
  system("simple -nox")
}
setwd(oldwd)
