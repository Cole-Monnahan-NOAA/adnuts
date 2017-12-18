# adnuts
An R package for NUTS sampling for ADMB and TMB models

For now, install using:

````devtools::install_github('colemonnahan/adnuts', build_vignettes=TRUE)````

The package contains a user guide: `vignette('adnuts')`. This will help get you started.

To try fitting your own models, see example code in the inst/demo.R file.

To use the ADMB functionality, you will first need to build the development
version of ADMB. Follow the ADMB instructions for how to compile ADMB from
source. Then, recompile your .tpl file with this modified version. Then you
model .exe file will contain the new MCMC functionality. Otherwise you will
not be able to use this package.

In the future this functionality will be included in the distributed ADMB
version, but for now it has not been merged in.
