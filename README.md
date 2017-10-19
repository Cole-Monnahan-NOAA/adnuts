# adnuts
An R package for NUTS sampling using TMB and ADMB

This package is still in development, but is fairly stable. I will be
submiting to CRAN hopefully in Ocotober, 2017. For now, install using:

````devtools::install_github("colemonnahan/adnuts", ref="master")````

This code comes with no garuantees and will change over the next few
months.

To try fitting your own models, see example code in the inst/demo.R file.

To use the ADMB functionality, you will first need to use my forked version
of ADMB at https://github.com/colemonnahan/admb. Follow the ADMB
instructions for how to compile ADMB from source. Then, recompile your .tpl
file with this modified version. Then you model .exe file will contain the
new MCMC functionality. Otherwise you will not be able to use
`sample_admb`.

In the future this functionality will be included in the distributed ADMB
version, but for now it has not been merged in.
