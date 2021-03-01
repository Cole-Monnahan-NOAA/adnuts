------------------------------------------------------------------------
adnuts 1.1.2 (2021-03-01get)
------------------------------------------------------------------------

* Improve console output for RStudio users. Was broken for NUTS chains and in serial.

* Add new option -verbose which suppresses almost all console output when set to FALSE

* Update demo file, vignette and README in preparation for submission to CRAN

* Add new function plot_uncertainties

* Expand continuous testing

* Add slot 'par_names' to objects of type adfit

------------------------------------------------------------------------
adnuts 1.1.1 (2021-02-19)
------------------------------------------------------------------------

* Add slot for par_names to adfit objects

* Add method `as.data.frame` for class `adfit`

* Improved and expanded testing via continuous integration

* Print the ADMB command to console when it fails to run properly to help user diagnose issues

* Improve console output for RStudio users. It will now print at conclusion of parallel runs.

* Fix bugs in model names for MacOS (use ./model insted of model internally)

* Fix small bug with mceval=TRUE for newest verison of stock synthesis

* Fix `sample_tmb` to work again for short-term use

------------------------------------------------------------------------
adnuts 1.1.0 (2020-07-13)
------------------------------------------------------------------------

* Change from `sample_admb` to `sample_nuts` and `sample_rwm` to
  run the NUTS and RWM algoriths, respectively.
  
* Rework metric options to allow user to access ADMB 12.2's new
  dense mass matrix adaptation scheme. Added new section
  demonstrating this in the vignette.
  
* Add more control to via 'skip_monitor',
  'skip_unbounded', and 'skip_optimization' arguments
  
* Remove TMB references from documentation and vignette, 
  instead pointing users to package 'tmbstan', and collate
  deprecated R code into a single file
  
* Migrate to new github repo: github.com/Cole-Monnahan-NOAA per
  NOAA's policy
  
* Add testing via testthat package

* Turn on calculation of ESS and Rhat manually, which get used in
  subsequent functions
  
* Created S3 class 'adfit' and generic methods print, summary,
  and plot
  
* Updated `pairs_admb` to have an 'order' argument for quickly
  plotting slow/fast parameters
  
* Add new function `plot_marginals` for quickly plotting posterior
  histograms
  
* Add new function `plot_sampler_params` to plot NUTS sampling

* Make parallel the default and deprecate the 'parallel'
  argument.
  
* Fix bug in parallel path which failed when it was absolute. Now
  can be relative or absolute. 
  
* Add check for valid version of ADMB

* Minor bug fixes and documentation updates

* Improve error handling and testing routines

------------------------------------------------------------------------
adnuts 1.0.1 (2019-03-15) 
------------------------------------------------------------------------

* Update ADMB algorithms to use "-maxfn 0 -phase 1000" instead of
  "-noest". This helps with Stock Synthesis and likely other
  models where some initialization is skipped with -noest which
  can lead to unusual and undesirable behavior. Also changed
  behavior with inits=NULL to pull MLE values from the
  admodel.hes file instead of pulling from the .par file for
  inits. This fixes some models when negative phases are used.

* Add function check_identifiable which examines a .hes file and
  reports which parameters are not well identified.

* Add function sample_inits to generate inits from a previous
  fitted object.

* Read in MLE values from the .hes file when inits=NULL, instead
  of from the .par file.

* Add informative errors for common issues.

* Minor bug fixes and updates.


------------------------------------------------------------------------
adnuts 1.0.0 (2018-02-04)
------------------------------------------------------------------------

Initial release.
