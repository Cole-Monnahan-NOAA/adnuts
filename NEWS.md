------------------------------------------------------------------------
adnuts 1.0.1.9000  
------------------------------------------------------------------------

* Rework metric options to allow user to access ADMB 13.0's new
  dense mass matrix adaptation scheme. Added new section
  demonstrating this in the vignette.
* Add more control to 'sample_admb' via 'skip_monitor',
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
* Updated pairs_admb to have an 'order' argument for quickly
  plotting slow/fast parameters
* Add new function plot_marginals for quickly plotting posterior
  histograms
* Add new function plot_sampler_params to plot NUTS sampling  
* Make parallel the default and deprecate the 'parallel'
  argument.
* Fix bug in parallel path which failed when it was absolute. Now
  can be relative or absolute. 
* Add check for valid version of ADMB
* Minor bug fixes and documentation updates

------------------------------------------------------------------------
adnuts 1.0.1 (2019-03-15) 
------------------------------------------------------------------------

* Update ADMB algorithms to use "-maxfn 0 -phase 1000" instead of "-noest". This helps with Stock Synthesis and likely other models where some initialization is skipped with -noest which can lead to unusual and undesirable behavior. Also changed behavior with inits=NULL to pull MLE values from the admodel.hes file instead of pulling from the .par file for inits. This fixes some models when negative phases are used.

* Add function check_identifiable which examines a .hes file and reports which parameters are not well identified.

* Add function sample_inits to generate inits from a previous fitted object.

* Read in MLE values from the .hes file when inits=NULL, instead of from the .par file. 

* Add informative errors for common issues.

* Minor bug fixes and updates.


------------------------------------------------------------------------
adnuts 1.0.0 (2018-02-04)
------------------------------------------------------------------------

Initial release.
