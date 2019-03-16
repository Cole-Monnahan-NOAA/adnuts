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
