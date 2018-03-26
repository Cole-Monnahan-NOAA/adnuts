
# adnuts 1.0.9000

* Update ADMB algorithms to use "-maxfn 0 -phase 1000" instead of "-noest". This helps with Stock Synthesis and likely other models where some initialization is skipped with -noest which can lead to unusual and undesirable behavior. Also change behavior with inits=NULL to pull MLE values from the admodel.hes file instead of putting to the .par file for inits. This fixes some models when negative phases are used.


 ------------------------------------------------------------------------
adnuts 1.0.0 (2018-02-04)
------------------------------------------------------------------------

Initial release.
