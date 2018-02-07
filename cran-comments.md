## Resubmission
This is a resubmission. In this version I have:

* Added URL to key paper cited in description.
* Edited the description to be more concise.
* Added new examples to important top level functions. In some cases this required providing a .RDS file which was the output of top level functions. This is due to the need to externally compile C++ models which can take >5s.
* Converted existing examples which were "dontrun" to be run.

## Test environments
* local Windows 10, R 3.4.3
* ubuntu 14.04 (on travis-ci), R 3.4.2
* macOS High Sierra 10.13.3, R 3.4.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. There was one NOTE because this is a new submission. 

## Downstream dependencies
There are currently no downstream dependencies for this package.
