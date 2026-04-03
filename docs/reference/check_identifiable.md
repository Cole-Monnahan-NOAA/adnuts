# Check identifiability from model Hessian

Check identifiability from model Hessian

## Usage

``` r
check_identifiable(model, path = getwd())
```

## Arguments

- model:

  Model name without file extension

- path:

  Path to model folder, defaults to working directory

## Value

Prints output of bad parameters and invisibly returns it.

## Details

Read in the admodel.hes file and check the eigenvalues to determine
which parameters are not identifiable and thus cause the Hessian to be
non-invertible. Use this to identify which parameters are problematic.
This function was converted from a version in the `FishStatsUtils`
package.
