# Check that the model is compiled with the right version of ADMB which is 12.0 or later

Check that the model is compiled with the right version of ADMB which is
12.0 or later

## Usage

``` r
.check_ADMB_version(model, path = getwd(), min.version = 12, warn = TRUE)
```

## Arguments

- model:

  Model name without file extension

- path:

  Path to model folder, defaults to working directory. NULL value
  specifies working directory (default).

- min.version:

  Minimum valid version (numeric). Defaults to 12.0.

- warn:

  Boolean whether to throw warnings or not

## Value

Nothing, errors out if either model could not be run or the version is
incompatible. If compatible nothing happens.

## Details

Some functionality of packages adnuts is imbedded in the ADMB source
code so that when a model is compiled it is contained in the model
executable. If this code does not exist adnuts will fail. The solution
is to update ADMB and recompile the model.
