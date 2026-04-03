# Check if the session is interactive or Rstudio which has implications for parallel output

Check if the session is interactive or Rstudio which has implications
for parallel output

## Usage

``` r
.check_console_printing(parallel)
```

## Arguments

- parallel:

  Boolean whether chain is executed in parallel mode or not.

## Value

Boolean whether output should be printed to console progressively, or
saved to file and printed at the end.

## Details

When using RStudio and RGui, the parallel output does not show on the
console. As a workaround it is captured in each cluster into a file and
then read in and printed.
