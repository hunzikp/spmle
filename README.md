spmle
================
Philipp Hunziker
June 21, 2020

# System Requirements

To use the (optional but recommended) Python speed-up, you need Python
3.6+ (and [conda](https://conda.io/docs/user-guide/install/windows.html)
if on Windows).

# Installation

Install `spmle` directly from this repository:

``` r
devtools::install_github("hunzikp/spmle")
```

In case you don’t have scipy and numpy installed (or don’t know what
these are) and you want to use the Python speed-up, run the following in
an R console:

``` r
spmle::install_pydep()
```

This command installs the Python dependencies in the `r-reticulate`
virtual environment.

# Usage

To use the Python speed-up you need to specify a python binary or
virtual environment where scipy and numpy are installed. If you ran the
`install_pydep()` command above, then you can use the following:

``` r
library(reticulate)
library(spmle)
reticulate::use_virtualenv("r-reticulate")  # In Unix
reticulate::use_condaenv("r-reticulate")  # In Windows
```

# Examples

See [here](scripts/stmodel_testing.md).
