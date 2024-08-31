# MITIGRE

<img src="img/mitigre.jpg" width="200" />

## Description

MITIGRE is an R package developed by Guillaume Cinkus and Vincent Marc at Avignon Université. It contains one function that launches an R shiny application for modelling isotopic tracers.

## Installation

MITIGRE requires an installation of R `4.0.0`. The instruction for the installation and the download of R can be found on the [CRAN website](https://cran.r-project.org/).

Once R is installed, MITIGRE can be installed from [GitHub](https://github.com/busemorose/MITIGRE).

``` r
if (!require("remotes")) install.packages("remotes") # install remotes package if needed
remotes::install_github("busemorose/MITIGRE") # install MITIGRE package
```
## Launch

Once the package is installed, the application can be loaded with the `MITIGRE()` function.

``` r
library(MITIGRE)
MITIGRE()
```

or

```
MITIGRE::MITIGRE()
```

## Data import

The data import field allows to load a file with the following conditions:

-   The data must be a plain text file
-   The file must have two columns representing precipitation and observed discharge, respectively named `obs` and `Cin`
-   The columns must be semicolon-separated and of the same length
-   Decimal separator is a dot `.`
