# MITIGRE

<img src="img/mitigre2.png" width="200" />

Modèle Intégré de Transfert Isotopique Global des Ressources en Eau

## Description

MITIGRE is an R package developed by Guillaume Cinkus and Vincent Marc at Avignon Université. It contains one function that launches an R shiny application for modelling isotopic tracers (³H, ¹⁴C, ¹⁸O, ²H) with lumped parameter models (EMM, EPM, PEM, DM, PFM and mixtures of up to 4 components), calibrated by Monte Carlo. Two computation engines are available: **Forward** (convolution of the whole input series, corrected MITIGRE scheme) and **Backward** (stepping back in time from each sample date, TracerLPM method; Jurgens et al., 2012).

## Installation

MITIGRE requires an installation of R `4.0.0` or later. The instruction for the installation and the download of R can be found on the [CRAN website](https://cran.r-project.org/).

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
-   The file must have two columns representing input concentration and observed concentration, respectively named `Cin` and `obs`, and an optional `date` column (`YYYY-MM-DD`, `DD/MM/YYYY` or decimal year)
-   The columns must be semicolon-separated (comma or tab also accepted) and of the same length
-   Decimal separator is a dot `.`
-   Missing values (`NA`) are allowed in `obs`

Example files are provided with the package:

``` r
list.files(system.file("extdata", package = "MITIGRE"), full.names = TRUE)
```

## Reference

Jurgens, B.C., Böhlke, J.K., Eberts, S.M. (2012). TracerLPM (Version 1): An Excel® workbook for interpreting groundwater age distributions from environmental tracer data. U.S. Geological Survey Techniques and Methods 4-F3, 60 p. https://pubs.usgs.gov/publication/tm4F3
