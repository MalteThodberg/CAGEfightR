# CAGEfightR

CAGEfightR is an R-package for analyzing Cap Analysis of Gene Expression (CAGE) data in Bioconductor. 

## Installation

Install the most recent stable version from Bioconductor:

```{r bioconductor, eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("CAGEfightR")
```

Or install the development version directly from GitHub using `devtools`:

``` r
# Install CRAN packages and CAGEfightR
devtools::install_github("MalteThodberg/CAGEfightR", build_vignettes = TRUE)
```

When using `devtools` you might need to manually install dependencies from Bioconductor (See `DESRIPTION` for a list of these).

## Examples

See the vignette for an in-depth guide to using `CAGEfightR`!
