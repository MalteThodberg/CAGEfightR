# CAGEfightR

CAGEfightR is an R-package for analyzing Cap Analysis of Gene Expression (CAGE) data in Bioconductor. 

## Installation

You can install CAGEfightR from github with `devtools`:

``` r
# Install Bioconductor packages:
BiocManager::install(c("GenomicRanges", "rtracklayer", "SummarizedExperiment", 
"BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicFeatures", 
"TxDb.Mmusculus.UCSC.mm9.knownGene", "org.Mm.eg.db,", "BiocParallel", 
"GenomicFiles", "Gviz"))

# Install CRAN packages and CAGEfightR
devtools::install_github("MalteThodberg/CAGEfightR")
```

## Examples

See the vignette for an in-depth guide to using `CAGEfightR`!
