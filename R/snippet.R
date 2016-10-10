library(pbapply)
library(GenomicRanges)
library(data.table)

source("R/IO.R")
source("R/coverage.R")
source("R/tag_clustering.R")
source("R/quantification.R")

ctss_files <- list.files(path="~/Desktop/to_be_zipped/", full.names = TRUE)
grl <- importAsGRangesList(ctss_files)

highMem <- coverageOfCTSS(grl, ctssCutoff=1)
lowMem <- coverageOfCTSS(ctss_files, ctssCutoff=1)

highTCs <- simpleTagClustering(highMem, tpmCutoff=0.01, mergeDist=25)
lowTCs <- simpleTagClustering(lowMem, tpmCutoff=0.01, mergeDist=25)

highEM <- quantifyFeatures(ctss=grl, features=highTCs)
lowEM <- quantifyFeatures(ctss=ctss_files, features=lowTCs)

