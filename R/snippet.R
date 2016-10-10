library(pbapply)
library(GenomicRanges)
library(GenomicFeatures)
library(data.table)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

source("R/IO.R")
source("R/coverage.R")
source("R/tag_clustering.R")
source("R/quantification.R")
source("R/annotation.R")

ctss_files <- list.files(path="~/Desktop/to_be_zipped/", full.names = TRUE)
grl <- importAsGRangesList(ctss_files)

highMem <- coverageOfCTSS(grl, ctssCutoff=1)
#lowMem <- coverageOfCTSS(ctss_files, ctssCutoff=1)

highTCs <- simpleTagClustering(highMem, tpmCutoff=0.1, mergeDist=20)
#lowTCs <- simpleTagClustering(lowMem, tpmCutoff=0.1, mergeDist=20)

highEM <- quantifyFeatures(ctss=grl, features=highTCs, cores=3)
#lowEM <- quantifyFeatures(ctss=ctss_files, features=highTCs)

highTCs$txType <- assignTxType(gr=highTCs, txdb=txdb)
highTCs$txID <- assignTxID(gr=highTCs, txdb=txdb)
highTCs$geneID <- assignGeneID(gr=highTCs, txdb=txdb)

GM <- sumByGene(em=highEM, genes=highTCs$geneID, keepUnannotated=TRUE)
