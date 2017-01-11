# library(CAGEfightR)
# library(GenomicFiles)
#
# # List of files
# tmp <- list.files(path="~/Desktop/bigwigs/", full.names=TRUE)
# tmp <- rep(tmp, 20)
# bwPlus <- BigWigFileList(grep(pattern=".plus.", x=tmp, value=TRUE))
# bwMinus <- BigWigFileList(grep(pattern=".minus.", x=tmp, value=TRUE))
# bedBoth <- list.files("~/Desktop/to_be_zipped/", full.names=TRUE)
#
# fastReadBigWig <- function(fnamePlus, fnameMinus){
# 	# Read files
# 	countsPlus <- import(fnamePlus)
# 	strand(countsPlus) <- "+"
# 	countsMinus <- import(fnameMinus)
# 	strand(countsMinus) <- "-"
#
# 	# To GR
# 	gr <- c(countsPlus, countsMinus)
#
# 	gr
# }
#
# TPM <- function(x){
# 	x / (sum(x) / 1e6)
# }
# BigWig2TPM <- function(fnamePlus, fnameMinus){
# 	gr <- fastReadBigWig(fnamePlus, fnameMinus)
# 	score(gr) <- TPM(score(gr))
# 	gr
# }
#
# x <- mapply(BigWig2TPM, bwPlus, bwMinus)
# x <- GRangesList(x)
#
# # Get covs
# countsPlus <- lapply(bwPlus, import, as="RleList")
# countsMinus <- lapply(bwMinus, import, as="RleList")
#
# # Calculate library sizes
# libSizes <- colSums(sapply(countsPlus, sum) + sapply(countsMinus, sum)) / 1e6
#
# # Calculate TPMs
# tpmPlus <- Map(function(x, y) x / y, countsPlus, libSizes)
# tpmMinus <- Map(function(x, y) x / y, countsMinus, libSizes)
#
# # Total Coverage
# sumPlus <- Reduce(function(x1, x2) x1 + x2, tpmPlus)
# sumMinus <- Reduce(function(x1, x2) x1 + x2, tpmMinus)
#
# # Merge to GRanges
# gr <- c(GRanges(sumPlus, strand="+"), GRanges(sumMinus, strand="-"))
#
# # Remove Zeroes
# gr <- subset(gr, score > 0)
#
# # Sort
# sort(gr)
#
# # TCs
# tcs <- simpleTagClustering(gr, tpmCutoff=1)
# sum(width(tcs) > 1000)
#
# # Quantify
# features_all <- BigWigSelection(tcs)
#
# quantifyBigWig <- function(bwPlus, bwMinus, features){
# 	# BigWig Selection
# 	features <- BigWigSelection(features)
#
# 	# Count via direct access
# 	vecPlus <- sum(import(bwPlus, selection=features_all, as="NumericList"))
# 	names(vecPlus) <- NULL
# 	vecMinus <- sum(import(bwMinus, selection=features_all, as="NumericList"))
# 	names(vecMinus) <- NULL
#
# 	# Black list strand
# 	vecPlus <- ifelse(strand(tcs) == "-", 0, vecPlus)
# 	vecMinus <- ifelse(strand(tcs) == "+", 0, vecMinus)
#
# 	# Merge strand
# 	vecBoth <- vecPlus + vecMinus
#
# 	# Return
# 	vecBoth
# }
#
# fastReadBigWig <- function(fnamePlus, fnameMinus){
# 	# Read files
# 	countsPlus <- import(fnamePlus)
# 	strand(countsPlus) <- "+"
# 	countsMinus <- import(fnameMinus)
# 	strand(countsMinus) <- "-"
#
# 	# To GR
# 	gr <- c(countsPlus, countsMinus)
#
# 	gr
# }
#
# quantifyGR <- function(bwPlus, bwMinus, features){
# 	# Load as GR
# 	gr <- fastReadBigWig(fnamePlus=bwPlus, fnameMinus=bwMinus)
#
# 	# Count overlaps
# 	o <- countScoredOverlaps(gr=gr, features=features)
#
# 	# return
# 	o
# }
#
# #quantifyBigWig(bwPlus[[1]], bwMinus[[1]], features=tcs)
# #quantifyGR(bwPlus[[1]], bwMinus[[1]], features=tcs)
#
#
#
# mapply(FUN=quantifyGR, bwPlus, bwMinus, list(features=tcs))
