# library(GenomicRanges)
# library(rtracklayer)
#
# # List of files
# tmp <- list.files(path="~/Desktop/bigwigs/", full.names=TRUE)
# tmp <- rep(tmp, 10)
# bwPlus <- BigWigFileList(grep(pattern=".plus.", x=tmp, value=TRUE))
# bwMinus <- BigWigFileList(grep(pattern=".minus.", x=tmp, value=TRUE))
#
# # Get lib sizes
# libSizesPlus <- sapply(bwPlus, function(x) sum(score(import.bw(x))))
# libSizesMinus <- sapply(bwMinus, function(x) sum(score(import.bw(x))))
# libSizes <- (libSizesPlus + libSizesMinus) / 1e6
#
# # Load into memory
# plus <- lapply(bwPlus, import, as="RleList")
# minus <- lapply(bwMinus, import, as="RleList")
#
# # TPM normalize
# plus <- Map(f=function(rl, i) rl / i, plus, libSizes)
# minus <- Map(f=function(rl, i) rl / i, minus, libSizes)
#
# # Sum
# plus <- Reduce(f=function(x1, x2) x1+x2, plus)
# minus <- Reduce(f=function(x1, x2) x1+x2, minus)
#
# # Smooth
# smoothCoverage <- function(rl, f=runmean, ...){
# 	rl - f(rl, ...)
# }
#
# width=801
# plus <-	smoothCoverage(plus, k=width, endrule="constant")
# minus <-	smoothCoverage(minus, k=width, endrule="constant")
#
# # Truncate
# truncateCoverage <- function(rle, x){
# 	rle[rle < x] <- x
# 	rle
# }
#
# truncateCoverage <- function(rle, x){
# 	o <- GRanges(rle)
# 	o <- subset(o, score >= x)
# 	o <- coverage(o, weight="score")
# 	o
# }
#
# # plus <- truncateCoverage(plus, 0)
# # minus <- truncateCoverage(minus, 0)
#
# # To Coverage
# gr <- c(GRanges(plus, strand="+"),
# 				GRanges(minus, strand="-"))
# gr <- subset(gr, score > 0)
