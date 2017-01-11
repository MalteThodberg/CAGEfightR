# library(CAGEfightR)
# library(data.table)
#
# # List of files
# tmp <- list.files(path="~/Desktop/bigwigs/", full.names=TRUE)
# bwPlus <- grep(pattern=".plus.", x=tmp, value=TRUE)
# bwMinus <- grep(pattern=".minus.", x=tmp, value=TRUE)
# bedBoth <- list.files("~/Desktop/to_be_zipped/", full.names=TRUE)
#
# ### Common functions
#
# TPM <- function(x){
# 	x / (sum(x) / 1e6)
# }
#
# preFilter <- function(dt, minCTSS=0, minTPM=0){
# 	# CTSS filtering
# 	if(minCTSS > 0){
# 		# CTSS filter
# 		dt <- dt[score >= minCTSS,]
# 	}
#
# 	# TPM filtering
# 	if(minTPM > 0){
# 		# CTSS filter
# 		dt <- dt[TPM(score) >= minTPM,]
# 	}
#
# 	# TPM normalize
# 	dt[, score := TPM(score)]
#
# 	# Return
# 	dt
# }
#
# addTPM <- function(dt1, dt2){
# 	stopifnot(all(colnames(dt1) == colnames(dt2)))
#
# 	# Bind all data.frame
# 	dt <- rbindlist(list(dt1, dt2))
#
# 	# Merge without preprocess
# 	dt <- dt[,.(score=sum(score)), by=c("seqnames", "start", "end", "strand")]
#
# 	# Return
# 	dt
# }
#
# ### BigWig Functions
#
# fastReadBigWig <- function(fnamePlus, fnameMinus){
# 	# Read files
# 	countsPlus <- import(fnamePlus)
# 	strand(countsPlus) <- "+"
# 	countsMinus <- import(fnameMinus)
# 	strand(countsMinus) <- "-"
#
# 	# Merge
# 	dt <- rbindlist(list(as.data.table(countsPlus), as.data.table(countsMinus)))
# 	rm(countsPlus, countsMinus)
#
# 	# Remove unused width column to save mem
# 	dt[,width:=NULL]
#
# 	# Return
# 	dt
# }
#
# singleBigWig2TPM <- function(fnamePlus, fnameMinus, ...){
# 	# Read files
# 	dt <- fastReadBigWig(fnamePlus=fnamePlus, fnameMinus=fnameMinus)
#
# 	# Filter
# 	dt <- preFilter(dt, ...)
#
# 	# Return
# 	dt
# }
#
# multipleBigWig2TPM <- function(fnamesPlus, fnamesMinus, ...){
# 	# Initial level
# 	dt1 <- singleBigWig2TPM(fnamesPlus[1], fnamesMinus[1], ...)
#
# 	# Necessary for-loop over indices
# 	is <- seq_along(fnamesPlus[-1])
#
# 	for(i in is){
# 		dt2 <- singleBigWig2TPM(fnamesPlus[i], fnamesMinus[i], ...)
# 		dt1 <- addTPM(dt1, dt2)
# 	}
#
# 	# Return
# 	dt1
# }
#
# BigWig2TPM <- function(plus, minus, ...){
# 	stopifnot(length(plus) == length(minus))
#
# 	if(length(plus) == 1){
# 		message("Single...")
# 		dt <- singleBigWig2TPM(fnamePlus=plus, fnameMinus=minus, ...)
# 	}else{
# 		message("Multiple...")
# 		dt <- multipleBigWig2TPM(fnamesPlus=plus, fnamesMinus=minus, ...)
# 	}
#
# 	# Return
# 	dt
# }
#
# # Tests
# head(BigWig2TPM(plus=bwPlus[1], minus=bwMinus[1]))
# head(BigWig2TPM(plus=bwPlus[1:3], minus=bwMinus[1:3]))
#
# ### Bed functions
#
# singleBED2TPM <- function(fname, ...){
# 	# Read file
# 	dt <- fastReadBED(fname=fname)
#
# 	# Filter
# 	dt <- preFilter(dt, ...)
#
# 	# Return
# 	dt
# }
#
# multipleBED2TPM <- function(fnames, ...){
# 	# Initial level
# 	dt1 <- singleBED2TPM(fnames[1], ...)
#
# 	# Necessary for-loop over indices
# 	for(fname in fnames[-1]){
# 		dt2 <- singleBED2TPM(fname=fname, ...)
# 		dt1 <- addTPM(dt1, dt2)
# 	}
#
# 	# Return
# 	dt1
# }
#
# BED2TPM <- function(ctss, ...){
# 	if(length(ctss) == 1){
# 		message("Single...")
# 		dt <- singleBED2TPM(fname=ctss, ...)
# 	}else{
# 		message("Multiple...")
# 		dt <- multipleBED2TPM(fnames=ctss, ...)
# 	}
#
# 	# Return
# 	dt
# }
#
# head(BED2TPM(ctss=bedBoth[5]))
# head(BED2TPM(ctss=bedBoth[1:3]))
#
# # BED2TPM <- function(x, ...){
# # 	if(length(x) == 1){
# # 		# Read file
# # 		dt <- fastReadBED(fname)
# #
# # 		# Filter
# # 		dt <- preFilter(dt, ...)
# # 	}else{
# # 		Reduce(addBED, fname)
# # 	}
# #
# # 	# Return
# # 	dt
# # }
# #
# # addBED <- function(dt1, fname, ...){
# # 	# Load new file
# # 	dt2 <- BED2TPM(fname, ...)
# #
# # 	# Merge
# # 	dt <- addTPM(dt1, dt2)
# #
# # 	# Return
# # 	dt
# # }
# #
# # addBigWig <- function(dt1, fnamePlus, fnameMinus, ...){
# # 	# Load new file
# # 	dt2 <- BigWig2TPM(fnamePlus, fnameMinus, ...)
# #
# # 	# Merge
# # 	dt <- addTPM(dt1, dt2)
# #
# # 	# Return
# # 	dt
# # }
#
#
# # Must add this as well
# # GRL2TPM <- function(grl, preFilterCTSS, preFilterTPM){
# #
# # }
#
# #
# # BEDsChunkedTPM <- function(dt1, fname, ...){
# # 	# Load new file
# # 	dt2 <- BED2TPM(fname, ...)
# #
# # 	# Merge
# # 	dt <- addTPM(dt1, dt2)
# #
# # 	# Return
# # 	dt
# # }
# #
# # runningBED2TPM <- function(fnames){
# #
# # }
# #
# # coverageOfCTSS <- function(ctss, ctssMinus, preFilterCTSS=0, preFilterTPM=0){
# # 	# Determine extension
# #
# # }
# #
# #
# #
# # d1 <- BigWig2TPM(bwPlus[[1]], bwMinus[[1]], minCTSS=2)
# # d2 <- BigWig2TPM(bwPlus[[2]], bwMinus[[2]], minCTSS=2)
# #
# # addTPM(d1, d2)
# #
# # multiCov <- Map(BigWig2TPM, bwPlus, bwMinus, minCTSS=2)
# #
# # chunkCoverage <-
# #
# # summedCov <- Reduce(mergeCoverage, multiCov)
#
#
