# library(data.table)
# library(GenomicRanges)
# library(rtracklayer)
# library(data.table)
# library(microbenchmark)
# library(BiocParallel)
#
# ### Class for pointing to files
#
# CTSSFile <- function(fname=NULL, fnamePlus=NULL, fnameMinus=NULL){
# 	# Create empty class
# 	o <- structure(list(combinedStrand=NULL,
# 											fileType=NULL, path=NULL,
# 											pathPlus=NULL, pathMinus=NULL),
# 								 class = "CTSSFile")
#
# 	supportedFileTypes <- c("bed", "bw")
#
# 	# Determine type of input
# 	combinedStrand <- !is.null(fname)
# 	splitStrand <- !is.null(fnamePlus) | !is.null(fnameMinus)
# 	stopifnot(sum(c(combinedStrand, splitStrand)) < 2)
#
# 	if(combinedStrand){
# 		# Check file for extension
# 		stopifnot(file.exists(fname))
# 		fileExt <- tools::file_ext(fname)
# 		stopifnot(fileExt %in% supportedFileTypes)
#
# 		# Set data
# 		o$combinedStrand <- TRUE
# 		o$fileType <- ".bed"
# 		o$path <- fname
# 	}else if(splitStrand){
# 		# Check file for extension
# 		# Add checks here
#
# 		# Set data
# 		o$combinedStrand <- TRUE
# 		o$fileType <- "bw"
# 		o$pathPlus <- fnamePlus
# 		o$pathMinus <- fnameMinus
# 	}else{
# 		warning("Returning empty CTSS file!")
# 	}
#
# 	# Return
# 	o
# }
#
# print.CTSSFile <- function(obj){
# 	if(is.null(obj$combinedStrand)){
# 		o <- "Empty CTSSFile"
# 	}else if(obj$combinedStrand){
# 		o <- paste0("CTSS data as ", obj$fileType)
# 	}else if(!obj$combinedStrand){
# 		o <- paste0("CTSS data as ", obj$fileType)
# 	}
# 	print(o)
# }
#
# ### Files
#
# # Tester files
# tmp <- list.files(path="~/Desktop/bigwigs/", full.names=TRUE)
# bwPlus <- grep(pattern=".plus.", x=tmp, value=TRUE)
# bwMinus <- grep(pattern=".minus.", x=tmp, value=TRUE)
# bedBoth <- list.files("~/Desktop/to_be_zipped/", full.names=TRUE)
# bedBoth <- rep(bedBoth, 20)
#
#
# ### Functions
#
# fastReadBED <- function(fname){
# 	fread(input=fname,
# 										data.table=TRUE,
# 										sep="\t",
# 										col.names=c("seqnames", "start", "end", "score", "strand"),
# 										select=c(1,2,3,5,6),
# 										verbose=FALSE,
# 										showProgress=FALSE,
# 										na.strings=NULL)
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
# fastReadTabixBED <- function(fname){
# 	import(fname)
# }
#
# DT2GR <- function(dt){
# 	GRanges(seqnames=Rle(dt$seqnames),
# 							 ranges=IRanges(dt$start+1, width=1),
# 							 strand=Rle(dt$strand),
# 							 score=dt$score)
# }
#
#
# BED2GR <- function(fname){
# 	DT2GR(fastReadBED(fname))
# }
#
# # microbenchmark(BED2GR(bedBoth[1]), times = 10)
# # microbenchmark(fastReadBigWig(bwPlus[[1]], bwMinus[[1]]), times=10)
# #
# # microbenchmark(lapply(bedBoth, BED2GR), times=10)
# # microbenchmark(mapply(fastReadBigWig, bwPlus, bwMinus), times=10)
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
# # Outer reduce loop
# stackDT <- function(dt1, dt2, ...){
# 	# TPM norm
# 	dt2 <- preFilter(dt2, ...)
#
# 	# Stack
# 	dt3 <- rbindlist(list(dt1, dt2), use.names=TRUE)
#
# 	# Coverage
# 	dt3 <- dt3[, .(score=sum(score)), by = c("seqnames", "start", "end", "strand")]
#
# 	# Return
# 	dt3
# }
#
# # Inner reduce loop
# stackFile <- function(dt1, fname, ...){
# 	# Read file
# 	dt2 <- fastReadBED(fname)
#
# 	# Call higher level
# 	dt3 <- stackDT(dt1, dt2, ...)
#
# 	# Return
# 	dt3
# }
#
# stackChunk <- function(fnames, ...){
# 	# Empty data.table for init
# 	et <- data.table(seqnames=character(),
# 									 start=numeric(),
# 									 end=numeric(),
# 									 score=numeric(),
# 									 strand=character())
#
# 	# Reduce with pass down args
# 	o <- Reduce(f=function(x1, x2) stackFile(dt1=x1, fname=x2, ...),
# 				 x=fnames, init=et)
#
# 	o
# }
#
# # Calculate tpm
# TPM <- function(x){
# 	x / (sum(x) / 1e6)
# }
#
# collapseDT <- function(listOfDTs, minCTSS=0, minTPM=0){
# 	# Merge into one giant data.table
# 	o <- rbindlist(listOfDTs, idcol="sample")
#
# 	# No longer need this big object
# 	rm(listOfDTs)
#
# 	# Filter
# 	if(minCTSS > 0){
# 		o <- o[score >= minCTSS,]
# 	}
#
# 	# Calculate TPM
# 	o <- o[, tpm := TPM(score), by=sample]
#
# 	# TPM trimming
# 	if(minTPM > 0){
# 		o <- o[tpm >= minTPM,]
#
# 		# Redo TPM calculations
# 		o <- o[, tpm := TPM(score), by=sample]
# 	}
#
# 	# Calculate coverage
# 	o <- o[, .(score=sum(tpm)), by = c("seqnames", "start", "end", "strand")]
#
# 	# Return
# 	o
# }
#
# ### Code
#
# # Get Chunks
# # cores <- 3
# # chunks <- split(bedBoth, factor(rep_len(seq_len(cores), length.out=length(bedBoth))))
#
# # # Empty frame for stacking
# # et <- data.table(seqnames=character(),
# # 											 start=numeric(),
# # 											 end=numeric(),
# # 											 score=numeric(),
# # 											 strand=character())
# #
# # # Test on single chunks
# # Reduce(f=function(x1, x2) stackFile(dt1=x1, fname=x2, minCTSS=2),
# # 			 x=chunks[[1]],
# # 			 init=et)
#
# coverageOfCTSS <- function(ctss, preFilterCTSS=0, preFilterTPM=0, highMem=FALSE, biocParallel=SerialParam()){
# 	# Check Files
# 	# ...
#
# 	# Find number of workers
# 	nWorkers <- biocParallel$workers
#
# 	if(!highMem){
# 		if(nWorkers == 1){
# 			message("Low Memory, 1 worker")
# 			o <- stackChunk(ctss, minCTSS=preFilterCTSS, minTPM=preFilterTPM)
# 		}else{
# 			message(paste0("Low Memory, ", nWorkers,  " workers"))
#
# 			# Chunks
# 			chunks <- split(ctss, factor(rep_len(seq_len(nWorkers), length.out=length(ctss))))
# 			stopifnot(min(sapply(chunks, length)) > 1)
#
# 			# First Reduce
# 			message(paste0("Processing chunks..."))
# 			o <- bplapply(X=ctss, FUN=stackChunk,
# 															minCTSS=preFilterCTSS, minTPM=preFilterTPM,
# 															BPPARAM=biocParallel)
#
# 			# Second reduce
# 			message(paste0("Merging chunks..."))
# 			o <- Reduce(f=stackDT, o)
# 		}
# 	}else if(highMem){
# 		message(paste0("High Memory, ", nWorkers,  " workers..."))
#
# 		# Load all BED files
# 		message("Reading all files...")
# 		o <- bplapply(X=ctss, FUN=fastReadBED, BPPARAM=biocParallel)
#
# 		# Collapse into coverage
# 		message("Processing all files...")
# 		o <- collapseDT(o, minCTSS=preFilterCTSS, minTPM=preFilterTPM)
# 	}
#
# 	# Output as GRanges
# 	message("Preparing GR")
# 	o <- DT2GR(o)
# 	o <- sort(o)
#
# 	# Return
# 	o
# }
#
# # Coverage
# #lowMem1 <- coverageOfCTSS(ctss=bedBoth)
# #lowMem3 <- coverageOfCTSS(ctss=bedBoth, biocParallel=MulticoreParam(workers=3))
# #highMem1 <- coverageOfCTSS(ctss=bedBoth, highMem=TRUE)
# highMem3 <- coverageOfCTSS(ctss=bedBoth, highMem=TRUE, biocParallel=MulticoreParam(workers=3))
# gc()
#
# # Tag clusters
# TCs <- CAGEfightR::simpleTagClustering(ctssCoverage=highMem3, tpmCutoff=1)
# tcs <- subset(TCs, score >= 3 & width > 1)
#
# # Quantify
# CM <- CAGEfightR::quantifyFeatures(ctss=bedBoth, features=tcs, cores=3)
#
# #
# # bplapply(X=bedBoth, FUN=stackChunk)
# #
# # # Serial on multiple chunks
# # firstReduce <- lapply(X=chunks, FUN=stackChunk, minCTSS=2)
# # secondReduce <- Reduce(f=stackDT, firstReduce)
# #
# # head(secondReduce)
# # # Empty data.table
# #
# # Reduce(stackFile, bedBoth[1:5], init=et)
# #
# # # Load all tables
# # d <- lapply(bedBoth[1:5], fastReadBED)
# # gc()
# # # Merge into one
# # d <- rbindlist(d, idcol="sample")
# # gc()
# # # Calculate tpm
# # TPM <- function(x){
# # 	x / (sum(x) / 1e6)
# # }
# #
# # d <- d[, tpm := TPM(score), by=sample]
# # gc()
# # # Summarize cov
# # d <- d[, .(score=sum(tpm)), by = c("seqnames", "start", "end", "strand")]
# # gc()
# # # Back to gr
# # DT2GR(d)
# # gc()
