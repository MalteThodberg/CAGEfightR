# #library(data.table)
# library(GenomicRanges)
# library(rtracklayer)
# #library(data.table)
# #library(microbenchmark)
# library(BiocParallel)
#
# ### Files
#
# # pombe bw files
# pombe_all <- list.files("~/example_data/pombe_bw/", full.names=TRUE)
# pombe_all <- grep(x=pombe_all, pattern="EMM|heat|nitro|H2O2", value=TRUE)
# pombe_plus <- BigWigFileList(grep(x=pombe_all, pattern=".plus.bw", value=TRUE))
# pombe_minus <- BigWigFileList(grep(x=pombe_all, pattern=".minus.bw", value=TRUE))
#
# # IBD
# ibd_all <- list.files("~/example_data/ibd_bw/", full.names=TRUE)
# ibd_all <- grep(x=ibd_all, pattern="hg19", value=TRUE)
# ibd_plus <- BigWigFileList(grep(x=ibd_all, pattern=".plus.bw", value=TRUE))
# ibd_minus <- BigWigFileList(grep(x=ibd_all, pattern=".minus.bw", value=TRUE))
#
# # Select
# bw_plus <- pombe_plus
# bw_minus <- pombe_minus
#
# ### BigWig manipulation
#
# # Check if genomes match
# checkBigWigGenomes <- function(bwPlus, bwMinus){
# 	# Get seqinfo
# 	seqInfoPlus <- lapply(bwPlus, seqinfo)
# 	seqInfoMinus <- lapply(bwMinus, seqinfo)
#
# 	# Check single seqinfo
# 	seqInfo <- unique(c(seqInfoPlus, seqInfoMinus))
#
# 	# Return the genome info
# 	seqInfo
# }
#
# # Merge different genomes to a common reference
# commonGenome <- function(bwPlus, bwMinus, method="intersect"){
# 	# Get seqinfo
# 	seqInfoPlus <- lapply(bwPlus, seqinfo)
# 	seqInfoMinus <- lapply(bwMinus, seqinfo)
#
# 	# Unique seqinfos
# 	seqInfo <- unique(c(seqInfoPlus, seqInfoMinus))
#
# 	# Sort every element
# 	seqInfo <- seqInfo[order(vapply(seqInfo, length, numeric(1)))]
#
# 	# Merge using either function
# 	if(method == "intersect"){
# 		o <- suppressWarnings(Reduce(f=intersect, seqInfo))
# 	}else if(method == "union"){
# 		o <- suppressWarnings(Reduce(f=merge, seqInfo))
# 	}
#
# 	# Sort
# 	o <- sortSeqlevels(o)
#
# 	# Return
# 	o
# }
#
# # Check if genome can be merged
# checkIfComparableGenomes <- function(bwPlus, bwMinus, seqInfo){
# 	# Get seqinfo
# 	seqInfoPlus <- lapply(bwPlus, seqinfo)
# 	seqInfoMinus <- lapply(bwMinus, seqinfo)
#
# 	# Unique seqinfos
# 	seqInfos <- unique(c(seqInfoPlus, seqInfoMinus))
#
# 	# Check if error is produces
# 	lapply(seqInfos, merge, y=seqInfo)
# }
#
#
# fastReadBigWig <- function(fnamePlus, fnameMinus, seqInfo=NULL){
# 	# Read files
# 	countsPlus <- import(fnamePlus)
# 	strand(countsPlus) <- "+"
# 	countsMinus <- import(fnameMinus)
# 	strand(countsMinus) <- "-"
#
# 	# To GR
# 	gr <- c(countsPlus, countsMinus)
#
# 	# Set new sequence info if necessary
# 	if(!is.null(seqInfo)){
# 		# Call error if genomes can't me merged with filename in question!
# 		chr_map <- match(seqlevels(seqInfo), seqlevels(gr))
# 		seqinfo(gr, new2old=chr_map, force=TRUE) <- seqInfo
# 	}
#
# 	gr
# }
#
# readCTSS <- function(bwPlus, bwMinus, genome=NULL, outFormat="GRangesList", biocParallel=bpparam(), ...){
# 	# Check genome if not provided
# 	if(is.null(genome)){
# 		# Check if all files have the same genome
# 		uniqueSeqInfos <- checkBigWigGenomes(bwPlus=bwPlus, bwMinus=bwMinus)
#
# 		if(length(uniqueSeqInfos) != 1){
# 			message("BigWig files do not have the same genome! Attempting to find common reference...")
# 			genome <- commonGenome(bwPlus=bwPlus, bwMinus=bwMinus, ...)
# 		}else{
# 			message("BigWig files have the same genome...")
# 			genome <- uniqueSeqInfos[[1]]
# 		}
#
# 	}else if(class(genome) == "character"){
# 		message("Retrieving genome from BSGenome...")
# 		genome <- SeqinfoForBSGenome(genome)
# 		checkIfComparableGenomes(bwPlus=bwPlus, bwMinus=bwMinus, seqInfo=genome)
#
# 	}else if(class(genome) == "Seqinfo"){
# 		message("Using supplied genome...")
# 		checkIfComparableGenomes(bwPlus=bwPlus, bwMinus=bwMinus, seqInfo=genome)
#
# 	}else{
# 		stop("Genome must be a Seqinfo-object or character!")
#
# 	}
#
# 	# Print used genome
# 	print(genome)
#
# 	#Read into R
# 	message("Reading in CTSS data from BigWig-files...")
# 	o <- bpmapply(FUN=fastReadBigWig, bwPlus, bwPlus,
# 								MoreArgs=list(seqInfo=genome),
# 								BPPARAM=biocParallel)
#
# 	if(outFormat == "GRangesList"){
# 		# Coerce to GRangesList
# 		message("Preparing final GRangesList...")
# 		o <- GRangesList(o)
# 	}else if(outFormat == "SimpleList"){
# 		message("Return list of GRanges...")
# 		o <- SimpleList(o)
# 	}
#
# 	# Return
# 	o
# }
#
# # TPM
# TPM <- function(x){
# 	x / (sum(x) / 1e6)
# }
#
# # Preprocessing function
# trimAndNormalizeGR <- function(gr, minCTSS=0, minTPM=0){
# 	# Filter on CTSS
# 	if(!is.null(minCTSS)){
# 		gr <- subset(gr, score >= minCTSS)
# 	}
#
# 	# Filter on TPM
# 	if(!is.null(minTPM)){
# 		gr <- subset(gr, TPM(score) >= minTPM)
# 	}
#
# 	# TPM normalize
# 	score(gr) <- TPM(score(gr))
#
# 	# Return
# 	gr
# }
#
# # Split by strad
# splitByStrand <- function(gr){
# 	split(gr, strand(gr))
# }
#
# smartUnlist <- function(grl){
# 	# Names
# 	if(is.null(names(grl))){
# 		i <- paste0("i", seq_len(length(grl)))
# 	}else{
# 		i <- names(grl)
# 		names(grl) <- NULL
# 	}
#
# 	# New name Rle-vector
# 	rleNames <- Rle(rep(i, elementNROWS(grl)))
#
# 	# Split
# 	gr <- unlist(grl)
# 	gr$i <- rleNames
#
# 	# Return
# 	gr
# }
#
# coverageOfCTSS <- function(grl, normalizationFunction=trimAndNormalizeGR, biocParallel=bpparam(), ...){
# 	# Checks
# 	# Should check if an apply loop produces only GRanges
#
# 	if(is.null(normalizationFunction)){
# 		message("Assuming GRangesList is already preprocessed")
# 	}else{
# 		# Preprocess
# 		message("Preprocessing...")
# 		grl <- bplapply(grl, normalizationFunction, ..., BPPARAM=biocParallel)
# 		grl <- GRangesList(grl)
# 		#grl <- endoapply(grl, normalizationFunction, ...)
# 	}
#
# 	# Flatten and split by strand
# 	message("Splitting...")
# 	grs <- splitByStrand(unlist(grl))
# 	rm(grl)
#
# 	# Plus cov
# 	message("Coverage...")
# 	plus <- coverage(grs$`+`, weight="score")
# 	minus <- coverage(grs$`-`, weight="score")
# 	rm(grs)
#
# 	# Merge as GR
# 	message("Output...")
# 	o <- c(GRanges(plus, strand="+"), GRanges(minus, strand="-"))
# 	rm(plus, minus)
#
# 	# Subset out zero ranges
# 	o <- subset(o, score > 0)
#
# 	# Return
# 	o
# }
#
# entropy.empirical
# freqs.empirical()
# entropy.plugin()
# ### Code
#
# # Check genome
# #si <- commonGenome(bwPlus=bw_plus, bwMinus=bw_minus, method="intersect")
#
# # Read in samples
# #register(SerialParam(), default=TRUE)
# register(MulticoreParam(workers=3), default=TRUE)
#
# ctssList <- readCTSS(bwPlus=bw_plus, bwMinus=bw_minus, outFormat="SimpleList")
#
# # Calculate coverage
# ctssCov <- coverageOfCTSS(grl=ctssList)
#
# # TCs
# tcs <- CAGEfightR::simpleTagClustering(ctssCov)
#
# # Quantify
# CM <- CAGEfightR::quantifyFeatures(ctss=ctssList, features=tcs)
#
#  rollmeanRle <- function (x, k){
#  	n <- length(x)
# 		cumsum(c(Rle(sum(window(x, 1, k))),
# 						 window(x, k + 1, n) - window(x, 1, n - k))) / k
# }
#
# ### Better preproc step
# #x <- smartUnlist(ctssList)
#
# # List of peaks
#
# shannonEntropy <-
#
# libSize <- function(gr){
# 	sum(score(gr))
# }
#
# approximateCounts <- function(gr, targetLibSize){
# 	# Normalize
# 	score(gr) <- score(gr) / sum(score(gr))
#
# 	# Approximate counts
# 	score(gr) <- round(score(gr) * targetLibSize)
#
# 	# Return
# 	gr
# }
#
# TCs <- splitByStrand(subset(tcs, width>=10 & score >= 3))
#
# tmp <- splitByStrand(ctssList[[1]])
# covplus <- coverage(tmp$`+`, weight="score")
# covminus <- coverage(tmp$`-`, weight="score")
#
# viewsPlus <- Views(covplus, as(TCs$`+`, "RangesList"))
# viewsPlus$I
#
# f1 <- function(x){rep(seq_along(x), x)}
# f2 <- function(y, z){f1(round(tpm*z))/z}
#
# ### Functions
# #
# # # Read BED
# # fastReadBED <- function(fname){
# # 	fread(input=fname,
# # 				data.table=TRUE,
# # 				sep="\t",
# # 				col.names=c("seqnames", "start", "end", "score", "strand"),
# # 				select=c(1,2,3,5,6),
# # 				verbose=FALSE,
# # 				showProgress=FALSE,
# # 				na.strings=NULL)
# # }
# #
# # # DT to gr
# # DT2GR <- function(dt, compr){
# # 	GRanges(seqnames=Rle(dt$seqnames),
# # 					ranges=IRanges(dt$start+1, width=1),
# # 					strand=Rle(dt$strand),
# # 					score=dt$score)
# # }
# #
# # DT2GR <- function(dt, compr){
# # 	GRanges(seqnames=Rle(dt$seqnames),
# # 					ranges=IRanges(dt$start+1, dt$end),
# # 					strand=Rle(dt$strand),
# # 					score=dt$score)
# # }
# #
# # # DT2GR <- function(dt, rleScore=FALSE, width=1){
# # # 	# Basic info
# # # 	tmp_seqnames <- Rle(dt$seqnames)
# # # 	tmp_strand <- Rle(dt$strand)
# # #
# # # 	# Build ranges
# # # 	if(is.null(width)){
# # # 		tmp_ranges <- IRanges(dt$start+1, end=dt$end)
# # # 	}else{
# # # 		tmp_ranges <- IRanges(dt$start+1, width=width)
# # # 	}
# # #
# # # 	# Build score
# # # 	if(rleScore){
# # # 		tmp_score <- Rle(dt$score)
# # # 	}else{
# # # 		tmp_score <- dt$score
# # # 	}
# # #
# # # 	# Build GRanges
# # # 	o <- GRanges(seqnames=tmp_seqnames,
# # # 					ranges=tmp_ranges,
# # # 					strand=tmp_strand,
# # # 					score=tmp_score)
# # #
# # # 	# Return
# # # 	o
# # # }
# #
# # TPM <- function(x){
# # 	x / (sum(x) / 1e6)
# # }
# #
# # collapseDT <- function(listOfDTs, minCTSS=0, minTPM=0){
# # 	# Merge into one giant data.table
# # 	o <- rbindlist(listOfDTs, idcol="sample")
# #
# # 	# No longer need this big object
# # 	rm(listOfDTs)
# #
# # 	# Filter
# # 	if(minCTSS > 0){
# # 		o <- o[score >= minCTSS,]
# # 	}
# #
# # 	# Calculate TPM
# # 	o <- o[, tpm := TPM(score), by=sample]
# #
# # 	# TPM trimming
# # 	if(minTPM > 0){
# # 		o <- o[tpm >= minTPM,]
# #
# # 		# Redo TPM calculations
# # 		o <- o[, tpm := TPM(score), by=sample]
# # 	}
# #
# # 	# Calculate coverage
# # 	o <- o[, .(score=sum(tpm)), by = c("seqnames", "start", "end", "strand")]
# #
# # 	# Return
# # 	o
# # }
# #
# # quantifyDT <- function(dt, features){
# # 	# DT to GR
# # 	gr <- DT2GR(dt)
# #
# # 	# Solution from bioconductor by M Morgan
# # 	hits <- as(findOverlaps(query=features, subject=gr), "List")
# # 	weightedCount <- sum(extractList(score(gr), hits))
# #
# # 	# Return
# # 	weightedCount
# # }
# #
# # mergeReadCTSS <- function(fnames, cores=3){
# # 	# Read list of GR
# # 	o <- mclapply(fnames, function(x) DT2GR(fastReadBED(x), rleScore=TRUE), mc.cores=cores)
# #
# # 	# Merge into GRangesList
# # 	dtList <- GRangesList(o)
# #
# # 	# Set names
# # 	names(o) <- fnames
# #
# # 	dtList
# # }
# #
# # splitReadCTSS <- function(fnames, cores=3){
# # 	# Read list of DTs
# # 	o <- mclapply(fnames, fastReadBED, mc.cores=cores)
# #
# # 	# Merge into DT
# # 	o <- rbindlist(o, idcol="sample")
# #
# # 	# To GR
# # 	o <- GRanges(seqnames=Rle(o$seqnames),
# # 							 ranges=IRanges(o$start+1, width=1),
# # 							 strand=Rle(o$strand),
# # 							 score=Rle(o$score),
# # 							 sample=Rle(o$sample))
# #
# # 	# Split into GRangesList
# # 	o <- split(o, o$sample)
# # 	names(o) <- basename(fnames)
# #
# # 	o
# # }
# #
# # # Split by strad
# # splitByStrand <- function(gr){
# # 	split(gr, strand(gr))
# # }
# #
# # normGR <- function(gr, minCTSS=NULL, minTPM=NULL){
# # 	# Filter on CTSS
# # 	if(!is.null(minCTSS)){
# # 		gr <- subset(gr, score >= minCTSS)
# # 	}
# #
# # 	# Filter on TPM
# # 	if(!is.null(minTPM)){
# # 		gr <- subset(gr, TPM(score) >= minTPM)
# # 	}
# #
# # 	# TPM normalize
# # 	score(gr) <- TPM(score(gr))
# #
# # 	# Return
# # 	gr
# # }
# #
# # # Low Mem loading
# # grl <- mclapply(bedBoth, function(x) DT2GR(fastReadBED(x)), mc.cores=3)
# # grl <- GRangesList(grl)
# #
# # # Normalize
# # grl1 <- endoapply(grl, normGR, minCTSS=2, minTPM=0.1)
# #
# # # Coverage
# # grs <- splitByStrand(unlist(grl1))
# #
# # # Plus cov
# # plus <- coverage(grs$`+`, weight="score")
# # minus <- coverage(grs$`-`, weight="score")
# #
# # # Final coverage
# # c(GRanges(plus, strand="+"), GRanges(minus, strand="-"))
# #
# # tpm1 <- lapply(grl1, normGR)
# # tpm2 <- endoapply(grl2, normGR)
# #
# #
# # grl2 <- mclapply(bedBoth, function(x) DT2GR(fastReadBED(x), rleScore=FALSE), mc.cores=3)
# # pryr::object_size(grl1)
# # pryr::object_size(grl2)
# #
# #
# # grl <- mergeReadCTSS(fnames=bedBoth, cores=2)
# # microbenchmark(splitReadCTSS(fnames=bedBoth), times=3)
# #
# #
# # ### GRanges code
# # dtList <- mclapply(bedBoth, function(x) DT2GR(fastReadBED(x)), mc.cores=3)
# # dtList <- GRangesList(dtList)
# # names(dtList) <- bedBoth
# # covgr <- CAGEfightR::coverageOfCTSS(dtList)
# #
# # # ### Code
# # #
# # # # Read list of DTs
# # # dtList <- mclapply(bedBoth, fastReadBED, mc.cores=3)
# # # names(dtList) <- bedBoth
# # #
# # # # GRlist
# # #
# # #
# # # Coverage
# # covgr <- DT2GR(collapseDT(dtList))
# #
# # #
# # # # TCs
# # # tcs <- CAGEfightR::simpleTagClustering(covgr)
# # #
# # # # Quantify
# # # CM <- mclapply(dtList, quantifyDT, features=tcs, mc.cores=3)
# # # CM <- do.call(cbind, CM)
# # #
# # # # Split GR
# # #
# # # splitGR <- function(x){
# # #
# # # }
# # #
# # # splitGRL <- function(x){
# # #
# # # }
# # #
# # # splitDT
