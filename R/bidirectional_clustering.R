# ### Deep enhancers
#
# library(data.table)
#
# # Set seqinfo to make strand matchable
# setDummySeqInfo <- function(gr){
# 	# Coerce to data.frame
# 	dt <- as.data.table(gr)
#
# 	# Maximum length
# 	max_length <- dt[, max(end), by=seqnames]
#
# 	# Match to seqinfo
# 	tmp <- match(max_length$seqnames, table=seqlevels(gr))
#
# 	# Modify seqinfo
# 	seqlengths(gr) <- max_length$V1[tmp]
#
# 	# Return object
# 	gr
# }
#
# lowMem <- setDummySeqInfo(lowMem)
#
# # Split into strand
# plus_strand <- subset(lowMem, strand=="+")
# minus_strand <- subset(lowMem, strand=="-")
#

# # Function for shifted summed coverage
# shiftedSummedCoverage <- function(gr, shift, k, weight="score", correct=TRUE){
# 	# Shifted coverage
# 	shifted_cov <- coverage(gr, weight=weight, shift=shift)
#
# 	# Running sum
# 	run_sum <- runsum(shifted_cov, k=k)
#
# 	if(correct){
# 		run_sum[run_sum < 0] <- 0
# 	}
#
# 	# Return
# 	run_sum
# }
#
# # Bidirectionality function
# simpleBalance <- function(FU, RD, RU, FD){
# 	E <- FU + RD
# 	B <- abs(FU - RD) / E
# 	B
# }
#
#
# bidirectionalLoci <- function(ctssCoverage, window, balanceThreshold, mergeDist, balanceFun=simpleBalance, ...){
# 	# Rescale window
# 	shift_val <- window / 2
# 	k_val <- window - 1
# s
# 	# Split into strand
# 	covByStrand <- splitByStrand(ctssCoverage)
#
# 	message("Summed Coverage on plus strand...")
# 	FU <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=shift_val, k=k_val)
# 	FD <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=-shift_val, k=k_val)
#
# 	message("Summed Coverage on minus strand...")
# 	RU <- shiftedSummedCoverage(gr=covByStrand$`-`, shift=shift_val, k=k_val)
# 	RD <- shiftedSummedCoverage(gr=covByStrand$`-`, shift=-shift_val, k=k_val)
# 	rm(covByStrand)
#
# 	message("Calculating balance score")
# 	# Calculate balance
# 	B <- balanceFun(FU=FU, RD=RD, FD=FD, RU=RU, ...)
# 	rm(FU, RD, FD, RU)
#
# 	# Set NAs to -1 to allow slicing
# 	B[is.na(B)] <- -1
#
# 	message("Slice-merge to find loci")
# 	# Slice
# 	bidirLoci <- GRanges(slice(B, lower=0, upper=balanceThreshold, rangesOnly=TRUE))
#
# 	# Merge
# 	bidirLoci <- reduce(bidirLoci, min.gapwidth=mergeDist)
#
# 	# Expand by window size
# 	start(bidirLoci) <- start(bidirLoci) - window
# 	end(bidirLoci) <- end(bidirLoci) + window
#
# 	message("Calculating statistics")
# 	# Coverage on both strands
# 	A <- coverage(ctssCoverage, weight="score")
#
# 	# Views on score and balance
# 	rl <- as(bidirLoci, "RangesList")
# 	viewsA <- Views(A, rl)
# 	viewsB <- Views(B, rl)
# 	rm(rl, A, B)
#
# 	# Midpoints and scores
# 	score(bidirLoci) <- unlist(viewSums(viewsA))
# 	bidirLoci$peak <- IRanges(unlist(viewWhichMaxs(viewsB)), width=1)
# 	rm(viewsA, viewsB)
#
# 	# Return
# 	bidirLoci
# }

#
# enhancers <- bidirectionalLoci(ctssCoverage=lowMem,
# 															window=200,
# 															balanceThreshold=0.25,
# 															mergeDist=10,
# 															balanceFun=simpleBalance)
#
# # Investigate enhancer
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#
# enhancers$txType <- assignTxType(gr=enhancers, txdb=txdb)
#
# table(enhancers$txType)
#
# # Enhancer windows search size
# w <- 200
# FS <- shiftedSummedCoverage(gr=plus_strand, weight="score", shift=w/2, k=w+1)
# FA <- shiftedSummedCoverage(gr=plus_strand, weight="score", shift=-w/2, k=w+1)
# RA <- shiftedSummedCoverage(gr=minus_strand, weight="score", shift=w/2, k=w+1)
# RS <- shiftedSummedCoverage(gr=minus_strand, weight="score", shift=-w/2, k=w+1)
#
# # balance
# expr <- FS + RS
# balance <- abs(FS - RS) / expr
#
# # Replace NAs with negative to allow slicing
# balance[is.na(balance)] <- -1
#
# # # Split into coverage
# # plus_cov_right <- coverage(plus_strand, weight="score", shift=100)
# # plus_cov_left <- coverage(plus_strand, weight="score", shift=-100)
# # minus_cov_right <- coverage(minus_strand, weight="score", shift=100)
# # minus_cov_left <- coverage(minus_strand, weight="score", shift=-100)
# #
# # # Running sum
# # plus_run_right <- runsum(plus_cov_right, k=201)
# # plus_run_left <- runsum(plus_cov_left, k=201)
# # minus_run_right <- runsum(minus_cov_right, k=201)
# # minus_run_left <- runsum(minus_cov_left, k=201)
# #
# # # Correct numerical instability
# # plus_run_right[plus_run_right < 0] <- 0
# # plus_run_left[plus_run_left < 0] <- 0
# # minus_run_right[minus_run_right < 0] <- 0
# # minus_run_left[minus_run_left < 0] <- 0
# #
# # # Balance score
# # expr <- plus_run_right + minus_run_left
# # balance <- abs(plus_run_right - minus_run_left) / expr
# #
# # # Replace NAs with negative to allow slicing
# # balance[is.na(balance)] <- -1
#
# # Slice-merge
# proto_enhancers <- GRanges(slice(balance, lower=0, upper=0.4, rangesOnly=TRUE))
# proto_enhancers <- reduce(proto_enhancers, min.gapwidth=10)
#
# # Expand
# enhancers <- proto_enhancers
# start(enhancers) <- start(enhancers) - w
# end(enhancers) <- end(enhancers) + w
#
# # Coverage on both strands
# both_cov <- coverage(lowMem, weight="score")
#
# # Correct numerical instability
# #both_cov[both_cov < 0] <- 0
#
# # Views
# cov_views <- Views(both_cov, as(enhancers, "RangesList"))
# bal_views <- Views(balance, as(enhancers, "RangesList"))
#
# # Midpoints and scores
# enhancer_scores <- unlist(viewSums(cov_views))
# enhancer_midpoints <- unlist(viewWhichMaxs(bal_views))
# enhancer_midpoints <- IRanges(start=enhancer_midpoints, width=1)
#
# # Merge output
# score(enhancers) <- enhancer_scores
# enhancers$peak <- enhancer_midpoints
#
# # Investigate enhancer
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#
# enhancers$txType <- assignTxType(gr=enhancers, txdb=txdb)
#
# table(enhancers$txType)
#
# subset(enhancers, txType %in% c("intergenic", "intron"))
#
# ggplot(as.data.frame(enhancers), aes(x=txType, y=log2(score))) +
# 	geom_boxplot()
#
# ggplot(as.data.frame(subset(enhancers, txType %in% c("intergenic", "intron"))),
# 			 aes(x=log10(width), y=log2(score), color=txType)) +
# 	geom_point(alpha=0.75)
#
# library(ggplot2)
