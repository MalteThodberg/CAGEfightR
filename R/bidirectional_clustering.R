### Helpers

#' @import S4Vectors IRanges GenomicRanges
shiftedSummedCoverage <- function(gr, shift, width, weight="score", correct=TRUE){
	# Shifted coverage
	shifted_cov <- coverage(gr, weight=weight, shift=shift)

	# Running sum
	run_sum <- runsum(shifted_cov, k=width, endrule="constant")

	if(correct){
		run_sum[run_sum < 0] <- 0
	}

	# Return
	run_sum
}

#' Bhattacharyya coefficient.
#'
#' Quantifies departure from the reference bidirectional site (50%/50% of expression in bidirectional arms) using the Bhattacharyya coefficient as a measure of statistical distance.
#'
#' @param plusUpstream RleList or other: Sum of plus strand upstream (reference=0.0)
#' @param plusDownstream RleList or other: Sum of plus strand downstream (reference=0.5)
#' @param minusUpstream RleList or other: Sum of minus strand upstream (reference=0.0)
#' @param minusDownstream RleList or other: Sum of minus strand downstream (reference=0.5)
#'
#' @return Bhattacharyya coefficient
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors
#' @export
BC <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
	# Check all input have the same class
	input_classes <- c(class(plusUpstream), class(plusDownstream), class(minusUpstream), class(minusDownstream))
	single_class <- unique(input_classes)
	stopifnot(length(single_class) == 1)

	# Sum of all
	S <- plusUpstream + plusDownstream + minusUpstream + minusDownstream

	# Only calcuate downstream arms - rest will be zero
	B <- sqrt((minusDownstream / S) * 0.5) + sqrt((plusDownstream / S) * 0.5)

	# Checks
	stopifnot(class(B) == single_class)

	# Return
	B
}

### Unused Balance functions

# # Simple balance
# simpleBalance <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
# 	abs(plusDownstream - minusDownstream) / (plusDownstream + minusDownstream)
#
# 	# # Balance
# 	# E <- plusDownstream + minusDownstream
# 	# B <- abs(plusDownstream - minusDownstream) / E
# 	#
# 	# # Return
# 	# B
# }
#
# # Bidirectionality function
# penalizedBalance <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
# 	# Penalize
# 	penalizedPlusDownstream <- plusDownstream - minusUpstream
# 	penalizedMinusDownstream <- minusDownstream - plusUpstream
#
# 	# Balance
# 	E <- penalizedPlusDownstream + penalizedMinusDownstream
# 	B <- abs(penalizedPlusDownstream - penalizedMinusDownstream) / E
#
# 	# Return
# 	B
# }
#
# # Bidirectionality function
# penalizedBalance <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
# 	B1 <- abs(plusDownstream - minusDownstream) / (plusDownstream + minusDownstream)
# 	B2 <- plusUpstream + minusUpstream / (plusUpstream+plusDownstream+minusUpstream+minusDownstream)
# 	B2/B1
# }
#
# oddsRatio <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
# 	(plusDownstream * minusDownstream) / (plusUpstream * minusUpstream)
# }
#
# oddsRatio <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
# 	(plusDownstream * minusDownstream) / (plusUpstream * minusUpstream)
# }
#
# HellingerDistance <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
# 	B <- sqrt(1-BC(plusUpstream=plusUpstream,
# 							plusDownstream=plusDownstream,
# 							minusUpstream=minusUpstream,
# 							minusDownstream=minusDownstream))
# 	names(B) <- names(plusUpstream)
# 	B
# }

### Slice merge

#' Bidirectional clustering of genome-wide coverage
#'
#' Finds sites with (balanced and divergent) bidirectional transcription using sliding windows of summed coverage. Custom functions for calculating balance can be used.
#'
#' @param ctssCoverage GRanges: CTSS coverage.
#' @param window integer: Width of sliding window used for calculating sumes.
#' @param balanceThreshold numeric: Minimum balance score to be considered a bidirectional site.
#' @param balanceFun function: Function used to calculate balance (must have arguments plusUpstream, plusDownstream, minusUpstream, minusDownstream).
#' @param advancedStats logical: Whether to return internal statistics for defining bidirectional sites.
#' @param ... additional arguments passed to balanceFun..
#'
#' @return GRanges with bidirectional sites: Minimum width is 1 + 2*window, TPM sum (on both strands) in the score column, and the maximal bidirectional site in the thick column. If advancedStats=TRUE, additional columns are added to the GRanges.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
bidirectionalClustering <- function(ctssCoverage, window=199, balanceThreshold=0.99, balanceFun=BC, advancedStats=FALSE, ...){
	# Checks
	stopifnot(window %% 2 == 1, window >= 3)
	stopifnot(is.function(BC))

	# Obtain shift
	shift_val <- ceiling(window / 2)

	# Split into strand
	covByStrand <- splitByStrand(ctssCoverage)

	message("Summed coverage on plus strand...")
	PD <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=-shift_val, width=window)
	PU <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=shift_val, width=window)
	#P <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=0, width=window)


	message("Summed coverage on minus strand...")
	MD <- shiftedSummedCoverage(gr=covByStrand$`-`, shift=shift_val, width=window)
	MU <- shiftedSummedCoverage(gr=covByStrand$`-`, shift=-shift_val, width=window)
	rm(covByStrand)

	message("Calculating balance score...")
	# Calculate balance
	B <- balanceFun(plusUpstream=PU, plusDownstream=PD, minusUpstream=MU, minusDownstream=MD, ...)
	rm(PD, PU, MD, MU)
	#BAL <<- GRanges(B)

	# Set NAs to 0 to allow slicing
	if(min(min(B, na.rm=TRUE)) < 0){
		stop("Balance function produced values below 0! Output values must be in range [0-Inf) to allow for slicing...")
	}else{
		B[is.na(B)] <- 0
	}
	#B[is.na(B)] <- -1
	#B[is.na(B)] <- 0

	message("Slice-merge to find bidirectional clusters...")
	# Slice
	#bidirLoci <- slice(B, lower=0, upper=balanceThreshold, rangesOnly=TRUE)
	bidirLoci <- slice(B, lower=balanceThreshold, upper=Inf, rangesOnly=TRUE)

	# Merge
	mergeDist <- window * 2
	bidirLoci <- reduce(bidirLoci, min.gapwidth=mergeDist)
	stopifnot(isDisjoint(bidirLoci))

	# Expand by window size
	#expLoci <- bidirLoci
	start(bidirLoci) <- methods::as(start(bidirLoci) - window, "IntegerList")
	end(bidirLoci) <- methods::as(end(bidirLoci) + window, "IntegerList")
	stopifnot(isDisjoint(bidirLoci))

	message("Calculating statistics...")
	# Coverage on both strands
	A <- coverage(ctssCoverage, weight="score")

	# Views on score and balance
	viewsA <- Views(A, bidirLoci)
	viewsB <- Views(B, bidirLoci)
	rm(A, B)

	#  Scores and midpoint
	scores <- unlist(viewSums(viewsA))
	#balance <- viewMins(viewsB)
	#balance <- viewMaxs(viewsB)
	peaks <-  viewRangeMaxs(viewsB)
	peaks <- resize(x=peaks, width=1, fix="center", use.names=FALSE)

	# Assemble output
	o <- GRanges(bidirLoci, score=scores, thick=unlist(peaks))
	rm(scores, peaks, viewsA)

	# Extra stats if needed?
	if(advancedStats){
		message("Calculating additional statistics...")
		o$balance <- unlist(viewMaxs(viewsB)) # Only balance for now
	}
	rm(viewsB)

	# Carry over seqinfo
	message("Preparing output...")
	seqinfo(o) <- seqinfo(ctssCoverage)
	o <- sort(o)

	# Names as IDs for both ranges and peaks
	o_ids <- paste0(seqnames(o), ":", start(o), "-", end(o))
	names(o) <- o_ids
	names(o$thick) <- o_ids
	rm(o_ids)

	# Return
	o
}

#' Calculate sample-wise bidirectional transcription
#'
#' For each bidirectional site, calculates how many individual samples shows transcription in both directions.
#'
#' @param ctss GRangesList or SimpleList: GRanges with CTSSs in the score column.
#' @param features GRanges: Bidirectional sites (with the maximal bidirectional site in the thick column).
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#'
#' @return Integer-vector of same length as features with the number of samples showing bidirectional transcription.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
calculateBidirectionality <- function(ctss, features, biocParallel=bpparam()){
	stopifnot(all(strand(features) == "*"))

	# Extract arms of enhances
	message("Extracting stranded arms...")
	arms_plus <- features
	start(arms_plus) <- end(arms_plus$thick)
	strand(arms_plus) <- "+"

	arms_minus <- features
	end(arms_minus) <- start(arms_minus$thick)
	strand(arms_minus) <- "-"

	# Quantify
	message(sprintf("Assessing bidirectionality of %d features in %d samples...",
									length(features), length(ctss)))
	mat_plus <- bplapply(ctss, overlapsAny, query=arms_plus, BPPARAM=biocParallel)
	mat_plus <- do.call(cbind, mat_plus)
	rm(arms_plus)

	mat_minus <- bplapply(ctss, overlapsAny, query=arms_minus, BPPARAM=biocParallel)
	mat_minus <- do.call(cbind, mat_minus)
	rm(arms_minus)

	# Compare and counts
	message("Counting sample-wise bidirectionality...")
	o <- mat_plus & mat_minus
	rm(mat_plus, mat_minus)
	o <- rowSums(o)

	# Checks
	stopifnot(length(o) == length(features))
	stopifnot(max(o) <= length(ctss))

	# Return
	o
}

### Enhancer tuning - NOT WORKING YET
#
# countBidir <- function(B, threshold, window){
# 	# Slice
# 	bidirLoci <- slice(B, lower=threshold, upper=Inf, rangesOnly=TRUE)
#
# 	# Merge
# 	mergeDist <- window * 2
# 	bidirLoci <- reduce(bidirLoci, min.gapwidth=mergeDist)
#
# 	# Count
# 	sum(elementNROWS(bidirLoci))
# }
#
# tuneBidirectionalClustering <- function(ctssCoverage, window=199, stepSize=10, balanceFun=BC, ...){
# 	# Checks
# 	stopifnot(window %% 2 == 1)
#
# 	# Obtain shift
# 	shift_val <- ceiling(window / 2)
#
# 	# Split into strand
# 	covByStrand <- splitByStrand(ctssCoverage)
#
# 	message("Summed Coverage on plus strand...")
# 	PD <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=-shift_val, width=window)
# 	PU <- shiftedSummedCoverage(gr=covByStrand$`+`, shift=shift_val, width=window)
#
# 	message("Summed Coverage on minus strand...")
# 	MD <- shiftedSummedCoverage(gr=covByStrand$`-`, shift=shift_val, width=window)
# 	MU <- shiftedSummedCoverage(gr=covByStrand$`-`, shift=-shift_val, width=window)
# 	rm(covByStrand)
#
# 	message("Calculating balance score...")
# 	# Calculate balance
# 	B <- balanceFun(plusUpstream=PU, plusDownstream=PD, minusUpstream=MU, minusDownstream=MD, ...)
# 	rm(PD, PU, MD, MU)
#
# 	# Set NAs to 0 to allow slicing
# 	if(min(min(B, na.rm=TRUE)) < 0){
# 		stop("Balance function produced values below 0! Output values must be in range [0-Inf) to allow for slicing...")
# 	}else{
# 		B[is.na(B)] <- 0
# 	}
#
# 	message("Slice-merge to find bidirectional clusters...")
# 	o <- data.frame(thresholds=seq(0.75, 1, 0.01))
# 	o$n <- unlist(bplapply(o$thresholds, countBidir, B=B, window=200))
#
# 	bidirLoci <- slice(B, lower=0.75, upper=Inf, rangesOnly=TRUE)
#
# 	# Merge
# 	mergeDist <- window * 2
# 	bidirLoci <- reduce(bidirLoci, min.gapwidth=mergeDist)
# 	stopifnot(isDisjoint(bidirLoci))
#
# 	}

### Snippet Code
#
# #bal <- bidirectionalClustering(ctssCoverage, balanceFun=penalizedBalance, balanceThreshold=0.01)
# #pen <- bidirectionalClustering(ctssCoverage, balanceFun=penalizedBalance, balanceThreshold=0.01)
# #or <- bidirectionalClustering(ctssCoverage, balanceFun=oddsRatio, balanceThreshold=40)
# #bc <- bidirectionalClustering(ctssCoverage, balanceFun=HellingerDistance, balanceThreshold=0.1)
# bc <- bidirectionalClustering(ctssCoverage, balanceThreshold=0.95, advancedStats =TRUE)
# bc$bidirectionality <- calcBidirectionality(ctss=ctssList, features=bc)
# bc$txType <- assignTxType(bc, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#
# library(Gviz)
# options(ucscChromosomeNames=FALSE)
# axisTrack <- GenomeAxisTrack(littleTicks=TRUE, exponent=0)
# geneTrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
# covTrack <- CAGEfightR::strandedCoverageTrack(gr=ctssCoverage, name="TPM")
#
# enhancers <- subset(bc, bidirectionality >= 1 & txType %in% c("intergenic", "intron"))
# s <- enhancers[order(score(enhancers), decreasing=TRUE)[1]]
# plotTracks(trackList=list(axisTrack, geneTrack, covTrack),
# 					 from=start(s)-1000, to=end(s)+1000, chromosome=seqnames(s))
#
#
#
# library(Gviz)
# options(ucscChromosomeNames=FALSE)
# axisTrack <- GenomeAxisTrack(littleTicks=TRUE, exponent=0)
# covTrack <- CAGEfightR::strandedCoverageTrack(gr=ctssCoverage, name="TPM")
# balanceTrack <- DataTrack(shift(BAL, 0), name="Balance", type="l", col="forestgreen")
# # PDTrack <- DataTrack(shift(GRanges(PD), 0), name="PlusDownstream", type="l")
# # MDTrack <- DataTrack(shift(GRanges(MD), 0), name="MinusDownstream", type="l")
# # PUTrack <- DataTrack(shift(GRanges(PU), 0), name="PlusUpstream", type="l")
# # MUTrack <- DataTrack(shift(GRanges(MU), 0), name="MinusUpstream", type="l")
# # PTrack <- DataTrack(shift(GRanges(P), 0), name="Plus", type="l")
#
# downstream <- DataTrack(bindAsGRanges(plus=PD, minus=MD), name="Downstream", type="l", groups=c("plus", "minus"), col=c("tomato", "cornflowerblue"))
# upstream <- DataTrack(bindAsGRanges(plus=PU, minus=MU), name="Upstream", type="l", groups=c("plus", "minus"), col=c("tomato", "cornflowerblue"))
#
# bc$txType <- assignTxType(bc, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
#
#
# # Plot with peak
# #s <- bc[5]
# s <- bc[order(score(bc), decreasing=TRUE)[1]]
#
# plotTracks(trackList=list(axisTrack, covTrack, downstream, upstream, balanceTrack),
# 					 from=start(s), to=end(s), chromosome=seqnames(s))
#
# plotTracks(trackList=list(axisTrack, covTrack, PDTrack, PTrack, MDTrack, PUTrack, MUTrack, balanceTrack),
# 					 from=start(s), to=end(s), chromosome=seqnames(s))

# ### Deep enhancers - ORIGINAL CODE
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
# ### Bidirectionality functions
# # Simple balance
# simpleBalance <- function(FU, RD, RU, FD){
# 	E <- FU + RD
# 	B <- abs(FU - RD) / E
# 	B
# }
#
# # Bidirectionality function
# simpleBalance <- function(FU, RD, RU, FD){
# 	E <- FU + RD
# 	B <- abs(FU - RD) / E
# 	B
# }
#
# bidirectionalClustering <- function(ctssCoverage, window, balanceThreshold=0.4, mergeDist=20, balanceFun=simpleBalance, ...){
# 	# Checks
# 	stopifnot(window %% 2 == 1)
#
# 	# Obtain shift
# 	shift_val <- ceiling(window / 2)
#
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
# 	#rm(covByStrand)
#
# 	message("Calculating balance score")
# 	# Calculate balance
# 	#B <- balanceFun(FU=FU, RD=RD, FD=FD, RU=RU, ...)
# 	B <- balanceFun(FU=FU, RD=RD, FD=FD, RU=RU)
# 	#rm(FU, RD, FD, RU)
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
# 	#rm(rl, A, B)
#
# 	# Midpoints and scores
# 	score(bidirLoci) <- unlist(viewSums(viewsA))
# 	bidirLoci$peak <- IRanges(unlist(viewWhichMaxs(viewsB)), width=1)
# 	#rm(viewsA, viewsB)
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
