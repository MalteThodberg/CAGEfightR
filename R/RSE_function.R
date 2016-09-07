#' Calculate TPM for a CTSS file
#'
#' Adds a column "tpm" to a CTSS files containing (naive) Tag-per-million score for each bp position
#'
#' @param gr GRanges: Single CTSS file.
#'
#' @return GRanges with column "tpm" added.
#'
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import GenomicRanges
#' @export
calcTPM <- function(gr){
	# Calculate total sum in millions
	s <- sum(score(gr)) / 1e6

	# Calculate tpm
	tpm <- score(gr) / s

	# Append to gr and return
	gr$tpm <- tpm

	gr
}

#' Import CTSS file
#'
#' Uses the `fread`-function from the data.table package to quickly import a CTSS file from disk into a GRanges object.
#'
#' @param fname character: Path to ctss file on disk.
#'
#' @return GRanges representation of CTSS file
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
readCTSS <- function(fname){
	# IDEAS:
	# Default to rtracklayer::import if file is zipped
	# Perfrom post checks of col classes

	# Quickly read data using fread
	df <- data.table::fread(input=fname,
													data.table=FALSE,
													col.names=c("seqname", "start", "end", "name", "score", "strand"))

	# Make into GR and return
	GenomicRanges::makeGRangesFromDataFrame(df=df,
																					keep.extra.columns=TRUE,
																					starts.in.df.are.0based=TRUE)
}

#' Count overlaps for a CTSS file
#'
#' Counts overlaps of a CTSS-GRanges with a set of GRanges.
#'
#' @param tcs GRanges: Genomic Ranges to count overlaps
#' @param ctss GRanges: CTSS-GRanges
#'
#' @return integer vector of same length as tcs: Number of tags in each genomic range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
overlapCTSS <- function(tcs, ctss){
	# TODO:
	# Might not be exported...

	stopifnot(class(ctss) == "GRanges")

	# Find overlaps
	fo <- findOverlaps(query=tcs, subject=ctss)

	# Dataframe to be merged
	df <- data.frame(tcIndex=queryHits(fo),
									 ctssIndex=subjectHits(fo))
	df$ctssScores <- score(ctss)[df$ctssIndex]
	df <- stats::aggregate(ctssScores ~ tcIndex, data=df, FUN=sum)

	# Allocate zero count vector
	tagCounts <- rep(0, length(tcs))
	tagCounts[df$tcIndex] <- df$ctssScores

	# Return
	tagCounts
}

#' Quantify CTSS file into an Expression Matrix (EM)
#'
#' For a set of CTSS-GRanges and a set of Tag Cluster (TCs), count number of CAGE-tags in each TC in each sample and return an expression matrix.
#'
#' @param tcs GRanges: TC coordinates.
#' @param ctss GRangesList: Set of CTSS-GRanges.
#'
#' @return Integer matrix of dimensions: Number of TCs x Number of samples.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
quantifyFeatures <- function(tcs, ctss){
	# Option for None GRangesList...

	# Pre-checks
	stopifnot(class(ctss) == "GRangesList")

	# Build EM
	EM <- do.call(cbind,lapply(ctss, overlapCTSS, tcs=tcs))

	# Return
	EM
}


#' Simple Tag-Clustering
#'
#' Generates simple Tag Clusters (TCs) from CTSS-GRanges: Calculates TPM-coverage across all samples and finds cluster above a certain cutoff, merging cluster within a specified distance.
#'
#' @param ctss GRangesList: CTSS-GRanges.
#' @param tpmCutoff numeric: BP-positions with this or less TPM will not be used for tag clustering.
#' @param mergeDist integer: TCs within this distance of eachother are merged.
#'
#' @return GRanges of TCs, with TC peak position and summed TPM expression.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
findTagClusters <- function(ctss, tpmCutoff=1, mergeDist=25){
	# Better name to describe clustering

	### Pre-checks
	stopifnot(class(ctss)=="GRangesList")
	# Check if scores are all integers
	# Check if ctss files has seqinfo objects

	### Prepare ranges
	message("Preparing CTSS-GRanges")

	# Add TPM column
	ctss_tpm <- endoapply(ctss, calcTPM)

	# Split by strand
	ctss_plus <- endoapply(ctss_tpm, function(x) subset(x, strand=="+"))
	ctss_minus <- endoapply(ctss_tpm, function(x) subset(x, strand=="-"))

	### Coverage and peaks
	message("Calculating coverage and finding peaks...")

	# Calculate Strand-wise coverage
	coverage_plus <- coverage(ctss_plus, weight="tpm")
	coverage_minus <- coverage(ctss_minus, weight="tpm")

	# Peak calling: Find peaks
	peaks_plus <- slice(coverage_plus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
	peaks_minus <- slice(coverage_minus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)

	# Peak calling: Merge nearby peaks
	reduced_plus <- reduce(peaks_plus, min.gapwidth=mergeDist)
	reduced_minus <- reduce(peaks_minus, min.gapwidth=mergeDist)

	### Coverage and peaks
	message("Calculating global peak characteristics...")

	# Peaks stats: summed TPM
	views_plus <- Views(coverage_plus, reduced_plus)
	views_minus <- Views(coverage_minus, reduced_minus)

	sum_plus <- unlist(viewSums(views_plus))
	sum_minus <- unlist(viewSums(views_minus))

	# Peaks stats: Dominant TSS
	peaks_plus <- unlist(viewWhichMaxs(views_plus))
	peaks_minus <- unlist(viewWhichMaxs(views_minus))

	ranges_plus <- IRanges(start=peaks_plus, width=1)
	ranges_minus <- IRanges(start=peaks_minus, width=1)

	# # Peaks shapes: IQR
	# vectors_plus <- unlist(viewApply(views_plus, as.vector))
	# dist_plus <- lapply(vectors_plus, function(x) rep(seq_len(length(x)), x))
	# iqr_plus <- sapply(dist_plus, IQR)
	# skewness_plus <- sapply(dist_plus, moments::skewness)
	# kurtosis_plus <- sapply(dist_plus, moments::kurtosis)
	#
	# vectors_minus <- unlist(viewApply(views_minus, as.vector))
	# dist_minus <- lapply(vectors_minus, function(x) rep(seq_len(length(x)), x))
	# iqr_minus <- sapply(dist_minus, IQR)
	# skewness_minus <- sapply(dist_minus, moments::skewness)
	# kurtosis_minus <- sapply(dist_minus, moments::kurtosis)

	### Merge into final TC set
	message("Merging Tag Cluster info...")

	# Assemble by strand
	TCs_plus <- GRanges(reduced_plus, strand="+", peak=ranges_plus, score=sum_plus)#,
											#iqr=iqr_plus, skewness=skewness_plus, kurtosis=kurtosis_plus)
	TCs_minus <- GRanges(reduced_minus, strand="-", peak=ranges_minus, score=sum_minus)#,
											 #iqr=iqr_minus, skewness=skewness_minus, kurtosis=kurtosis_minus)

	# Merge into final objects
	TCs_both <- c(TCs_plus, TCs_minus)

	# Add ids
	TC_ids <- paste0(seqnames(TCs_both), ":", start(TCs_both), "-", end(TCs_both), ";", strand(TCs_both))

	names(TCs_both) <- TC_ids

	# Return
	TCs_both
}

#' Collect data into SummarizedExperiment
#'
#' Collects CTSS-GRanges, Tag Cluster (TCs), Expression Matrix (EM) and study design into a SummarizedExperiment. Additional options for documenting the object are available.
#'
#' @param ctss GRangesList: CTSS-GRanges.
#' @param tcs GRanges: TCs
#' @param em matrix, integer: EM
#' @param design data.frame: Study note
#'
#' @note This function does not sort or match the input, so arguments must be provided with matching columns, rows, etc.
#'
#' @return RangedSummarizedExperiment
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
assembleRSE <- function(ctss, tcs, em, design){
	### TO DO
	# Option for keepnig ctss in final object
	# Get column names from design rownames
	# Option for time stamping object
	# Option for author stamp

	# Preserve CTSS files
	sample_data <- S4Vectors::DataFrame(ctss=ctss, design)
	rownames(sample_data) <- rownames(design)
	colnames(em) <- rownames(design)

	# Assemble
	RSE <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=em),
															rowRanges=tcs,
															colData=sample_data)

	# Return
	RSE
}


#' Summarize CAGE experiment.
#'
#' This function is a user-friendly wrapper multiple other functions. Given a set of CTSS-GRanges, will find Tag Cluster (TCs) and quantify expression within them. Results are returned as a RangedSummarizedExperiment (RSE).
#'
#' @param ctss GRangesList: CTSS-GRanges.
#' @param design data.frame: Study design.
#' @param fun function: Function for finding TCs (See Details for additional information on customizing.).
#' @param ... additional arguments passed to fun
#'
#' @details TBA: Description of custom clustering
#' @note This function does not sort or match the input, so arguments must be provided with matching columns, rows, etc.
#' @return RangedSummarizedExperiment.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
summarizeCAGE <- function(ctss, design, fun=findTagClusters, ...){
	# Complete series
	TCs <- findTagClusters(ctss=ctss, ...)
	EM <- quantifyFeatures(tcs=TCs, ctss=ctss)
	SE <- assembleRSE(ctss=ctss, tcs=TCs, em=EM, design=design)

	# Return
	SE
}


# # Code heap
# library(IRanges)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(SummarizedExperiment)
#
# # Pombe files
# load("ctssFiles.rda")
#
# # Dummy design
# design <- tibble::tibble(tmp=names(ctss_grs))
# design <- as.data.frame(tidyr::separate(design, col=tmp,
# 													into=c("Genome", "Assembly", "Lane", "Medium", "Treatment", "Replicate"), sep="_"))
# rownames(design) <- apply(design, 1, function(x) paste(x[4:6], collapse=""))
#
#
# #TCs <- findTagClusters(ctssFiles=ctss_grs, tpmCutoff=0)
# #EM <- quantifyFeatures(tcs=TCs, ctss=ctss_grs)
# #SE <- assembleRSE(ctss=ctss_grs, tcs=TCs, em=EM, design=design)
# SE <- summarizeCAGE(ctss=ctss_grs, design=design)
#
# # Make into RSE
