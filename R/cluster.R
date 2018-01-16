### Main S4 functions

#' Tag Clustering of pooled CTSSs.
#'
#' Finds unidirectional Tag Cluster (TCs) with a pooled TPM above a certain threshold using a slice-reduce approach. Addtionally calculates the sum and peak position of the TCs.
#'
#' @param object GRanges or RangedSummarizedExperiment: Basepair-wise pooled CTSS.
#' @param pooledCutoff numeric: Minimum pooled value to be considered as TC.
#' @param mergeDist integer: Merge TCs within this distance.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with TPM sum as the score column, and TC peak as the thick column.
#'
#' @family Tag-clustering functions
#' @export
setGeneric("clusterUnidirectionally", function(object, ...) {
	standardGeneric("clusterUnidirectionally")
})

#' @import assertthat S4Vectors IRanges GenomicRanges
#' @rdname clusterUnidirectionally
setMethod("clusterUnidirectionally", signature(object="GenomicRanges"), function(object, pooledCutoff=0, mergeDist=20){
	# Pre-checks
	assert_that(isDisjoint(object),
							!is.null(score(object)),
							is.numeric(score(object)),
							not_empty(seqlengths(object)),
							is.number(pooledCutoff),
							is.count(mergeDist))

	# Split coverage by strand
	message("Splitting by strand...")
	coverage_stranded <- splitByStrand(object)
	coverage_plus <- coverage(coverage_stranded$`+`, weight="score")
	coverage_minus <- coverage(coverage_stranded$`-`, weight="score")
	rm(coverage_stranded)

	# Peak calling: Find peaks
	message("Finding tag clusters...")
	peaks_plus <- slice(coverage_plus, lower=pooledCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
	peaks_minus <- slice(coverage_minus, lower=pooledCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)

	# Peak calling: Merge nearby peaks
	reduced_plus <- reduce(peaks_plus, min.gapwidth=mergeDist)
	reduced_minus <- reduce(peaks_minus, min.gapwidth=mergeDist)

	### Coverage and peaks
	message("Calculating tag cluster statistics...")
	TCs <- TCstats(coverage_plus=coverage_plus, coverage_minus=coverage_minus,
								 tcs_plus=reduced_plus, tcs_minus=reduced_minus)

	# Carry over seqinfo and sort
	message("Preparing output...")
	seqinfo(TCs) <- seqinfo(object)
	TCs <- sort(TCs)

	# Print some basic stats
	message("Tag clustering summary:")
	summarizeWidths(TCs)

	# Return
	TCs
})

#' @import assertthat S4Vectors IRanges GenomicRanges
#' @rdname clusterUnidirectionally
setMethod("clusterUnidirectionally", signature(object="RangedSummarizedExperiment"), function(object, ...){
	clusterUnidirectionally(rowRanges(object), ...)
})

#' @import assertthat S4Vectors IRanges GenomicRanges
#' @rdname clusterUnidirectionally
setMethod("clusterUnidirectionally", signature(object="GPos"), function(object, ...){
	warning("Using temporary GPos-method in clusterUnidirectionally!")
	clusterUnidirectionally(methods::as(object, "GRanges"), ...)
})

### Helper functions

#' @import S4Vectors IRanges GenomicRanges
TCstats <- function(coverage_plus, coverage_minus, tcs_plus, tcs_minus){
	# Check classes
	stopifnot(class(coverage_plus)=="SimpleRleList",
						class(coverage_minus)=="SimpleRleList",
						class(tcs_plus)=="CompressedIRangesList",
						class(tcs_minus)=="CompressedIRangesList")

	# Check seqlevels
	stopifnot(length(unique(list(names(coverage_plus),
															 names(tcs_plus),
															 names(coverage_minus),
															 names(tcs_minus)))) == 1)

	# Obtain views
	views_plus <- Views(coverage_plus, tcs_plus)
	views_minus <- Views(coverage_minus, tcs_minus)

	# Calculate Sums
	sum_plus <- unlist(viewSums(views_plus))
	sum_minus <- unlist(viewSums(views_minus))

	# Find peaks
	ranges_plus <- viewRangeMaxs(views_plus)
	ranges_minus <- viewRangeMaxs(views_minus)
	ranges_plus <- resize(unlist(ranges_plus), width=1, fix="center")
	ranges_minus <- resize(unlist(ranges_minus), width=1, fix="center")

	# Merge into GRanges
	TCs <- c(GRanges(tcs_plus, strand="+", score=sum_plus, thick=ranges_plus),
					 GRanges(tcs_minus, strand="-", score=sum_minus, thick=ranges_minus))

	# Names as IDs for both ranges and peaks
	TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
	names(TCs) <- TC_ids
	names(TCs$thick) <- TC_ids

	# Return
	TCs
}

#' @import S4Vectors IRanges GenomicRanges
summarizeWidths <- function(gr){
	# Checks
	stopifnot(class(gr) == "GRanges")

	# Cut up widths
	x <- cut(width(gr),
					 breaks=c(1, 10, 100, 1000, Inf),
					 labels=c(">=1", ">=10", ">=100", ">=1000"),
					 include.lowest=TRUE)

	# Get freqs and props
	y <- table(Width=x)
	z <- prop.table(y)

	# Format to data.frame
	w <- merge(as.data.frame(y, responseName="Count"),
						 as.data.frame(z, responseName="Percent"))

	# Add Total row
	w <- rbind(data.frame(Width="Total",
												Count=sum(w$Count),
												Percent=sum(w$Percent)),
						 w)

	# Reformat to percent
	w$Percent <- paste0(format(w$Percent * 100, digits=1), " %")

	# To string and message
	s <- paste(utils::capture.output(print(w, row.names=FALSE)), collapse="\n")
	message(s)
}


