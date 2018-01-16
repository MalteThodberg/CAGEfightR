#### Main S4 functions ####

#' Determine the optimal pooled threshold for unidirectional tag clustering.
#'
#' This function counts the number of Tag Clusters (TCs) for an exponentially increasing series of expression cutoffs.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Pooled CTSS.
#' @param steps integer: Number of thresholds to analyze (in addition to treshold=0).
#' @param maxExponent numeric: The maximal threshold to analyse is obtained as min(score)*2^maxExponent.
#' @param mergeDist integer: Merge TCs within this distance.
#' @param searchMethod character: For advanced user only, see details.
#' @param ... additional arguments passed to methods.
#'
#' @return data.frame with two columns: threshold and nTCs (number of Tag Clusters)
#' @family Tag-clustering functions
#' @export
setGeneric("tuneTagClustering", function(object, ...) {
	standardGeneric("tuneTagClustering")
})

#' @import assertthat S4Vectors IRanges GenomicRanges
#' @rdname tuneTagClustering
setMethod("tuneTagClustering", signature(object="GenomicRanges"), function(object, steps=10L, maxExponent=1, mergeDist=20L, searchMethod="exponential"){
	# Pre-checks
	assert_that(isDisjoint(object),
							!is.null(score(object)),
							is.numeric(score(object)),
							not_empty(seqlengths(object)),
							is.count(steps),
							is.numeric(maxExponent),
							is.count(mergeDist))

	# Setup series
	message("Finding thresholds to be tested...")
	if(searchMethod == "exponential"){
		exp_series <- min(score(object)) * 2^seq(from=0, to=maxExponent, length.out=steps)
	}else if(searchMethod == "minUnique"){
		steps_series <- seq_len(steps)
		exp_series <- sort(unique(score(object)), partial=steps_series)[steps_series]
		rm(steps_series)
	}else{
		stop("searchMethod must be either exponential or minUnique")
	}

	# Add zero
	exp_series <- c(0, exp_series)

	# Split coverage by strand
	message("Splitting coverage by strand...")
	by_strand <- splitByStrand(object)
	coverage_plus <- coverage(by_strand$`+`, weight="score")
	coverage_minus <- coverage(by_strand$`-`, weight="score")
	rm(by_strand)

	# Print some info
	message("Iterating over ", length(exp_series), " thresholds using ", BiocParallel::bpworkers(), " worker(s)...")

	# Count
	message("Analyzing plus strand...")
	n_plus <- BiocParallel::bpvec(exp_series, countClusters, cv=coverage_plus, mergeDist=mergeDist)
	message("Analyzing minus strand...")
	n_minus <- BiocParallel::bpvec(exp_series, countClusters, cv=coverage_minus, mergeDist=mergeDist)

	# Assemble results
	message("Preparing output...")
	o <- data.frame(threshold=exp_series, nTCs=n_plus+n_minus)

	# Return
	o
})

#' @import SummarizedExperiment
#' @rdname tuneTagClustering
setMethod("tuneTagClustering", signature(object="RangedSummarizedExperiment"), function(object, ...){
	tuneTagClustering(rowRanges(object), ...)
})

#' @import SummarizedExperiment
#' @rdname tuneTagClustering
setMethod("tuneTagClustering", signature(object="GPos"), function(object, ...){
	warning("Using temporary GPos-method in tuneTagClustering!")
	tuneTagClustering(methods::as(object, "GRanges"), ...)
})

#### Helpers ####

#' @import S4Vectors IRanges GenomicRanges
countClusters <- function(thresholds, cv, mergeDist=20){
	# Slice
	o <- lapply(thresholds, slice, x=cv, includeLower=FALSE, upper=Inf, rangesOnly=TRUE)

	# Reduce
	o <- lapply(o, reduce, min.gapwidth=mergeDist)

	# Count
	o <- lapply(o, elementNROWS)

	# Sum
	o <- vapply(o, sum, integer(1))

	# Return
	o
}
