#' Internal function: countOverlaps with scores
#'
#' Similar to countOverlaps, but takes into account the scores column.
#'
#' @param gr GRanges: Ranges with scores.
#' @param features GRanges: Ranges to be quantified.
#'
#' @return numeric vector of same length as features.
#'
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Quantification functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
countScoredOverlaps <- function(gr, features){
	# Solution from bioconductor by M Morgan
	hits <- methods::as(findOverlaps(query=features, subject=gr), "List")
	weightedCount <- sum(extractList(score(gr), hits))

	# Return
	weightedCount
}

#' Quantify CTSS across features
#'
#' Counts the number of CTSS-tags in a set of features
#'
#' @param ctss GRangesList or SimpleList: GRanges with CTSSs in the score column.
#' @param features GRanges: Ranges to be quantified.
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#'
#' @return matrix with rows as features and columns as samples.
#'
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Quantification functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
quantifyFeatures <- function(ctss, features, biocParallel=bpparam()){
	# Get some basic info for printing
	message(sprintf("Quantifying expression of %d features in %d samples...",
									length(features), length(ctss)))

	# Quantify each GRange
	m <- bplapply(ctss, countScoredOverlaps, features=features, BPPARAM=biocParallel)

	# Merge into matrix
	message("Building Expression Matrix...")
	m <- do.call(cbind, m)

	# Attempt to set names
	rownames(m) <- names(features)
	colnames(m) <- names(ctss)

	# Return
	m
}
