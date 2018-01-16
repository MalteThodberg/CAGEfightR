#### Main shape function ####

#' Quantify Tag Cluster shapes
#'
#' Apply a shape-function to the pooled CTSS of every TC.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Tag clusters.
#' @param pooled GenomicRanges or RangedSummarizedExperiment: Pooled CTSS as the score column.
#' @param outputColumn character: Name of column to hold shape statistics.
#' @param shapeFunction function: Function to apply to each TC (See details).
#' @param ... additional arguments passed to shapeFunction.
#'
#' @return Adds a column
#' @family Shape functions
#' @export
setGeneric("calcShape", function(object, pooled, ...) {
	standardGeneric("calcShape")
})

#' @import assertthat S4Vectors IRanges GenomicRanges
#' @rdname calcShape
setMethod("calcShape", signature(object="GRanges", pooled="GenomicRanges"), function(object, pooled, outputColumn="IQR", shapeFunction=shapeIQR, ...){
	# Pre-checks
	assert_that(!is.null(score(pooled)),
							is.numeric(score(pooled)),
							isDisjoint(pooled),
							is.character(outputColumn),
							is.function(shapeFunction),
							identical(seqlengths(object), seqlengths(pooled)))

	# Warnings
	if(outputColumn %in% colnames(mcols(object))){
		warning("object already has a column named ", outputColumn," in mcols: It will be overwritten!")
	}

	# Names need to be set for sorting later
	if(is.null(names(object))){
		message("Adding names...")
		names(object) <- paste("TC", seq_along(object))
	}

	# Split by strand
	message("Splitting coverage by strand...")
	covByStrand <- splitByStrand(pooled)
	TCsByStrand <- splitByStrand(object)

	# Coverage by strand
	coverage_plus <- coverage(covByStrand$`+`, weight="score")
	coverage_minus <- coverage(covByStrand$`-`, weight="score")
	rm(covByStrand)

	# Views
	message("Applying function to each TC...")
	views_plus <- Views(coverage_plus, methods::as(TCsByStrand$`+`, "RangesList"))
	views_minus <- Views(coverage_minus, methods::as(TCsByStrand$`-`, "RangesList"))
	rm(coverage_plus, coverage_minus)

	# Applying functions to views
	stat_plus <- viewApply(views_plus, shapeFunction, ...)
	stat_minus <- viewApply(views_minus, shapeFunction, ...)

	# Reset names
	message("Assembling output...")
	stat_plus <- as.numeric(unlist(stat_plus, use.names=FALSE))
	names(stat_plus) <- names(TCsByStrand$`+`)
	stat_minus <- as.numeric(unlist(stat_minus, use.names=FALSE))
	names(stat_minus) <- names(TCsByStrand$`-`)
	rm(TCsByStrand)

	# Reassemble in same order
	o <- c(stat_plus, stat_minus)
	o <- o[match(names(object), names(o))]
	names(o) <- NULL
	rm(stat_plus, stat_minus)

	# Add to object
	mcols(object)[,outputColumn] <- o

	# Return
	object
})

#' @import SummarizedExperiment
#' @rdname calcShape
setMethod("calcShape", signature(object="RangedSummarizedExperiment", pooled="GenomicRanges"), function(object, pooled, ...){
	rowRanges(object) <- calcShape(rowRanges(object), pooled, ...)
	object
})

#' @import SummarizedExperiment
#' @rdname calcShape
setMethod("calcShape", signature(object="GRanges", pooled="RangedSummarizedExperiment"), function(object, pooled, ...){
	calcShape(object, rowRanges(pooled), ...)
})

#' @import SummarizedExperiment
#' @rdname calcShape
setMethod("calcShape", signature(object="RangedSummarizedExperiment", pooled="RangedSummarizedExperiment"), function(object, pooled, ...){
	rowRanges(object) <- calcShape(rowRanges(object), rowRanges(pooled), ...)
	object
})

#' @import SummarizedExperiment
#' @rdname calcShape
setMethod("calcShape", signature(object="GRanges", pooled="GPos"), function(object, pooled, ...){
	warning("Using temporary GPos-method in calcShape!")
	calcShape(object, methods::as(pooled, "GRanges"), ...)

})

#### Individual shape functions ####

#' Shape statitic: Interquartile range
#'
#' Calculates the interquartile range of a vector.
#'
#' @param x numeric Rle vector: Coverage series.
#' @param lower numeric: Lower quartile.
#' @param upper numeric: Upper quartile.
#'
#' @return Numeric
#' @family Shape functions
#' @import S4Vectors
#' @export
shapeIQR <- function(x, lower=0.25, upper=0.75){
	# To normal vector
	x <- as.vector(x)

	# Scale by sum
	x <- x / sum(x)

	# Cumulate sum
	x <- cumsum(x)

	# Find the threshold
	lowerPos <- Position(function(y) y >= lower, x)
	upperPos <- Position(function(y) y >= upper, x)

	# Return difference
	upperPos - lowerPos
}

#' Shape statistic: Shannon Entropy
#'
#' Calculates the Shannon Entropy (base log2) for a vector. Zeros are removed before calculation.
#'
#' @param x numeric Rle vector: Coverage series.
#'
#' @return Numeric.
#' @family Shape functions
#' @import S4Vectors
#' @export
shapeEntropy <- function(x){
	# To normal vector
	#x <- as.vector(x[x > 0])
	x <- as.vector(x)

	# Scale by sum
	x <- x / sum(x)

	# Calculate entropy
	o <- suppressWarnings(-sum(ifelse(x > 0, x * log2(x), 0)))

	# Return
	o
}

isobreak <- function(i, x){
	# Split into segments
	x1 <- x[1:i]
	x2 <- x[i:length(x)]

	# Reverse second
	x2 <- -x2

	# Fit isotonic
	fit1 <- stats::isoreg(x1)
	fit2 <- stats::isoreg(x2)

	# Extract RSEM
	res1 <- sum(stats::residuals(fit1)^2)
	res2 <- sum(stats::residuals(fit2)^2)

	# Return
	res1 + res2
}

#' Shape statistic: Multimodality
#'
#' @param x numeric Rle vector: Coverage series.
#'
#' @return Numeric.
#' @family Shape functions
#' @import S4Vectors
#' @export
shapeMultimodality <- function(x){
	# Convert from Rle
	x <- as.vector(x)

	# Isobreak regression for all breakpoints
	o <- vapply(X=seq_along(x), FUN=isobreak, FUN.VALUE=numeric(1), x=x)

	# Best solution
	o <- min(o)

	# Return
	o
}
