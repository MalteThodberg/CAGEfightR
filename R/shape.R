#' Calculate Shannon Entropy
#'
#' Calculates the Shannon Entropy (base log2) for a vector. Zeros are removed before calculation.
#'
#' @param x numeric Rle vector: Coverage series.
#'
#' @return Numeric.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Shape functions
#' @export
ShannonEntropy <- function(x){
	# To normal vector
	#x <- as.vector(x[x > 0])
	x <- as.vector(x)

	# Scale by sum
	x <- x / sum(x)

	# Calculate entropy
	o <- -sum(ifelse(x > 0, x * log2(x), 0))
	#o <- -sum(x * log2(x))

	# Return
	o
}

#' Calculate InterQuartile Range
#'
#' Calculates the interquartile range of a vector.
#'
#' @param x numeric Rle vector: Coverage series.
#' @param lower numeric: Lower quartile.
#' @param upper numeric: Upper quartile.
#'
#' @return Numeric
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Shape functions
#' @export
InterQuartileRange <- function(x, lower=0.25, upper=0.75){
	# To normal vector
	x <- as.vector(x)

	# Scale by sum
	x <- x / sum(x)

	# Cumulate sum
	x <- cumsum(x)

	# Find the threshold
	lowerPos <- Position(function(y) y >= lower, x)
	upperPos <- Position(function(y) y >= upper, x) + 1

	# Return difference
	upperPos - lowerPos
}

#' Round Rle vector to integers while preserving sum
#'
#' Round an Rle vector, while preserving the sum of the vector as closely as possible. Slightly slower than runningRounding. Can be used with the HartigansDip function.
#'
#' @param x numeric Rle vector: Coverage series.
#'
#' @return Integer vector.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @references \url{https://stackoverflow.com/a/32544987}
#' @import S4Vectors
#' @family Shape functions
#' @export
sumRounding <- function(x) {
	x <- as.vector(x)
	y <- floor(x)
	indices <- tail(order(x-y), round(sum(x)) - sum(y))
	y[indices] <- y[indices] + 1
	y
}

#' Round Rle vector to integers while preserving running sum
#'
#' Round an Rle vector, while preserving the running sum. Slightly faster than sumRounding. Can be used with the HartigansDip function.
#'
#' @param x numeric Rle vector: Coverage series.
#'
#' @return Integer vector.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @references \url{https://stackoverflow.com/a/36970695}
#' @import S4Vectors
#' @export
runningRounding <- function(x) {
	x <- as.vector(x)
	diff(c(0, round(cumsum(x))))
}

#' Round Rle vector to integers
#'
#' Wrapper for round(). Can be used with the HartigansDip function.
#'
#' @param x numeric Rle vector.
#'
#' @return Integer vector.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors
#' @export
simpleRounding <- function(x){
	round(as.vector(x))
}

#' Hartigans' dip statistics
#'
#' Calculates the Hartigans Dip Statistics for a vector, by first rounding values to integers. A dip statistics of 0 indicates a unimodal distribtuion and a dip statistics of 1 indicates a bi- or multimodal distribution.
#'
#' @param x numeric Rle vector: Coverage series.
#' @param roundingFun function: Function for rounding an Rle vector.
#'
#' @return numeric.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
HartigansDip <- function(x, roundingFun=runningRounding){
	# Round
	x <- roundingFun(x)

	# Expand
	x <-  rep(1:length(x), x )

	# Dip
	diptest::dip(x, full.result=FALSE, min.is.0=TRUE)
}

#' Calculate TC shape
#'
#' Apply a shape-function to the global coverage of every TC.
#'
#' @param ctssCoverage GRanges: GRanges with CTSS as score column.
#' @param TCs GRanges: TCs to be analyzed.
#' @param shapeFunction function: Function to apply to each TC (See details).
#' @param ... additional arguments passed to shapeFunction.
#'
#' @return Numeric vector of same length as TCs holding output from shapeFunction.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Shape functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
calculateShape <- function(ctssCoverage, TCs, shapeFunction, ...){
	# Names are used for merging
	stopifnot(!is.null(names(TCs)))

	# Split by strand
	message("Splitting coverage by strand...")
	covByStrand <- splitByStrand(ctssCoverage)
	TCsByStrand <- splitByStrand(TCs)

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
	o <- o[match(names(TCs), names(o))]
	names(o) <- NULL

	# Return
	o
}
