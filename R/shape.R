#### Main shape function ####

#' Calculate Tag Cluster shapes
#'
#' Apply a shape-function to the pooled CTSS signal of every Tag Cluster (TC).
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: TCs.
#' @param pooled GenomicRanges or RangedSummarizedExperiment: Pooled CTSS as the
#'   score column.
#' @param outputColumn character: Name of column to hold shape statistics.
#' @param shapeFunction function: Function to apply to each TC (See details).
#' @param ... additional arguments passed to shapeFunction.
#'
#' @return object with calculated shape statistics added as a column in rowData
#'   (or mcols).
#' @family Calculation functions
#' @family Shape functions
#' @export
#' @examples
#' data(exampleCTSSs)
#' data(exampleUnidirectional)
#'
#' # Calculate pooled CTSSs using pre-calculated number of total tags:
#' exampleCTSSs <- calcTPM(exampleCTSSs, totalTags='totalTags')
#' exampleCTSSs <- calcPooled(exampleCTSSs)
#'
#' # Calculate shape statistics
#' calcShape(exampleUnidirectional, pooled=exampleCTSSs,
#'     outputColumn='entropy', shapeFunction=shapeEntropy)
#' calcShape(exampleUnidirectional, pooled=exampleCTSSs, outputColumn='IQR',
#'     shapeFunction=shapeIQR, lower=0.2, upper=0.8)
#'
#' # See the vignette for how to implement custom shape functions!
setGeneric("calcShape", function(object, pooled, ...) {
    standardGeneric("calcShape")
})

#' @rdname calcShape
setMethod("calcShape", signature(object = "GRanges", pooled = "GRanges"),
          function(object, pooled, outputColumn = "IQR",
                   shapeFunction = shapeIQR, ...) {
    # Pre-checks
    assert_that(checkPooled(pooled),
                is.character(outputColumn),
                is.function(shapeFunction),
                identical(seqlengths(object), seqlengths(pooled)))

    # Warnings
    if (outputColumn %in% colnames(mcols(object))) {
        warning("object already has a column named ",
                outputColumn, " in mcols: It will be overwritten!")
    }

    # Names need to be set for sorting later
    if (is.null(names(object))) {
        message("Adding names...")
        names(object) <- paste("TC", seq_along(object))
    }

    # Split by strand
    message("Splitting by strand...")
    covByStrand <- splitPooled(pooled)
    TCsByStrand <- splitByStrand(object)

    #covByStrand <- splitByStrand(pooled)
    #TCsByStrand <- splitByStrand(object)

    # Coverage by strand
    #coverage_plus <- coverage(covByStrand$`+`, weight = "score")
    #coverage_minus <- coverage(covByStrand$`-`, weight = "score")
    #rm(covByStrand)

    # Views
    message("Applying function to each cluster...")
    views_plus <- Views(covByStrand$`+`, methods::as(TCsByStrand$`+`, 'IRangesList'))
    views_minus <- Views(covByStrand$`-`, methods::as(TCsByStrand$`-`, 'IRangesList'))
    rm(covByStrand)

    # Tmp solution circumventing direct use of RangesList
    #views_plus <- Views(coverage_plus, split(ranges(TCsByStrand$`+`), seqnames(TCsByStrand$`+`)))
    #views_minus <- Views(coverage_minus, split(ranges(TCsByStrand$`-`), seqnames(TCsByStrand$`-`)))

    # Applying functions to views
    stat_plus <- viewApply(views_plus, shapeFunction, ...)
    stat_minus <- viewApply(views_minus, shapeFunction, ...)

    # Reset names
    message("Preparing output output...")
    stat_plus <- as.numeric(unlist(stat_plus, use.names = FALSE))
    names(stat_plus) <- names(TCsByStrand$`+`)
    stat_minus <- as.numeric(unlist(stat_minus, use.names = FALSE))
    names(stat_minus) <- names(TCsByStrand$`-`)
    rm(TCsByStrand)

    # Reassemble in same order
    o <- c(stat_plus, stat_minus)
    o <- o[match(names(object), names(o))]
    names(o) <- NULL
    rm(stat_plus, stat_minus)

    # Add to object
    mcols(object)[, outputColumn] <- o

    # Return
    object
})

#' @rdname calcShape
setMethod("calcShape", signature(object = "RangedSummarizedExperiment", pooled = "GRanges"),
    function(object, pooled, ...) {
        rowRanges(object) <- calcShape(rowRanges(object), pooled, ...)
        object
    })

#' @rdname calcShape
setMethod("calcShape", signature(object = "GRanges", pooled = "RangedSummarizedExperiment"),
    function(object, pooled, ...) {
        calcShape(object, rowRanges(pooled), ...)
    })

#' @rdname calcShape
setMethod("calcShape", signature(object = "GRanges", pooled = "GPos"),
					function(object, pooled, ...) {
						calcShape(object, methods::as(pooled, "GRanges"), ...)
					})

#' @rdname calcShape
setMethod("calcShape", signature(object = "RangedSummarizedExperiment", pooled = "RangedSummarizedExperiment"),
    function(object, pooled, ...) {
        rowRanges(object) <- calcShape(rowRanges(object), rowRanges(pooled), ...)
        object
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
#' @export
#' @examples
#' # Hypothetical shard/broad clusters:
#' x_sharp <- Rle(c(1,1,1,4,5,2,1,1))
#' x_broad <- Rle(c(1,2,3,5,4,3,2,1))
#'
#' # Calculate IQR
#' shapeIQR(x_sharp)
#' shapeIQR(x_broad)
#'
#' # See calcShape for more usage examples
shapeIQR <- function(x, lower = 0.25, upper = 0.75) {
    # To normal vector
    x <- as.vector(x)

    # Scale by sum
    x <- x/sum(x)

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
#' Calculates the Shannon Entropy (base log2) for a vector. Zeros are removed
#' before calculation.
#'
#' @param x numeric Rle vector: Coverage series.
#'
#' @return Numeric.
#' @family Shape functions
#' @export
#' @examples
#' # Hypothetical shard/broad clusters:
#' x_sharp <- Rle(c(1,1,1,4,5,2,1,1))
#' x_broad <- Rle(c(1,2,3,5,4,3,2,1))
#'
#' # Calculate Entropy
#' shapeEntropy(x_sharp)
#' shapeEntropy(x_broad)
#'
#' # See calcShape for more usage examples
shapeEntropy <- function(x) {
    # To normal vector x <- as.vector(x[x > 0])
    x <- as.vector(x)

    # Scale by sum
    x <- x/sum(x)

    # Calculate entropy
    o <- suppressWarnings(-sum(ifelse(x > 0, x * log2(x), 0)))

    # Return
    o
}

isobreak <- function(i, x) {
    # Split into segments x1 <- x[1:i]
    x1 <- x[seq_len(i)]
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
#' @examples
#' # See calcShape for usage examples
shapeMultimodality <- function(x) {
    # Convert from Rle
    x <- as.vector(x)

    # Isobreak regression for all breakpoints
    o <- vapply(X = seq_along(x), FUN = isobreak, FUN.VALUE = numeric(1), x = x)

    # Best solution
    o <- min(o)

    # Return
    o
}
