#' Swap ranges in a GRanges.
#'
#' Swap out the range of a GRanges-object with another IRanges-object stored inside the same object. I.e., swapping cluster widths with cluster peaks.
#'
#' @param object GRanges or RangedSummarizedExperiment: Primary ranges to be swapped out.
#' @param inputColumn character: Name of column holding IRanges to be swapped in.
#' @param outputColumn character or NULL: Name of column to hold swapped out ranges, if NULL original ranges are not saved.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with inputColumn swapped in as ranges.
#'
#' @family Swapping functions
#' @export
setGeneric("swapRanges", function(object, ...){
	standardGeneric("swapRanges")
})

#' @import assertthat S4Vectors IRanges GenomicRanges
#' @rdname swapRanges
setMethod("swapRanges", signature(object="GenomicRanges"), function(object, inputColumn="thick", outputColumn=NULL){
	# Pre-checks
	assert_that(is.string(inputColumn),
							is.element(inputColumn, colnames(mcols(object))),
							class(mcols(object)[,inputColumn]) == "IRanges",
							is.string(outputColumn) | is.null(outputColumn))

	# Warnings
	if(!is.null(outputColumn)){
		if(outputColumn %in% colnames(mcols(object))){
			warning("object already has a column named ", outputColumn," in mcols: It will be overwritten!")
		}
	}

	# Save old range as new column
	mcols(object)[,outputColumn] <- ranges(object)

	# Switch in the new column
	ranges(object) <- mcols(object)[,inputColumn]

	# Use original names
	names(object) <- names(mcols(object)[,outputColumn])

	# Return
	object
})

#' @import SummarizedExperiment
#' @rdname swapRanges
setMethod("swapRanges", signature(object="RangedSummarizedExperiment"), function(object, ...){
	rowRanges(object) <- swapRanges(rowRanges(object), ...)
	object
})

#' Swap scores in SummarizedExperiment
#'
#' Take scores for a specific sample and a specific assay and put them into rowData.
#'
#' @param object SummarizedExperiment: CAGE-data
#' @param outputColumn character: Column in rowData to to hold swapped in scores.
#' @param inputAssay character: Name of assay to take scores from.
#' @param sample character: Name of sample to take scores from.
#'
#' @return SummarizedExperiment with sample scores from inputAssay in rowRata.
#' @import assertthat SummarizedExperiment
#' @export
swapScores <- function(object, outputColumn="score", inputAssay, sample){
	assert_that(methods::is(object, "SummarizedExperiment"),
							is.string(outputColumn),
							is.string(inputAssay),
							is.element(inputAssay, assayNames(object)),
							is.string(sample),
							is.element(sample, colnames(object)))

	# Warnings
	if(outputColumn %in% colnames(rowData(object))){
		warning("object already has a column named ", outputColumn," in rowData: It will be overwritten!")
	}

	# Swap in new column from assay
	rowData(object)[,outputColumn] <- assay(object, inputAssay)[,sample]

	# Return
	object
}
