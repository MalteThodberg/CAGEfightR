#' Subset by support across samples
#'
#' A convienient wrapper around calcSupport and subset.
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS,
#'   cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputColumn character: Name of column in rowRanges to hold support
#'   values.
#' @param unexpressed numeric: Support will be calculated based on features
#'   larger than this cutoff.
#' @param minSamples numeric: Only features with support in more than this
#'   number of samples will be kept.
#'
#' @return RangedSummarizedExperiment with support added as a column in
#'   rowRanges and features with less support than minSamples removed.
#'
#' @family Subsetting functions
#' @family Calculation functions
#'
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' data(exampleBidirectional)
#'
#' # Keep clusters with at least one tag in two samples
#' subsetBySupport(exampleBidirectional)
#'
#' # Keep clusters with at least two tags in four samples
#' subsetBySupport(exampleBidirectional, unexpressed=1, minSamples=2)
subsetBySupport <- function(object, inputAssay="counts", outputColumn="support", unexpressed=0, minSamples=1){
	# Pre-checks
	assert_that(is.numeric(minSamples))

	# Call function
	message("Calculating support...")
	object <- calcSupport(object=object, inputAssay=inputAssay, outputColumn=outputColumn, unexpressed=unexpressed)
	before <- nrow(object)

	# Subset
	message("Subsetting...")
	object <- object[rowData(object)[,outputColumn] > minSamples,]
	after <- nrow(object)
	removed <- before-after

	# Print some info
	message("Removed ", removed, " out of ", before, " regions (", round(removed/before*100, digits=1), "%)")

	# Return
	object
}

#' Subset by composition across samples
#'
#' A convenient wrapper around calcComposition and subset.
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS,
#'   cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputColumn character: Name of column in rowRanges to hold
#'   composition values.
#' @param unexpressed numeric: Composition will be calculated based on features
#'   larger than this cutoff.
#' @param genes character: Name of column in rowData holding genes (NAs are not
#'   allowed.)
#' @param minSamples numeric: Only features with composition in more than this
#'   number of samples will be kept.
#'
#' @return RangedSummarizedExperiment with composition values added as a column
#'   in rowData and features with less composition than minSamples removed.
#' @family Subsetting functions
#' @family Calculation functions
#'
#' @export
#' @examples
#' data(exampleUnidirectional)
#'
#' # Annotate clusters with geneIDs:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'
#' exampleUnidirectional <- assignGeneID(exampleUnidirectional,
#'                                       geneModels=txdb,
#'                                       outputColumn="geneID")
#' exampleUnidirectional <- subset(exampleUnidirectional, !is.na(geneID))
#'
#' # Keep only clusters more than 10% in more than one sample:
#' calcComposition(exampleUnidirectional)
#'
#' # Keep only clusters more than 5% in more than 2 samples:
#' subsetByComposition(exampleUnidirectional, unexpressed = 0.05, minSamples=2)
subsetByComposition <- function(object, inputAssay="counts", outputColumn="composition", unexpressed=0.1, genes="geneID", minSamples=1){
	# Pre-checks
	assert_that(is.numeric(minSamples))

	# Call function
	message("Calculating composition...")
	object <- calcComposition(object=object, inputAssay=inputAssay, outputColumn=outputColumn, unexpressed=unexpressed, genes=genes)
	before <- nrow(object)

	# Subset
	message("Subsetting...")
	object <- object[rowData(object)[,outputColumn] > minSamples,]
	after <- nrow(object)
	removed <- before-after

	# Print some info
	message("Removed ", removed, " out of ", before, " regions (", round(removed/before*100, digits=1), "%)")

	# Return
	object
}
