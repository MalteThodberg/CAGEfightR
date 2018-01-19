#' Calculate support of CAGE data.
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS, cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputColumn character: Name of column in rowRanges to hold support values.
#' @param unexpressed numeric: Support will be calculated based on features larger than this cutoff.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @import assertthat SummarizedExperiment
#' @return RangedSummarizedExperiment with support added as a column in rowRanges.
#' @export
calcSupport <- function(object, inputAssay="counts", outputColumn="support", unexpressed=0){
	# Prechecks
	assert_that(methods::is(object, "SummarizedExperiment"),
							is.string(inputAssay),
							inputAssay %in% assayNames(object),
							is.string(outputColumn),
							is.number(unexpressed))

	if(outputColumn %in% colnames(rowData(object))){
		warning("object already has a column named ", outputColumn," in rowData: It will be overwritten!")
	}

	# Calculate support
	rowData(object)[,outputColumn] <- as.integer(Matrix::rowSums(assay(object, inputAssay) > unexpressed))

	# Return
	object
}

#' Calculate the total number of CAGE tags across samples.
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS, cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputColumn character: Name of column in colData to hold number of total tags.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @import assertthat SummarizedExperiment
#' @return RangedSummarizedExperiment with total tags added as a column in colData.
#' @export
calcTotalTags <- function(object, inputAssay="counts", outputColumn="totalTags"){
	# Prechecks
	assert_that(class(object) == "RangedSummarizedExperiment",
							not_empty(object),
							inputAssay %in% assayNames(object),
							is.string(inputAssay),
							is.string(outputColumn))

	if(outputColumn %in% colnames(colData(object))){
		warning("object already has a column named ", outputColumn," in colData: It will be overwritten!")
	}

	# Calculate colSums
	colData(object)[,outputColumn] <- Matrix::colSums(assay(object, inputAssay))

	# Return
	object
}

#' Calculate CAGE Tags-Per-Million (TPM)
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS, cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputAssay character: Name of assay to hold TPM values.
#' @param totalTags character or NULL: Column in colData holding the total number of tags for each samples. If NULL, this will be calculated using calcTotalTags.
#' @param outputColumn character: Name of column in colData to hold number of total tags, only used if totalTags is NULL.
#'
#' @return RangedSummarizedExperiment with TPM-values as a new assay. If totalTags is NULL, total tags added as a column in colData.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @import assertthat SummarizedExperiment
#' @export
calcTPM <- function(object, inputAssay="counts", outputAssay="TPM", totalTags=NULL, outputColumn="totalTags"){
	# Prechecks
	assert_that(class(object) == "RangedSummarizedExperiment",
							not_empty(object),
							is.string(inputAssay),
							inputAssay %in% assayNames(object),
							is.string(outputAssay))

	if(is.null(totalTags)){
		message("Calculating library sizes...")
		object <- calcTotalTags(object=object, inputAssay=inputAssay, outputColumn=outputColumn)
		totalTags <- outputColumn
	}else if(is.string(totalTags)){
		message("Using supplied library sizes...")
		assert_that(totalTags %in% colnames(colData(object)),
								is.numeric(colData(object)[,totalTags]),
								all(colData(object)[,totalTags] >= 0))
	}else{
		stop("totalTags should NULL or a column in colData!")
	}

	if(outputAssay %in% assayNames(object)){
		warning("object already has an assay named ", outputAssay,": It will be overwritten!")
	}

	# Scale counts to TPM
	message("Calculating TPM...")
	assay(object, outputAssay) <- Matrix::t(Matrix::t(assay(object, inputAssay)) / (colData(object)[,totalTags] / 1e6))

	# Return
	object
}

#' Calculate pooled expression across all samples.
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS, cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputColumn character: Name of column in rowRanges to hold pooled expression.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @import assertthat SummarizedExperiment
#' @return RangedSummarizedExperiment with pooled expression added as a column in rowRanges.
#' @export
calcPooled <- function(object, inputAssay="TPM", outputColumn="score"){
	# Prechecks
	assert_that(class(object) == "RangedSummarizedExperiment",
							not_empty(object),
							is.string(inputAssay),
							inputAssay %in% assayNames(object),
							is.string(outputColumn))

	if(outputColumn %in% colnames(rowData(object))){
		warning("object already has a column named ", outputColumn," in rowData: It will be overwritten!")
	}

	# Calculate colSums
	rowData(object)[,outputColumn] <- Matrix::rowSums(assay(object, inputAssay))

	# Return
	object
}

#' Calculate composition of CAGE data.
#'
#' For every feature (e.g. TSSs) count in how many samples it is expressed above a certain fraction (e.g. 10 percent) within a grouping, usually genes. This counts is refered to as a composition value.
#'
#' @param object RangedSummarizedExperiment: CAGE data quantified at CTSS, cluster or gene-level.
#' @param inputAssay character: Name of assay holding input expression values.
#' @param outputColumn character: Name of column in rowRanges to hold composition values.
#' @param unexpressed numeric: Composition will be calculated based on features larger than this cutoff.
#' @param genes character: Name of column in rowData holding genes (NAs are not allowed.)
#'
#' @import assertthat S4Vectors SummarizedExperiment
#' @return RangedSummarizedExperiment with composition added as a column in rowData.
#' @export
calcComposition <- function(object, inputAssay="counts", outputColumn="composition", unexpressed=0.1, genes="geneID"){
	assert_that(methods::is(object, "SummarizedExperiment"),
							is.string(inputAssay),
							inputAssay %in% assayNames(object),
							is.string(outputColumn),
							is.number(unexpressed),
							unexpressed >= 0 & unexpressed <= 1,
							is.string(genes),
							is.element(genes, colnames(rowData(object))),
							is.character(rowData(object)[,genes]),
							noNA(rowData(object)[,genes]))

	if(outputColumn %in% colnames(rowData(object))){
		warning("object already has a column named ", outputColumn," in rowData: It will be overwritten!")
	}

	# Extract gene-wise matrices
	L <- splitAsList(x=assay(object, inputAssay), f=rowData(object)[,genes], drop = FALSE)

	# Scale and find high compositions
	L <- endoapply(L, function(x) scale(x, center=FALSE, scale=Matrix::colSums(x)) > unexpressed) # Scale always coerces to a base::matrix!

	# Calculate count high compositions
	L <- endoapply(L, rowSums, na.rm=TRUE)

	# Back to integer vector
	L <- unsplit(L, f=rowData(object)[,genes], drop = FALSE)
	L <- as.integer(L)

	# Post-checks
	stopifnot(length(L) == nrow(object),
						noNA(L))

	# Append
	rowData(object)[,outputColumn] <- L

	# Return
	object
}