#' Calculate TPM.
#'
#' Calculate Tags-Per-Million (TPM) for a numeric vector.
#'
#' @param x Integer: Counts to be TPM normalized.
#'
#' @return Numeric vector of same length as x.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Helper functions
#' @export
TPM <- function(x){
	x / (sum(x) / 1e6)
}

#' Split GRanges by strand.
#'
#' Split a GRanges into a GRangesList by strand.
#'
#' @param gr GRanges-object.
#'
#' @return GRangesList with three elements: +, - and *.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Helper functions
#' @import S4Vectors GenomicRanges
#' @export
splitByStrand <- function(gr){
	split(gr, strand(gr))
}

#' Basic preprocessing of CTSS-files.
#'
#' Function that trims away weakly expressed CTSSs by discarding based on CTSS counts and/or TPM-values, before performing a final TPM-calculation.
#'
#' @param gr GRanges: GRanges holding CTSS counts in the score column.
#' @param preFilterCTSS Integer: Minimum count for retained CTSSs.
#' @param preFilterTPM Numeric: Minimum TPM-value for retained CTSSs.
#'
#' @note The order of trimming is: Optional trimming on CTSS counts, optional trimming of TPM-values, final calculation of TPM.
#'
#' @return GRanges with TPM-values in the score column.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Coverage functions
#' @import GenomicRanges
#' @export
trimAndNormalizeGR <- function(gr, preFilterCTSS=NULL, preFilterTPM=NULL){
	# Filter on CTSS
	if(!is.null(preFilterCTSS)){
		gr <- subset(gr, score >= preFilterCTSS)
	}

	# Filter on TPM
	if(!is.null(preFilterTPM)){
		gr <- subset(gr, TPM(score) >= preFilterTPM)
	}

	# TPM normalize
	score(gr) <- TPM(score(gr))

	# Return
	gr
}

#' Calculate global CTSS coverage.
#'
#' Given a series of GRanges, preprocesses and normalizes the CTSS-counts before calculating the global coverage by summing over all samples.
#'
#' @param grl GRangesList or SimpleList: GRanges with CTSSs in the score column.
#' @param normalizationFunction function or NULL: Function to preprocess the individal GRanges (See Details), or NULL if data should be considered as already normalized.
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#' @param ... Additional arguments passed to normalizationFunction.
#'
#' @return GRanges with global summed-coverage as the score column.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Coverage functions
#' @import S4Vectors IRanges GenomicRanges BiocParallel
#' @export
coverageOfCTSS <- function(grl, normalizationFunction=trimAndNormalizeGR, biocParallel=bpparam(), ...){
	# Checks
	# Should check if an apply loop produces only GRanges

	if(is.null(normalizationFunction)){
		message("Assuming CTSS is already preprocessed...")
		grl <- GRangesList(grl)
	}else{
		# Preprocess
		message("Preprocessing...")
		grl <- bplapply(grl, normalizationFunction, ..., BPPARAM=biocParallel)
		grl <- GRangesList(grl)
	}

	# Flatten and split by strand
	message("Splitting...")
	grs <- splitByStrand(unlist(grl))
	rm(grl)

	# Plus cov
	message("Coverage...")
	plus <- coverage(grs$`+`, weight="score")
	minus <- coverage(grs$`-`, weight="score")
	rm(grs)

	# Merge as GR
	message("Output...")
	o <- c(GRanges(plus, strand="+"), GRanges(minus, strand="-"))
	rm(plus, minus)

	# Subset out zero ranges
	o <- subset(o, score > 0)

	# Return
	o
}
