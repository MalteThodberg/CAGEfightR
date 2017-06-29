#' Complete basic CAGE analysis
#'
#' This functions wraps a number of functions: readCTSS, coverageOfCTSS, tuneTagclusting, tagClustering & quantifyFeatures. Returns results as a RangedSummarizedExperiment
#'
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param design DataFrame or data.frame: Information about individual samples.
#' @param genome Seqinfo or NULL: If not NULL, use the supplied Seqinfo as genome.
#' @param keepCoverage logical: Whether to keep global coverage in the output object.
#' @param mergeDist integer: Merge TCs within this distance.
#' @param stepSize integer: Size of steps for increasing threshold.
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#'
#' @return RangedSummarizedExperiment.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import GenomicRanges SummarizedExperiment
#' @export
quickCAGE <- function(bwPlus, bwMinus, design, genome=NULL, keepCoverage=FALSE, mergeDist=20, stepSize=10, biocParallel=bpparam()){
	### Pre-checks
	# Classes
	stopifnot(class(bwPlus) == "BigWigFileList",
						class(bwMinus) == "BigWigFileList",
						class(design) == "DataFrame" | is.data.frame(design))
	stopifnot(identical(length(bwPlus), length(bwMinus), nrow(design)))

	# Rownames to enable matching
	stopifnot(!is.null(names(bwPlus)),
						!is.null(names(bwMinus)),
						!is.null(rownames(design)))

	stopifnot(all(names(bwPlus) == names(bwMinus)))
	stopifnot(all(rownames(design) == names(bwMinus)))

	### Wrap functions

	# Read in CTSS files
	message("### Running readCTSS...")
	ctss_grl <- readCTSS(bwPlus=bwPlus, bwMinus=bwMinus, genome = genome,
											 outFormat = "GRangesList", biocParallel = biocParallel)

	# Coverage
	message("\n### Running coverageOfCTSS...")
	global_cov <- coverageOfCTSS(grl=ctss_grl)

	# Tune
	message("\n### Running tuneTagClustering...")
	bestTC <- tuneTagClustering(global_cov, mergeDist=mergeDist, stepSize=stepSize)
	#optimalTC <- min(subset(bestTC, nTCs == max(nTCs), select=threshold))
	optimalTC <- min(bestTC[bestTC$nTCs == max(bestTC$nTCs),])

	# TCs
	message("\n### Running tagClustering...")
	TCs <- tagClustering(ctssCoverage=global_cov,
											 tpmCutoff=optimalTC,
											 mergeDist=mergeDist)

	# Quantify TCs
	message("\n### Running quantifyFeatures...")
	CM <- quantifyFeatures(ctss=ctss_grl, features=TCs,
												 biocParallel = biocParallel)

	### Assemble RSE
	message("\n### Assembling RangedSummarizedExperiment...")

	# As summarized experiment
	RSE <- SummarizedExperiment(assays=list(counts=CM),
																	rowRanges=TCs,
																	colData=design)

	# Stamp with TPM cutoff
	metadata(RSE)$optimalTPM <- optimalTC

	# Stamp with metaData
	metadata(RSE)$user <- Sys.getenv("USER")
	metadata(RSE)$data <- date()
	metadata(RSE)$sessionInfo <- utils::sessionInfo()


	# Check if coverage should be kept in output
	if(keepCoverage){
		metadata(RSE)$globalCoverage <- global_cov
		stopifnot(identical(seqinfo(RSE), seqinfo(global_cov)))
	}

	### Post-checks

	# Post-checks
	stopifnot(class(RSE) == "RangedSummarizedExperiment",
						!is.null(rownames(RSE)),
						!is.null(colnames(RSE)))

	# Return
	RSE
}
