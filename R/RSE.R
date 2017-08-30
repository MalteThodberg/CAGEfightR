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
	### Mark deprecation
	.Deprecated(c("quickTSSs", "quickEnhancers"))

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

#' @import GenomicRanges SummarizedExperiment
quickSetter <- function(design, bwPlus, bwMinus, ctss, ctssCoverage, genome, biocParallel){
	### Pre checks
	stopifnot(class(design) == "DataFrame",
						!is.null(rownames(design)))

	### Read in samples
	if(any(!is.null(bwPlus), !is.null(bwMinus))){ # BigWigs

		stopifnot(class(bwPlus) == "BigWigFileList",
							class(bwMinus) == "BigWigFileList",
							identical(length(bwPlus), length(bwMinus), nrow(design)),
							!is.null(names(bwPlus)),
							!is.null(names(bwMinus)),
							all(names(bwPlus) == names(bwMinus)),
							all(rownames(design) == names(bwMinus)))

		message("### Running readCTSS...")
		ctss <- readCTSS(bwPlus=bwPlus,
										 bwMinus=bwMinus,
										 genome=genome,
										 outFormat="GRangesList",
										 biocParallel=biocParallel)
	}else if(!is.null(ctss)){ # GRangesList
		message("Using supplied CTSS...")
		stopifnot(class(ctss) == "GRangesList",
							length(ctss) == nrow(design),
							!is.null(names(ctss)),
							all(rownames(design) == names(ctss)))
	}else{
		stop("CTSS must be supplied as two BigWigFileList or a GRangesList")
	}

	### Coverage
	if(is.null(ctssCoverage)){
		message("\n### Running coverageOfCTSS...")
		ctssCoverage <- coverageOfCTSS(grl=ctss)
	}else{
		message("Using supplied CTSS coverage...")
		stopifnot(class(ctssCoverage) == "GRanges")
	}

	### Return
	list(ctss=ctss, ctssCoverage=ctssCoverage)
}

#' Complete pipeline for finding and quantifying TSSs.
#'
#' This functions wraps a number of functions: readCTSS, coverageOfCTSS,
#' tuneTagclusting, tagClustering & quantifyFeatures. Returns results as a
#' RangedSummarizedExperiment
#'
#' @param design DataFrame or data.frame: Information about individual samples.
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param ctss GRangesList: CTSS data for samples.
#' @param ctssCoverage GRanges: CTSS-coverage across all samples.
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
quickTSSs <- function(design, bwPlus=NULL, bwMinus=NULL, ctss=NULL, ctssCoverage=NULL, genome=NULL, keepCoverage=FALSE, mergeDist=20, stepSize=10, biocParallel=bpparam()){
	### Prechecks
	stopifnot(is.logical(keepCoverage),
						mergeDist >= 0,
						stepSize > 0)

	### Run the quicksetter
	tmp <- quickSetter(design=design, bwPlus=bwPlus, bwMinus=bwMinus, ctss=ctss, ctssCoverage=ctssCoverage, genome=genome, biocParallel=biocParallel)
	ctss <- tmp$ctss
	ctssCoverage <- tmp$ctssCoverage
	rm(tmp)

	### Tag Clustering

	# Tune
	message("\n### Running tuneTagClustering...")
	bestTC <- tuneTagClustering(ctssCoverage, mergeDist=mergeDist, stepSize=stepSize)
	optimalTC <- min(bestTC[bestTC$nTCs == max(bestTC$nTCs),])

	# TCs
	message("\n### Running tagClustering...")
	TCs <- tagClustering(ctssCoverage=ctssCoverage,
											 tpmCutoff=optimalTC,
											 mergeDist=mergeDist)

	# Quantify TCs
	message("\n### Running quantifyFeatures...")
	CM <- quantifyFeatures(ctss=ctss,
												 features=TCs,
												 biocParallel=biocParallel)

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
		metadata(RSE)$ctssCoverage <- ctssCoverage
		stopifnot(identical(seqinfo(RSE), seqinfo(ctssCoverage)))
	}

	### Post-checks

	# Post-checks
	stopifnot(class(RSE) == "RangedSummarizedExperiment",
						!is.null(rownames(RSE)),
						!is.null(colnames(RSE)))

	# Return
	RSE
}

#' Complete pipeline for finding and quantifying enhancers.
#'
#' @param design DataFrame or data.frame: Information about individual samples.
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param ctss GRangesList: CTSS data for samples.
#' @param ctssCoverage GRanges: CTSS-coverage across all samples.
#' @param genome Seqinfo or NULL: If not NULL, use the supplied Seqinfo as genome.
#' @param keepCoverage logical: Whether to keep global coverage in the output object.
#' @param window integer: Width of sliding window used for calculating sums.
#' @param balanceThreshold numeric: Minimum balance score to be considered a bidirectional site.
#' @param minBidirectionality integer: Minimum number of samples that show bidirectional transcription.
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#'
#' @return RangedSummarizedExperiment.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import GenomicRanges SummarizedExperiment
#' @export
quickEnhancers <- function(design, bwPlus=NULL, bwMinus=NULL, ctss=NULL, ctssCoverage=NULL, genome=NULL, keepCoverage=FALSE, window=199, balanceThreshold=0.95, minBidirectionality=1, biocParallel=bpparam()){
	### Prechecks
	stopifnot(is.logical(keepCoverage))

	### Run the quicksetter
	tmp <- quickSetter(design=design, bwPlus=bwPlus, bwMinus=bwMinus,
										 ctss=ctss, ctssCoverage=ctssCoverage,
										 genome=genome, biocParallel=biocParallel)

	ctss <- tmp$ctss
	ctssCoverage <- tmp$ctssCoverage
	rm(tmp)

	### Tag Clustering

	# Tune
	message("\n### Running bidirectionalClustering...")
	BCs <- bidirectionalClustering(ctssCoverage=ctssCoverage, window=window, balanceThreshold=balanceThreshold)

	message("\n### Running calculateBidirectionality...")
	BCs$bidirectionality <- calculateBidirectionality(ctss=ctss, features=BCs, biocParallel=biocParallel)
	BCs <- subset(BCs, bidirectionality >= minBidirectionality)

	# Quantify TCs
	message("\n### Running quantifyFeatures...")
	CM <- quantifyFeatures(ctss=ctss,
												 features=BCs,
												 biocParallel=biocParallel)

	### Assemble RSE
	message("\n### Assembling RangedSummarizedExperiment...")

	# As summarized experiment
	RSE <- SummarizedExperiment(assays=list(counts=CM),
															rowRanges=BCs,
															colData=design)

	# Stamp with metaData
	metadata(RSE)$user <- Sys.getenv("USER")
	metadata(RSE)$data <- date()
	metadata(RSE)$sessionInfo <- utils::sessionInfo()

	# Check if coverage should be kept in output
	if(keepCoverage){
		metadata(RSE)$ctssCoverage <- ctssCoverage
		stopifnot(identical(seqinfo(RSE), seqinfo(ctssCoverage)))
	}

	### Post-checks

	# Post-checks
	stopifnot(class(RSE) == "RangedSummarizedExperiment",
						!is.null(rownames(RSE)),
						!is.null(colnames(RSE)))

	# Return
	RSE
}

#' Stack two RangedSummarizedExperiments
#'
#' A convienient wrapper around rbind for stacking rows RangedSummarizedExperiment-objecs (RSEs). Automatically keeps only shared rowData and assays, and has the option of discarding overlapping features. Useful for combing RSEs containing TSSs and enhancers.
#'
#' @param a RangedSummarizedExperiment: RSE to be stacked.
#' @param b RangedSummarizedExperiment: RSE to be stacked.
#' @param removeIfOverlapping character: Whether to discard features in either "a" or "b" if they overlap. "both" does not remove overlapping features.
#'
#' @return Stacked and sorted RangedSummarizedExperiment.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors GenomicRanges SummarizedExperiment
#' @export
stackRSEs <- function(a, b, removeIfOverlapping="both"){
	# Check if RSEs are combatible
	stopifnot(class(a) == "RangedSummarizedExperiment",
						class(b) == "RangedSummarizedExperiment",
						!is.null(rownames(a)),
						!is.null(rownames(b)),
						all(colnames(a) == colnames(b)),
						identical(seqinfo(a), seqinfo(b)),
						identical(colData(a), colData(b)),
						removeIfOverlapping %in% c("a", "b", "both"))

	# Remove overlapping features
	if(removeIfOverlapping == "a"){
		# Find overlaps
		os <- overlapsAny(a, b)

		# Remove
		message("Removing overlapping features from a: ", sum(os))
		a <- subset(a, !os)
	}else if(removeIfOverlapping == "b"){
		# Find overlaps
		os <- overlapsAny(b, a)

		# Remove
		message("Removing overlapping features from b: ", sum(os))
		b <- subset(b, !os)
	}else{
		os <- intersect(rowRanges(a), rowRanges(b))
		message("Overlapping features were not removed: ", length(os))
	}
	rm(os)

	# Only retain shared assays
	shared_assays <- intersect(assayNames(a), assayNames(b))
	message("Keeping assays: ", paste(shared_assays, collapse=", "))
	assays(a) <- assays(a)[shared_assays]
	assays(b) <- assays(b)[shared_assays]
	rm(shared_assays)

	# Only retain shared rowcolumns
	shared_cols <- intersect(colnames(mcols(a)), colnames(mcols(b)))
	message("Keeping columns: ", paste(shared_cols, collapse=", "))
	mcols(a) <- mcols(a)[,shared_cols]
	mcols(b) <- mcols(b)[,shared_cols]
	rm(shared_cols)

	# Attempt to merge metadata
	message("Merging metadata...")
	new_meta <- c(metadata(a), metadata(b))
	metadata(a) <- list()
	metadata(b) <- list()
	new_meta <- new_meta[!duplicated(new_meta)]

	# Bind, add meta and re-sort
	message("Stacking RSEs and re-sorting...")
	o <- rbind(a, b)
	metadata(o) <- new_meta
	rm(new_meta)
	o <- sort(o)

	# Return
	o
}
