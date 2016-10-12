#' Quantify CTSS across features
#'
#' Counts the number of CTSS-tags in a set of features
#'
#' @param ctss GRangesList or character: CTSS data stores as a GRangesList or a vector of paths to CTSS-files.
#' @param features GRanges: Ranges to be quantified.
#' @param cores NULL or integer: If not NULL, use multiple cores to quantify features.
#'
#' @return matrix with rows as features and columns as samples.
#'
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
quantifyFeatures <- function(ctss, features, cores) UseMethod("quantifyFeatures")

#' @describeIn quantifyFeatures default method (Always throws an error)
#' @export
quantifyFeatures.default <- function(ctss, features, cores){
	stop("ctss must be either a character or GRangesList")
}

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
#' @import S4Vectors IRanges GenomicRanges
#' @export
countScoredOverlaps <- function(gr, features){
		# Solution from bioconductor by M Morgan
		hits <- as(findOverlaps(query=features, subject=gr), "List")
		weightedCount <- sum(extractList(score(gr), hits))

		# Return
		weightedCount
}

#' @describeIn quantifyFeatures high memory method
#' @export
quantifyFeatures.GRangesList <- function(ctss, features, cores=NULL){
	# Get some basic info for printing
	nFeatures <- length(features)
	nSamples <- length(ctss)
	m <- sprintf("Quantifying expression of %d features in %d samples",
							 nFeatures, nSamples)
	message(m)

	# Quantify each GRange in GRanges list
	if(is.null(cores)){
		x <- pbapply::pblapply(ctss, countScoredOverlaps, features=features)
	}else{
		message(paste0("Using cores: ", cores))
		x <- parallel::mclapply(ctss, countScoredOverlaps, features=features,
									mc.preschedule=FALSE, mc.cores=cores)
	}

	# Merge into matrix
	message("Building Expression Matrix")
	x <- do.call(cbind, x)
	rownames(x) <- names(features)

	# Return
	x
}

#' @describeIn quantifyFeatures low memory method
#' @export
quantifyFeatures.character <- function(ctss, features, cores=NULL){
	# Multiple cores not allowed
	if(!is.null(cores)){
		stop("Low memory quantification does not currently support using multiple cores")
	}

	# Quantify each GRange in GRanges list
	tmp_fun <- function(x, f){
		# Load file
		d <- fastReadBED(x)

		# Turn into GRanges
		d <- makeGRangesFromDataFrame(df=d,
																	keep.extra.columns=TRUE,
																	starts.in.df.are.0based=TRUE)

		# Overlap
		d <- countScoredOverlaps(gr=d, features=f)

		# Return
		d
	}

	x <- pbapply::pblapply(ctss, tmp_fun, f=features)

	# Merge into matrix
	x <- do.call(cbind, x)
	rownames(x) <- names(features)

	# Return
	x
}
