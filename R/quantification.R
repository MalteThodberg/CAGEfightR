# #' Quantify CTSS across features
# #'
# #' Counts the number of CTSS-tags in a set of features
# #'
# #' @param ctss GRangesList or character: CTSS data stores as a GRangesList or a vector of paths to CTSS-files.
# #' @param features GRanges: Ranges to be quantified.
# #' @param cores NULL or integer: If not NULL, use multiple cores to quantify features.
# #'
# #' @return matrix with rows as features and columns as samples.
# #'
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import S4Vectors IRanges GenomicRanges
# #' @export
# quantifyFeatures <- function(ctss, features, cores=NULL) UseMethod("quantifyFeatures")
#
# #' @describeIn quantifyFeatures default method (Always throws an error)
# #' @export
# quantifyFeatures.default <- function(ctss, features, cores=NULL){
# 	stop("ctss must be either a character or GRangesList")
# }
#
# #' Internal function: countOverlaps with scores
# #'
# #' Similar to countOverlaps, but takes into account the scores column.
# #'
# #' @param gr GRanges: Ranges with scores.
# #' @param features GRanges: Ranges to be quantified.
# #'
# #' @return numeric vector of same length as features.
# #'
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import S4Vectors IRanges GenomicRanges
# #' @export
# countScoredOverlaps <- function(gr, features){
# 		# Solution from bioconductor by M Morgan
# 		hits <- as(findOverlaps(query=features, subject=gr), "List")
# 		weightedCount <- sum(extractList(score(gr), hits))
#
# 		# Return
# 		weightedCount
# }
#
# #' @describeIn quantifyFeatures high memory method
# #' @export
# quantifyFeatures.GRangesList <- function(ctss, features, cores=NULL){
# 	# Get some basic info for printing
# 	nFeatures <- length(features)
# 	nSamples <- length(ctss)
# 	m <- sprintf("Quantifying expression of %d features in %d samples",
# 							 nFeatures, nSamples)
# 	message(m)
#
# 	# Quantify each GRange in GRanges list
# 	if(is.null(cores)){
# 		x <- pbapply::pblapply(ctss, countScoredOverlaps, features=features)
# 	}else{
# 		message(paste0("Using cores: ", cores))
# 		x <- parallel::mclapply(ctss, countScoredOverlaps, features=features,
# 									mc.preschedule=FALSE, mc.cores=cores)
# 	}
#
# 	# Merge into matrix
# 	message("Building Expression Matrix")
# 	x <- do.call(cbind, x)
# 	rownames(x) <- names(features)
#
# 	# Return
# 	x
# }
#
# #' Internal function: Low memory CTSS overlap
# #'
# #' Quickly read in a CTSS file prior to overlap with countScoredOverlaps
# #'
# #' @param fname character: CTSS-BED file
# #' @param features GRanges: Features to be quantified.
# #'
# #' @return numeric vector with the same length as gr
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import S4Vectors IRanges GenomicRanges
# #' @export
# lowMemOverlap <- function(fname, features){
# 	# Load file
# 	d <- fastReadBED(fname)
#
# 	# Turn into GRanges
# 	d <- GRanges(seqnames=Rle(d$seqnames),
# 							 ranges=IRanges(d$start+1, width=1),
# 							 strand=Rle(d$strand),
# 							 score=Rle(d$score))
#
# 	# Overlap
# 	d <- countScoredOverlaps(gr=d, features=features)
#
# 	# Return
# 	d
# }
#
# #' @describeIn quantifyFeatures low memory method
# #' @export
# quantifyFeatures.character <- function(ctss, features, cores=NULL){
# 	# Get some basic info for printing
# 	nFeatures <- length(features)
# 	nSamples <- length(ctss)
# 	m <- sprintf("Quantifying expression of %d features in %d samples",
# 							 nFeatures, nSamples)
# 	message(m)
#
# 	# Choose between serial an paralel
# 	if(is.null(cores)){
# 		x <- pbapply::pblapply(ctss, lowMemOverlap, features=features)
# 	}else{
# 		message(paste0("Using cores: ", cores))
# 		x <- parallel::mclapply(ctss, lowMemOverlap, features=features,
# 														mc.preschedule=FALSE, mc.cores=cores)
# 	}
#
# 	# Merge into matrix
# 	x <- do.call(cbind, x)
# 	rownames(x) <- names(features)
#
# 	# Return
# 	x
# }
