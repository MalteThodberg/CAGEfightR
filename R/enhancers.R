#### Helper functions ####

padlag <- function(r, k=0){
	stopifnot(class(r) == "Rle")
	stopifnot(k <= length(r))

	if(k == 0){
		o <- r
	}else if(k > 0){
		o <- c(window(r, 1+k), Rle(values=NaN, lengths=abs(k)))
	}else if(k < 0){
		o <- c(Rle(values=NaN, lengths=abs(k)), window(r, 1, length(r)-abs(k)))
	}

	stopifnot(length(r) == length(o))
	o
}

coverageWindows <- function(pooled, window, balanceFun, ...){
	# Obtain shift
	shift_val <- ceiling(window / 2)

	# Split into strand
	covByStrand <- splitByStrand(pooled)
	rm(pooled)

	message("Calculating windowed coverage on plus strand...")
	P <- coverage(covByStrand$`+`, weight="score")
	P <- runsum(P, k=window, endrule="constant")
	P[P < 0] <- 0
	PD <- endoapply(P, padlag, k=shift_val)
	PU <- endoapply(P, padlag, k=-shift_val)
	rm(P)

	message("Calculating windowed coverage on minus strand...")
	M <- coverage(covByStrand$`-`, weight="score")
	M <- runsum(M, k=window, endrule="constant")
	M[M < 0] <- 0
	MD <- endoapply(M, padlag, k=-shift_val)
	MU <- endoapply(M, padlag, k=shift_val)
	rm(covByStrand, M)

	if(!is.null(balanceFun)){
		message("Calculating balance score...")
		B <- mendoapply(BC, PU, PD, MU, MD)
	}else{
		B <- NULL
	}

	# Build output
	o <- list(PD=PD, MD=MD, PU=PU, MU=MU, B=B)

	# Return
	o
}

BC <- function(plusUpstream, plusDownstream, minusUpstream, minusDownstream){
	# Check all input have the same class
	input_classes <- c(class(plusUpstream), class(plusDownstream),
										 class(minusUpstream), class(minusDownstream))
	single_class <- unique(input_classes)
	stopifnot(length(single_class) == 1)

	# Sum of all
	S <- plusUpstream + plusDownstream + minusUpstream + minusDownstream

	# Only calcuate downstream arms - rest will be zero
	B <- sqrt((minusDownstream / S) * 0.5) + sqrt((plusDownstream / S) * 0.5)

	# Checks
	stopifnot(class(B) == single_class)

	# Return
	B
}

#### Main functions ####

#' Bidirectional clustering of pooled CTSSs.
#'
#' Finds sites with (balanced and divergent) bidirectional transcription using
#' sliding windows of summed coverage: The Bhattacharyya coefficient (BC) is
#' used to quantify depature from a perfectly balanced site, and a slice-reduce
#' is used to identify sites.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Pooled CTSSs
#'   stored in the score column.
#' @param window integer: Width of sliding window used for calculating window
#'   sums.
#' @param balanceThreshold numeric: Minimum value of the BC to use for
#'   slice-reduce, a value of 1 corresponds to perfectly balanced sites.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with bidirectional sites: Minimum width is 1 + 2*window, TPM
#'   sum (on both strands) in the score column, maximal bidirectional site in
#'   the thick column and maximum balance in the balance column.
#' @family Clustering functions
#'
#' @import assertthat S4Vectors IRanges GenomicRanges
#' @export
#' @examples
#' \dontrun{
#' data(exampleCTSSs)
#'
#' # Calculate pooledTPM, using supplied number of total tags
#' exampleCTSSs <- calcTPM(exampleCTSSs,
#'                         inputAssay="counts",
#'                         outputAssay="TPM",
#'                         totalTags="totalTags")
#' exampleCTSSs <- calcPooled(exampleCTSSs, inputAssay="TPM")
#'
#' # Cluster using defaults: balance-treshold of 199 and window of 199 bp:
#' clusterBidirectionally(exampleCTSSs)
#'
#' # Use custom thresholds:
#' clusterBidirectionally(exampleCTSSs, balanceThreshold=0.99, window=101)
#' }
setGeneric("clusterBidirectionally", function(object, ...) {
	standardGeneric("clusterBidirectionally")
})

#' @rdname clusterBidirectionally
setMethod("clusterBidirectionally", signature(object="GenomicRanges"),
					function(object, window=199, balanceThreshold=0.95){
	# Pre-checks
	assert_that(isDisjoint(object),
							noNA(seqlengths(object)),
							!is.null(score(object)),
							is.numeric(score(object)),
							is.count(window),
							window %% 2 == 1,
							window >= 3,
							is.number(balanceThreshold),
							balanceThreshold >= 0 & balanceThreshold <= 1)

	# Get windows
	cw <- coverageWindows(pooled=object, window=window, balanceFun=BC)
	B <- cw$B
	rm(cw)

	# Prepare for slicing
	if(any(any(B < 0, na.rm=TRUE))){
		stop("Balance function produced values below 0!",
				 "Output values must be in range [0-Inf) to allow for slicing...")
	}else{
		B[is.na(B)] <- 0
	}

	message("Slice-reduce to find bidirectional clusters...")
	# Slice
	bidirLoci <- slice(B, lower=balanceThreshold, upper=Inf,
										 includeLower=FALSE, rangesOnly=TRUE)

	# Merge
	mergeDist <- window * 2
	bidirLoci <- reduce(bidirLoci, min.gapwidth=mergeDist)
	stopifnot(isDisjoint(bidirLoci))

	# Expand by window size
	start(bidirLoci) <- methods::as(start(bidirLoci) - window, "IntegerList")
	end(bidirLoci) <- methods::as(end(bidirLoci) + window, "IntegerList")
	stopifnot(isDisjoint(bidirLoci))

	if(!sum(elementNROWS(bidirLoci), na.rm=TRUE) > 0){
		warning("No bidirectional sites found at the given balance threshold!")
	}

	message("Calculating statistics...")
	# Coverage on both strands
	A <- coverage(object, weight="score")

	# Views on score and balance
	viewsA <- Views(A, bidirLoci)
	viewsB <- Views(B, bidirLoci)
	rm(A, B)

	#  Scores and midpoint
	scores <- viewSums(viewsA)
	peaks <-  viewRangeMaxs(viewsB)
	peaks <- resize(x=peaks, width=1, fix="center", use.names=FALSE)

	# Maximum balance
	balance <- viewMaxs(viewsB)

	# Assemble output
	o <- GRanges(bidirLoci,
							 score=unlist(scores),
							 thick=unlist(peaks),
							 balance=unlist(balance))
	rm(scores, peaks, balance, viewsA, viewsB)

	# Carry over seqinfo
	message("Preparing output...")
	seqinfo(o) <- seqinfo(object)
	o <- sort(o)

	# Names as IDs for both ranges and peaks
	o_ids <- paste0(seqnames(o), ":", start(o), "-", end(o))
	names(o) <- o_ids
	names(o$thick) <- o_ids
	rm(o_ids)

	# Return
	o
})

#' @rdname clusterBidirectionally
setMethod("clusterBidirectionally",
					signature(object="RangedSummarizedExperiment"),
					function(object, ...){
	clusterBidirectionally(rowRanges(object), ...)
})

#' Calculate sample-wise bidirectionally of clusters.
#'
#' For each cluster, calculate how many individual samples shows transcription
#' in both directions. This is refered to as the "bidirectionality". Clusters
#' must be unstranded (*) and have a midpoint stored in the thick column
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Unstranded
#'   clusters with midpoints stored in the "thick" column.
#' @param samples RangedSummarizedExperiment: Sample-wise CTSSs stored as an
#'   assay.
#' @param inputAssay character: Name of assay in samples holding input CTSS
#'   values.
#' @param outputColumn character: Name of column in object to hold
#'   bidirectionality values.
#' @param ... additional arguments passed to methods.
#'
#' @return object returned with bidirectionality scores added in rowData (or
#'   mcols).
#' @family Calculation functions
#' @export
#' @examples
#' data(exampleCTSSs)
#' data(exampleBidirectional)
#'
#' calcBidirectionality(exampleBidirectional, samples=exampleCTSSs)
setGeneric("calcBidirectionality", function(object, ...) {
	standardGeneric("calcBidirectionality")
})

#' @rdname calcBidirectionality
#' @export
setMethod("calcBidirectionality", signature(object="GRanges"),
					function(object, samples, inputAssay="counts",
									 outputColumn="bidirectionality"){
	# Pre-checks
	assert_that("thick" %in% colnames(mcols(object)),
							methods::is(mcols(object)[,"thick"], "IRanges"),
							all(poverlaps(mcols(object)$thick,
														ranges(object),
														type = "within")),
							all(strand(object) == "*"),
							methods::is(samples, "RangedSummarizedExperiment"),
							isDisjoint(samples),
							not_empty(seqlengths(samples)),
							noNA(seqlengths(samples)),
							#!is.null(score(rowRanges(samples))),
							#is.numeric(score(rowRanges(samples))),
							is.string(inputAssay),
							inputAssay %in% assayNames(samples),
							is.string(outputColumn))

	# Warnings
	if(outputColumn %in% colnames(mcols(object))){
		warning("object already has a column named ",
						outputColumn," in mcols: It will be overwritten!")
	}

	# Extract arms
	arms_plus <- object
	start(arms_plus) <- end(arms_plus$thick)
	strand(arms_plus) <- "+"

	arms_minus <- object
	end(arms_minus) <- start(arms_minus$thick)
	strand(arms_minus) <- "-"

	# Quantify
	mat_plus <- suppressMessages(quantifyClusters(samples, arms_plus,
																								inputAssay=inputAssay))
	mat_plus <- assay(mat_plus, inputAssay) > 0
	mat_minus <- suppressMessages(quantifyClusters(samples, arms_minus,
																								 inputAssay=inputAssay))
	mat_minus <- assay(mat_minus, inputAssay) > 0

	# Compare and count
	o <- mat_plus & mat_minus
	rm(mat_plus, mat_minus)
	o <- rowSums(o)

	# Post-checks
	stopifnot(length(o) == length(object),
						max(o) <= ncol(samples))

	# Add to object
	mcols(object)[,outputColumn] <- o

	# Return
	object
})

#' @rdname calcBidirectionality
#' @export
setMethod("calcBidirectionality",
					signature(object="RangedSummarizedExperiment"), function(object, ...){
	rowRanges(object) <- calcBidirectionality(rowRanges(object), ...)
	object
})

#' Subset by sample-wise bidirectionality of clusters.
#'
#' A convenient wrapper around calcBidirectionality and subset.
#'
#' @param object GRanges or RangedSummarizedExperiment: Unstranded clusters with
#'   peaks stored in the "thick" column.
#' @param samples RangedSummarizedExperiment: Sample-wise CTSSs stored as an
#'   assay.
#' @param inputAssay character: Name of assay in samples holding input CTSS
#'   values.
#' @param outputColumn character: Name of column in object to hold
#'   bidirectionality values.
#' @param minSamples integer: Only regions with bidirectionality above this
#'   value are retained.
#' @param ... additional arguments passed to methods.
#'
#' @return object with bidirectionality values added as a column, and low
#'   bidirectionaly regions removed.
#'
#' @family Subsetting functions
#' @family Calculation functions
#'
#' @export
#' @examples
#' data(exampleCTSSs)
#' data(exampleBidirectional)
#'
#' # Keep only clusters that are bidirectional in at least one sample:
#' subsetByBidirectionality(exampleBidirectional, samples=exampleCTSSs)
setGeneric("subsetByBidirectionality", function(object, ...) {
	standardGeneric("subsetByBidirectionality")
})

#' @rdname subsetByBidirectionality
#' @export
setMethod("subsetByBidirectionality", signature(object="GRanges"),
					function(object, samples, inputAssay="counts",
									 outputColumn="bidirectionality", minSamples=0){
	assert_that(is.number(minSamples))

	# Call function
	message("Calculating bidirectionality...")
	object <- calcBidirectionality(object=object,
																 samples=samples,
																 inputAssay=inputAssay,
																 outputColumn=outputColumn)
	before <- length(object)

	# Subset
	message("Subsetting...")
	object <- object[mcols(object)[,outputColumn] > minSamples,]
	after <- length(object)
	removed <- before-after

	# Print some info
	message("Removed ", removed, " out of ", before,
					" regions (", round(removed/before*100, digits=1), "%)")

	# Return
	object
})

#' @rdname subsetByBidirectionality
#' @export
setMethod("subsetByBidirectionality",
					signature(object="RangedSummarizedExperiment"), function(object, ...){
	# Force call of GRanges generic
	methods::selectMethod("subsetByBidirectionality",
												"GRanges")@.Data(object=object, ...)
})
