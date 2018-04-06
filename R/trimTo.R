
symmetricPercentiles <- function(r, prop) {
    # Convert to normal vector
    r <- as.vector(r)
    
    # Calculate proportions
    r_sum <- sum(r)
    r_prop <- r/r_sum
    
    # Cumulate sum from either side
    f_cum <- cumsum(r_prop)
    r_cum <- rev(cumsum(rev(r_prop)))
    
    # Trim both sides
    half_prop <- prop/2
    left <- Position(function(x) x >= half_prop, f_cum, right = FALSE, nomatch = 1)
    right <- Position(function(x) x >= half_prop, r_cum, right = TRUE, nomatch = length(r))
    
    # Return range
    c(left, right)
}

asymmetricPercentiles <- function(r, prop) {
    # Convert to normal vector
    r <- as.vector(r)
    
    # Calculate proportions
    r_sum <- sum(r)
    r_prop <- r/r_sum
    
    # Sort be minimum cumsum from either side
    f_cum <- cumsum(r_prop)
    r_cum <- rev(cumsum(rev(r_prop)))
    min_cum <- pmin(f_cum, r_cum)
    o_cum <- order(min_cum)
    
    # Order original proportions and cut
    r_ord <- cumsum(r_prop[o_cum])
    left <- Position(function(x) x > prop, r_ord, right = FALSE, nomatch = 1)
    o_out <- o_cum[left:length(r_ord)]
    
    # Return range
    range(o_out)
}

#' Trim width of TCs to expression percentiles
#'
#' Given a set of TCs and genome-wide CTSS coverage, reduce the width of TC
#' until a certain amount of expression has been removed.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: TCs to be trimmed.
#' @param pooled GenomicRanges or RangedSummarizedExperiment: CTSS coverage.
#' @param percentile numeric: Fraction of expression to remove from TCs.
#' @param symmetric logical: Whether to trim the same amount from both edges of
#'   the TC (TRUE) or always trim from the least expressed end (FALSE).
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with trimmed TCs, including recalculated peaks and scores.
#' @family Clustering functions
#' @family Trimming functions
#' @export
#' @examples
#' data(exampleCTSSs)
#' data(exampleBidirectional)
#'
#' # Calculate pooled CTSSs
#' exampleCTSSs <- calcTPM(exampleCTSSs, totalTags='totalTags')
#' exampleCTSSs <- calcPooled(exampleCTSSs)
#'
#' # Choose a few wide clusters:
#' TCs <- subset(exampleUnidirectional, width >= 100)
#'
#' # Symmetric trimming (same percentage from each side):
#' TCs_sym <- trimToPercentiles(TCs, pooled=exampleCTSSs, symmetric=FALSE)
#'
#' # Assymmetric trimming (always trim from lowest side):
#' TCs_asym <- trimToPercentiles(TCs, pooled=exampleCTSSs, symmetric=TRUE)
#'
#' # Compare the two results sets of widths:
#' summary(width(TCs_sym) - width(TCs_asym))
setGeneric("trimToPercentiles", function(object, pooled, ...) {
    standardGeneric("trimToPercentiles")
})

#' @rdname trimToPercentiles
setMethod("trimToPercentiles", signature(object = "GRanges", pooled = "GenomicRanges"), 
    function(object, pooled, percentile = 0.1, symmetric = FALSE) {
        # Pre-checks
        assert_that(!is.null(score(pooled)), is.numeric(score(pooled)), is.numeric(percentile), 
            percentile > 0 & percentile <= 1, is.logical(symmetric), identical(seqlengths(object), 
                seqlengths(pooled)))
        
        # Split by strand
        message("Splitting by strand...")
        coverage_stranded <- splitByStrand(pooled)
        coverage_plus <- coverage(coverage_stranded$`+`, weight = "score")
        coverage_minus <- coverage(coverage_stranded$`-`, weight = "score")
        tcs_stranded <- splitByStrand(object)
        rm(coverage_stranded)
        
        # Convert to IRangesList irl_plus <- methods::as(tcs_stranded$`+`,'RangesList')
        # irl_minus <- methods::as(tcs_stranded$`-`,'RangesList')
        irl_plus <- split(ranges(tcs_stranded$`+`), seqnames(tcs_stranded$`+`))
        irl_minus <- split(ranges(tcs_stranded$`-`), seqnames(tcs_stranded$`-`))
        rm(tcs_stranded)
        
        # Obtain views
        views_plus <- Views(coverage_plus, irl_plus)
        views_minus <- Views(coverage_minus, irl_minus)
        
        # Extract Adjustments
        if (symmetric) {
            message(paste0("Symmetric trimming to percentile: ", percentile * 100, 
                "%"))
            adjust_plus <- viewApply(views_plus, symmetricPercentiles, prop = percentile, 
                simplify = FALSE)
            adjust_minus <- viewApply(views_minus, symmetricPercentiles, prop = percentile, 
                simplify = FALSE)
        } else if (!symmetric) {
            message(paste0("Asymmetric trimming to percentile: ", percentile * 100, 
                "%"))
            adjust_plus <- viewApply(views_plus, asymmetricPercentiles, prop = percentile, 
                simplify = FALSE)
            adjust_minus <- viewApply(views_minus, asymmetricPercentiles, prop = percentile, 
                simplify = FALSE)
        } else {
            stop("Additional percentile functions not yet implemented!")
        }
        rm(views_plus, views_minus)
        
        # Adjustments as IntegerLists
        message("Adjusting ranges...")
        left_plus <- methods::as(lapply(adjust_plus, function(x) lapply(x, function(x) x[1])), 
            "IntegerList")
        right_plus <- methods::as(lapply(adjust_plus, function(x) lapply(x, function(x) x[2])), 
            "IntegerList")
        left_minus <- methods::as(lapply(adjust_minus, function(x) lapply(x, function(x) x[1])), 
            "IntegerList")
        right_minus <- methods::as(lapply(adjust_minus, function(x) lapply(x, function(x) x[2])), 
            "IntegerList")
        
        # Narrow ranges
        irl_plus <- narrow(irl_plus, start = left_plus, end = right_plus)
        irl_minus <- narrow(irl_minus, start = left_minus, end = right_minus)
        rm(left_plus, right_plus, left_minus, right_minus)
        
        # Calculate new stats
        message("Calculating new stats...")
        trimmedTCs <- TCstats(coverage_plus = coverage_plus, coverage_minus = coverage_minus, 
            tcs_plus = irl_plus, tcs_minus = irl_minus)
        rm(coverage_plus, coverage_minus, irl_plus, irl_minus)
        
        # Carry over seqinfo and sort
        message("Preparing output...")
        seqinfo(trimmedTCs) <- seqinfo(object)
        trimmedTCs <- sort(trimmedTCs)
        
        # Print some basic stats
        message("Tag clustering summary:")
        summarizeWidths(trimmedTCs)
        
        # Return
        trimmedTCs
    })

#' @rdname trimToPercentiles
setMethod("trimToPercentiles", signature(object = "RangedSummarizedExperiment", pooled = "GenomicRanges"), 
    function(object, pooled, ...) {
        trimToPercentiles(rowRanges(object), pooled, ...)
    })

#' @rdname trimToPercentiles
setMethod("trimToPercentiles", signature(object = "GRanges", pooled = "RangedSummarizedExperiment"), 
    function(object, pooled, ...) {
        trimToPercentiles(object, rowRanges(pooled), ...)
    })

#' @rdname trimToPercentiles
setMethod("trimToPercentiles", signature(object = "RangedSummarizedExperiment", pooled = "RangedSummarizedExperiment"), 
    function(object, pooled, ...) {
        trimToPercentiles(rowRanges(object), rowRanges(pooled), ...)
    })

#' Trim width of TCs by distance from TC peak
#'
#' Trim the width of TCs by distance from the TC peaks.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Tag clusters.
#' @param pooled GenomicRanges or RangedSummarizedExperiment: Basepair-wise
#'   pooled CTSS (stored in the score column).
#' @param upstream integer: Maximum upstream distance from TC peak.
#' @param downstream integer: Maximum downstream distance from TC peak.
#' @param peaks character: Name of column in TCs holding TC peaks as an IRanges.
#' @param ... additional arguments passed to methods.
#'
#' @return data.frame with two columns: threshold and nTCs (number of Tag
#'   Clusters)
#' @family Clustering functions
#' @family Trimming functions
#' @export
#' @examples
#' data(exampleCTSSs)
#' data(exampleBidirectional)
#'
#' # Calculate pooled CTSSs
#' exampleCTSSs <- calcTPM(exampleCTSSs, totalTags='totalTags')
#' exampleCTSSs <- calcPooled(exampleCTSSs)
#'
#' # Choose a few wide clusters:
#' TCs <- subset(exampleUnidirectional, width >= 100)
#'
#' # Trim to +/- 10 bp of TC peak
#' trimToPeak(TCs, pooled=exampleCTSSs, upstream=10, downstream=10)
setGeneric("trimToPeak", function(object, pooled, ...) {
    standardGeneric("trimToPeak")
})

#' @rdname trimToPeak
setMethod("trimToPeak", signature(object = "GRanges", pooled = "GenomicRanges"), 
    function(object, pooled, upstream, downstream, peaks = "thick") {
        # Pre-Checks
        assert_that(is.count(upstream), is.count(downstream), identical(seqlengths(object), 
            seqlengths(pooled)))
        
        # Extract peaks
        message("Trimming TCs around peaks...")
        peaks <- swapRanges(object)
        
        # Expand peaks
        expanded <- promoters(peaks, upstream = upstream, downstream = downstream)
        
        # Intersect to obtain filtered peaks.
        ranges(peaks) <- pintersect(ranges(object), ranges(expanded))
        
        # Split by strand
        message("Splitting by strand...")
        coverage_stranded <- splitByStrand(pooled)
        coverage_plus <- coverage(coverage_stranded$`+`, weight = "score")
        coverage_minus <- coverage(coverage_stranded$`-`, weight = "score")
        tcs_stranded <- splitByStrand(peaks)
        rm(coverage_stranded, peaks)
        
        # Convert to IRangesList irl_plus <- methods::as(tcs_stranded$`+`,'RangesList')
        # irl_minus <- methods::as(tcs_stranded$`-`,'RangesList')
        irl_plus <- split(ranges(tcs_stranded$`+`), seqnames(tcs_stranded$`+`))
        irl_minus <- split(ranges(tcs_stranded$`-`), seqnames(tcs_stranded$`-`))
        rm(tcs_stranded)
        
        # Calculate new stats
        message("Calculating new stats...")
        trimmedTCs <- TCstats(coverage_plus = coverage_plus, coverage_minus = coverage_minus, 
            tcs_plus = irl_plus, tcs_minus = irl_minus)
        rm(coverage_plus, coverage_minus, irl_plus, irl_minus)
        
        # Carry over seqinfo and sort
        message("Preparing output...")
        seqinfo(trimmedTCs) <- seqinfo(object)
        trimmedTCs <- sort(trimmedTCs)
        
        # Print some basic stats
        message("Tag clustering summary:")
        summarizeWidths(trimmedTCs)
        
        # Return
        trimmedTCs
    })

#' @rdname trimToPeak
setMethod("trimToPeak", signature(object = "RangedSummarizedExperiment", pooled = "GenomicRanges"), 
    function(object, pooled, ...) {
        trimToPeak(rowRanges(object), pooled, ...)
    })

#' @rdname trimToPeak
setMethod("trimToPeak", signature(object = "GRanges", pooled = "RangedSummarizedExperiment"), 
    function(object, pooled, ...) {
        trimToPeak(object, rowRanges(pooled), ...)
    })

#' @rdname trimToPeak
setMethod("trimToPeak", signature(object = "RangedSummarizedExperiment", pooled = "RangedSummarizedExperiment"), 
    function(object, pooled, ...) {
        trimToPeak(rowRanges(object), rowRanges(pooled), ...)
    })

