#' Identify and quantify Transcription Start Sites (TSSs).
#'
#' A convienient wrapper around calcTPM, calcPooled, tuneTagClustering,
#' clusterUnidirectionally and quantifyClusters.
#'
#' @param object RangedSummarizedExperiment: Location and counts of CTSSs,
#'   usually found by calling quantifyCTSSs.
#'
#' @return RangedSummarizedExperiment containing location and counts of TSSs
#' @family Wrapper functions
#' @export
#' @examples
#' # See the CAGEfightR vignette for an overview!
quickTSSs <- function(object) {
    # Pre-checks
    assert_that(methods::is(object, "RangedSummarizedExperiment"), isDisjoint(object),
        not_empty(seqlengths(object)), noNA(seqlengths(object)), is.element("counts",
            assayNames(object)))

    if (is.null(score(rowRanges(object)))) {
        message(" - Running calcTPM and calcPooled:")
        object <- calcTPM(object)
        invisible(gc())
        object <- calcPooled(object)
        invisible(gc())
    } else {
        message("Using existing score column!")
        assert_that(is.numeric(score(rowRanges(object))))
    }

    message("\n - Running tuneTagClustering:")
    tuned <- tuneTagClustering(object, searchMethod = "exponential")
    invisible(gc())
    pooledCutoff <- tuned[which.max(tuned$nTCs), 1]
    message("Optimal pooled cutoff: ", pooledCutoff)

    message("\n - Running clusterUnidirectionally:")
    TCs <- clusterUnidirectionally(object, pooledCutoff = pooledCutoff)
    invisible(gc())

    message("\n - Running quantifyClusters:")
    TCs <- quantifyClusters(object, TCs)
    invisible(gc())

    # Return
    TCs
}

#' Identify and quantify enhancers.
#'
#' A convienient wrapper around clusterBidirectionally, subsetByBidirectionality
#' and quantifyClusters.
#'
#' @param object RangedSummarizedExperiment: Location and counts of CTSSs,
#'   usually found by calling quantifyCTSSs.
#'
#' @return RangedSummarizedExperiment containing location and counts of
#'   enhancers.
#' @family Wrapper functions
#' @export
#' @examples
#' # See the CAGEfightR vignette for an overview!
quickEnhancers <- function(object) {
    # Pre-checks
    assert_that(methods::is(object, "RangedSummarizedExperiment"),
                isDisjoint(object),
                not_empty(seqlengths(object)),
                noNA(seqlengths(object)),
                is.element("counts", assayNames(object)))

    if (is.null(score(rowRanges(object)))) {
        message(" - Running calcTPM and calcPooled:")
        object <- calcTPM(object)
        invisible(gc())
        object <- calcPooled(object)
        invisible(gc())
    } else {
        message("Using existing score column!")
        assert_that(is.numeric(score(rowRanges(object))))
    }

    message("\n - Running clusterBidirectionally:")
    enhancers <- clusterBidirectionally(object)
    invisible(gc())

    message("\n - Running subsetByBidirectionality:")
    enhancers <- subsetByBidirectionality(enhancers, object)
    invisible(gc())

    message("\n - Running quantifyClusters:")
    enhancers <- quantifyClusters(object, enhancers)
    invisible(gc())

    # Return
    enhancers
}
