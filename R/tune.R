#### Main S4 functions ####

#' Determine the optimal pooled threshold for unidirectional tag clustering.
#'
#' This function counts the number of Tag Clusters (TCs) for an series of small
#' incremental pooled cutoffs
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Pooled CTSS.
#' @param steps integer: Number of thresholds to analyze (in addition to
#'   treshold=0).
#' @param mergeDist integer: Merge TCs within this distance.
#' @param searchMethod character: For advanced user only, see details.
#' @param maxExponent numeric: The maximal threshold to analyse is obtained as
#'   min(score)*2^maxExponent (only used if searchMethod='exponential').
#' @param ... additional arguments passed to methods.
#'
#' @return data.frame with two columns: threshold and nTCs (number of Tag
#'   Clusters)
#' @family Clustering functions
#' @export
#' @examples
#' data(exampleCTSSs)
#'
#' # Calculate pooledTPM, using supplied number of total tags
#' exampleCTSSs <- calcTPM(exampleCTSSs,
#'                         inputAssay='counts',
#'                         outputAssay='TPM',
#'                         totalTags='totalTags')
#' exampleCTSSs <- calcPooled(exampleCTSSs, inputAssay='TPM')
#'
#' # Set backend
#' library(BiocParallel)
#' register(SerialParam())
#'
#' # Find optimal slice-threshold for reduce distance of 20:
#' tuneTagClustering(object=exampleCTSSs)
setGeneric("tuneTagClustering", function(object, ...) {
    standardGeneric("tuneTagClustering")
})

#' @rdname tuneTagClustering
setMethod("tuneTagClustering", signature(object = "GRanges"),
          function(object, steps = 10L, mergeDist = 20L,
                   searchMethod = "minUnique", maxExponent = 1) {
    # Pre-checks
    assert_that(checkPooled(object),
                is.count(steps),
                is.numeric(maxExponent),
                is.count(mergeDist))

    # Setup series
    message("Finding thresholds to be tested...")
    if (searchMethod == "exponential") {
        exp_series <- min(score(object)) * 2^seq(from = 0,
                                                 to = maxExponent,
                                                 length.out = steps)
    } else if (searchMethod == "minUnique") {
        steps_series <- seq_len(steps)
        exp_series <- sort(unique(score(object)),
                           partial = steps_series)[steps_series]
        rm(steps_series)
    } else {
        stop("searchMethod must be either exponential or minUnique")
    }

    # Add zero
    exp_series <- c(0, exp_series)

    # Split coverage by strand
    message("Splitting by strand...")
    covByStrand <- splitPooled(object)

    # by_strand <- splitByStrand(object)
    # coverage_plus <- coverage(by_strand$`+`, weight = "score")
    # coverage_minus <- coverage(by_strand$`-`, weight = "score")
    # rm(by_strand)

    # Print some info
    message("Iterating over ", length(exp_series),
            " thresholds using ", BiocParallel::bpworkers(), " worker(s)...")

    # Count
    message("Analyzing plus strand...")
    n_plus <- BiocParallel::bpvec(exp_series, countClusters,
                                  cv = covByStrand$`+`, mergeDist = mergeDist)
    message("Analyzing minus strand...")
    n_minus <- BiocParallel::bpvec(exp_series, countClusters,
                                   cv = covByStrand$`-`, mergeDist = mergeDist)

    # Assemble results
    message("Preparing output...")
    o <- data.frame(threshold = exp_series, nTCs = n_plus + n_minus)

    # Return
    o
})

#' @rdname tuneTagClustering
setMethod("tuneTagClustering", signature(object = "RangedSummarizedExperiment"),
    function(object, ...) {
        tuneTagClustering(rowRanges(object), ...)
    })

#' @rdname tuneTagClustering
setMethod("tuneTagClustering",
          signature(object = "GPos"), function(object, ...) {
    tuneTagClustering(methods::as(object, "GRanges"), ...)
})

#### Helpers ####

countClusters <- function(thresholds, cv, mergeDist = 20L) {
    # Slice
    o <- lapply(thresholds, slice,
                x = cv, includeLower = FALSE,
                upper = Inf, rangesOnly = TRUE)

    # Reduce
    o <- lapply(o, reduce, min.gapwidth = mergeDist)

    # Count
    o <- lapply(o, elementNROWS)

    # Sum
    o <- vapply(o, sum, integer(1))

    # Return
    o
}
