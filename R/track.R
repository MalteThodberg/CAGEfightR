#' Create Genome Browser track of CTSSs.
#'
#' Create a Gviz-track of CTSSs, where Plus/minus strand signal is shown
#' positive/negative. This representation makes it easy to identify
#' bidirectional peaks.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Ranges with CTSSs
#'   in the score column.
#' @param plusColor character: Color for plus-strand coverage.
#' @param minusColor character: Color for minus-strand coverage.
#' @param ... additional arguments passed on to DataTrack.
#'
#' @return DataTrack-object.
#' @family Genome Browser functions
#' @export
#' @examples
#' library(Gviz)
#' data(exampleCTSSs)
#' data(exampleUnidirectional)
#' data(exampleBidirectional)
#'
#' # Example uni- and bidirectional clusters
#' TC <- rowRanges(subset(exampleUnidirectional, width>=100)[3,])
#' BC <- rowRanges(exampleBidirectional[3,])
#'
#' # Create pooled trakc
#' subsetOfCTSSs <- subsetByOverlaps(exampleCTSSs, c(BC, TC, ignore.mcols=TRUE))
#' pooledTrack <- trackCTSS(subsetOfCTSSs)
#'
#' # Plot
#' plotTracks(pooledTrack, from=start(TC)-100, to=end(TC)+100,
#'            chromosome=seqnames(TC), name='TC')
#' plotTracks(pooledTrack, from=start(BC)-100, to=end(BC)+100,
#'            chromosome=seqnames(BC), name='BC')
#'
#' # See vignette for examples on how to combine multiple Gviz tracks
setGeneric("trackCTSS", function(object, ...) {
    standardGeneric("trackCTSS")
})

#' @rdname trackCTSS
#' @export
setMethod("trackCTSS", signature(object = "GRanges"), function(object, plusColor = "cornflowerblue",
    minusColor = "tomato", ...) {
    # Pre-checks
    assert_that(isDisjoint(object), !is.null(score(object)), is.numeric(score(object)),
        not_empty(seqlengths(object)), noNA(seqlengths(object)), is.string(plusColor),
        is.string(minusColor))

    # Vector by strand
    message("Splitting pooled signal by strand...")
    by_strand <- splitByStrand(object)
    plus_coverage <- coverage(by_strand$`+`, weight = "score")
    minus_coverage <- 0 - coverage(by_strand$`-`, weight = "score")
    rm(by_strand)

    # Back to GRanges
    message("Preparing track...")
    names(minus_coverage) <- names(plus_coverage)
    o <- bindAsGRanges(plus = plus_coverage, minus = minus_coverage)

    # Build track
    o <- Gviz::DataTrack(o, type = "histogram", groups = c("plus", "minus"), col = c(minusColor,
        plusColor), ...)

    # Return
    o
})

#' @rdname trackCTSS
#' @export
setMethod("trackCTSS", signature(object = "RangedSummarizedExperiment"), function(object,
    ...) {
    trackCTSS(rowRanges(object), ...)
})

#' @rdname trackCTSS
setMethod("trackCTSS", signature(object = "GPos"), function(object, ...) {
    trackCTSS(methods::as(object, "GRanges"), ...)
})

#' Create genome browser track of clusters.
#'
#' Create a Gviz-track of clusters (unidirectional TCs or bidirectional
#' enhancers), where cluster strand and peak is indicated.
#'
#' @param object GRanges: GRanges with peaks in the thick-column.
#' @param plusColor character: Color for plus-strand features.
#' @param minusColor character: Color for minus-strand features.
#' @param unstrandedColor character: Color for unstranded features.
#' @param ... additional arguments passed on to GeneRegionTrack.
#'
#' @return GeneRegionTrack-object.
#'
#' @family Genome Browser functions
#' @export
#' @examples
#' library(Gviz)
#' data(exampleUnidirectional)
#'
#' # Find some wide unidirectional clusters:
#' TCs <- subset(exampleUnidirectional, width >= 100)
#'
#' # Create track
#' clusters_track <- trackClusters(TCs[1:2,], name='Tag clusters', col=NULL)
#'
#' # Plot
#' plotTracks(clusters_track)
#'
#' # See vignette for examples on how to combine multiple Gviz tracks
setGeneric("trackClusters", function(object, ...) {
    standardGeneric("trackClusters")
})

#' @rdname trackClusters
#' @export
setMethod("trackClusters", signature(object = "GRanges"), function(object,
    plusColor = "cornflowerblue", minusColor = "tomato", unstrandedColor = "hotpink",
    ...) {
    # Pre-checks
    assert_that("thick" %in% colnames(mcols(object)), methods::is(mcols(object)[,
        "thick"], "IRanges"), all(poverlaps(mcols(object)$thick, ranges(object),
        type = "within")), is.string(plusColor), is.string(minusColor), is.string(unstrandedColor))

    # Extract peaks
    message("Setting thick and thin features...")
    insideThick <- swapRanges(object)

    # Remove mcols and add features for thin feature
    names(insideThick) <- NULL
    mcols(insideThick) <- NULL
    insideThick$feature <- ifelse(strand(insideThick) == "+", "thickPlus", "thickMinus")
    insideThick$feature <- ifelse(strand(insideThick) == "*", "thickUnstranded",
        insideThick$feature)

    # Remove peaks from TCs
    outsideThick <- setdiff(object, insideThick)

    # Remove mcols and add features
    mcols(outsideThick) <- NULL
    outsideThick$feature <- ifelse(strand(outsideThick) == "+", "thinPlus", "thinMinus")
    outsideThick$feature <- ifelse(strand(outsideThick) == "*", "thinUnstranded",
        outsideThick$feature)

    # Temporary to GRangesList for easy sorting
    message("Merging and sorting...")
    o <- sort(c(insideThick, outsideThick))
    fo <- findOverlaps(o, object, select = "arbitrary")
    o <- split(o, fo)
    names(o) <- names(object)
    o <- unlist(o)
    rm(object)


    # Add necessary columns for track
    o$transcript <- names(o)
    o$gene <- o$transcript
    o$symbol <- o$transcript

    # Build track
    message("Preparing track...")
    o <- Gviz::GeneRegionTrack(o, thinBoxFeature = c("thinPlus", "thinMinus", "thinUnstranded"),
        min.distance = 0, collapse = FALSE, thinPlus = plusColor, thickPlus = plusColor,
        thinMinus = minusColor, thickMinus = minusColor, thinUnstranded = unstrandedColor,
        thickUnstranded = unstrandedColor, ...)

    # Return
    o
})

#' @rdname trackClusters
#' @export
setMethod("trackClusters", signature(object = "RangedSummarizedExperiment"), function(object,
    ...) {
    trackClusters(rowRanges(object), ...)
})

#' Create Genome Browser Track of bidirectional balance scores
#'
#' Visualize balance scores used for detectiong of bidirectional sites. Mainly
#' intended as diagnostic tools for expert user.
#'
#' @param object GenomicRanges or RangedSummarizedExperiment: Ranges with CTSSs
#'   in the score column.
#' @param window integer: Width of sliding window used for calculating windowed
#'   sums.
#' @param plusColor character: Color for plus-strand coverage.
#' @param minusColor character: Color for minus-strand coverage.
#' @param balanceColor character: Color for bidirectional balance.
#' @param ... additional arguments passed to DataTrack.
#'
#' @note Potentially consumes a large amount of memory!
#' @return list of 3 DataTracks for upstream, downstream and balance.
#' @family Genome Browser functions
#' @export
#' @examples
#' library(Gviz)
#' data(exampleCTSSs)
#' data(exampleBidirectional)
#'
#' # Calculate pooled CTSSs
#' exampleCTSSs <- calcTPM(exampleCTSSs, totalTags='totalTags')
#' exampleCTSSs <- calcPooled(exampleCTSSs)
#'
#' # Find a bidirectional cluster to plot:
#' BC <- rowRanges(exampleBidirectional[10,])
#' start(BC) <- start(BC) - 250
#' end(BC) <- end(BC) + 250
#' subsetOfCTSSs <- subsetByOverlaps(exampleCTSSs, BC)
#'
#' # Build balance track
#' balance_track <- trackBalance(subsetOfCTSSs)
#'
#' # Plot
#' plotTracks(balance_track, from=start(BC), to=end(BC),
#'            chromosome=seqnames(BC))
setGeneric("trackBalance", function(object, ...) {
    standardGeneric("trackBalance")
})

#' @rdname trackBalance
#' @export
setMethod("trackBalance", signature(object = "GRanges"), function(object, window = 199,
    plusColor = "cornflowerblue", minusColor = "tomato", balanceColor = "forestgreen",
    ...) {
    # Pre-checks
    assert_that(isDisjoint(object), !is.null(score(object)), is.numeric(score(object)),
        not_empty(seqlengths(object)), noNA(seqlengths(object)), is.string(plusColor),
        is.string(minusColor), is.string(balanceColor))

    # Get windows
    cw <- CAGEfightR:::coverageWindows(pooled = object, window = window, balanceFun = BC)

    # Assemble tracks
    message("Building tracks...")
    o <- list(downstream = Gviz::DataTrack(bindAsGRanges(plus = cw$PD, minus = cw$MD),
        name = "Downstream", type = "l", groups = c("plus", "minus"), col = c(minusColor,
            plusColor)), upstream = Gviz::DataTrack(bindAsGRanges(plus = cw$PU, minus = cw$MU),
        name = "Upstream", type = "l", groups = c("plus", "minus"), col = c(minusColor,
            plusColor)))

    if (!is.null(cw$B)) {
        o$balance <- Gviz::DataTrack(GRanges(cw$B), name = "Balance", type = "l",
            col = balanceColor)
    }

    # Return
    o
})

#' @rdname trackBalance
#' @export
setMethod("trackBalance", signature(object = "GPos"), function(object, ...) {
    trackBalance(methods::as(object, "GRanges"), ...)
})

#' @rdname trackBalance
#' @export
setMethod("trackBalance", signature(object = "RangedSummarizedExperiment"), function(object,
    ...) {
    trackBalance(rowRanges(object), ...)
})
