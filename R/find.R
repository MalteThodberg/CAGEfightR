#### Stretches ####

#' Find stretches of clusters
#'
#' Finds stretches or groups of clusters along the genome, where each cluster is
#' within a certain distance of the next. Once stretches have been identified,
#' the average pairwise correlation between all clusters in the stretch is
#' calculated. A typical use case is to look for stretches of enhancers, often
#' refered to as "super enhancers".
#'
#' @param object GRanges or RangedSummarizedExperiment: Clusters, possibly with
#'   expression for calculating correlations.
#' @param inputAssay character: Name of assay holding expression values (if
#'   object is a RangedSummarizedExperiment)
#' @param mergeDist integer: Maximum distance between clusters to be merged into
#'   stretches.
#' @param minSize integer: Minimum number of clusters in stretches.
#' @param corFun function: Function for calculating correlations. Should behave
#'   and produce output similar to cor().
#' @param ... additional arguments passed to methods or ultimately corFun.
#'
#' @return A GRanges containing stretches with number of clusters and average
#'   pairwise correlations calculated. The revmap can be used to retrieve the
#'   original clusters (see example below.)
#' @family Spatial functions
#' @export
#'
#' @examples
#' # Calculate TPM values for bidirectional clusters
#' data(exampleBidirectional)
#' BCs <- calcTPM(exampleBidirectional)
#'
#' # Find stretches
#' pearson_stretches <- findStretches(BCs, inputAssay="TPM")
#'
#' # Use Kendall instead of pearson and require bigger stretches
#' kendall_stretches <- findStretches(BCs, inputAssay="TPM",
#'                                    minSize=5, method="kendall")
#'
#' # Use the revmap to get stretches as a GRangesList
#' grl <- extractList(rowRanges(BCs), kendall_stretches$revmap)
#' names(grl) <- names(kendall_stretches)
setGeneric("findStretches", function(object, ...) standardGeneric("findStretches"))

#' @rdname findStretches
setMethod("findStretches", signature = "GRanges",
          definition = function(object, mergeDist=10000L, minSize=3L){
              # Pre-checks
              assert_that(is.count(mergeDist),
                          is.count(minSize))

              # Find stretches
              stretches <- GenomicRanges::reduce(object,
                                                 min.gapwidth=mergeDist,
                                                 ignore.strand=TRUE,
                                                 with.revmap=TRUE)

              # Subset to minimum size
              stretches <- stretches[elementNROWS(stretches$revmap) >= minSize]
              names(stretches) <- as.character(stretches)

              # Count number of clusters
              stretches$nClusters <- elementNROWS(stretches$revmap)

              # Return
              stretches
          })

#' @rdname findStretches
setMethod("findStretches", signature = "RangedSummarizedExperiment",
          definition = function(object, inputAssay, mergeDist=10000L, minSize=3L, corFun=cor, ...){
              assert_that(is.string(inputAssay),
                          inputAssay %in% assayNames(object),
                          is.function(corFun))

              message("Finding stretches...")
              stretches <- findStretches(rowRanges(object),
                                         mergeDist=mergeDist,
                                         minSize=minSize)


              message("Calculating correlations...")
              # Extract and transpose matrices
              aveCors <- extractList(x=assay(object, inputAssay),
                                     i=stretches$revmap)
              aveCors <- endoapply(aveCors, t)

              # Calculate pairwise correlations
              aveCors <- endoapply(aveCors, corFun, ...)

              # Calculate mean (lower part of matrix)
              aveCors <- vapply(aveCors,
                                FUN=function(x) mean(x[lower.tri(x)]),
                                FUN.VALUE = numeric(1))

              # Append and remove
              stretches$aveCor <- aveCors
              rm(aveCors)

              message("# Stretch summary:")
              message("Number of stretches: ", length(stretches))
              message("Total number of clusters inside stretches: ",
                      sum(stretches$nClusters), " / ", length(object))
              message("Minimum clusters: ", min(stretches$nClusters))
              message("Maximum clusters: ", max(stretches$nClusters))
              message("Minimum width: ", min(width(stretches)))
              message("Maximum width: ", max(width(stretches)))
              message("Summary of average pairwise correlations: ")
              tmp <- utils::capture.output(summary(stretches$aveCor))
              message(tmp[1])
              message(tmp[2])

              # Return
              stretches
          })

#### Links ####

# Determine orientation
determineOrientation <- function(gi){
    # Get positions as midpoint
    a <- InteractionSet::anchors(gi)
    a <- lapply(a, granges)
    a <- lapply(a, resize, width=1, fix="center")

    # Switch depending on first strand
    o <- ifelse(strand(a$first) == "+",
                start(a$first) >= start(a$second),
                start(a$first) <= start(a$second))
    o <- ifelse(o, "upstream", "downstream")

    # Return
    o
}

# Helper for bpvec
corHelper <- function(ii, m, corFun, ...){
    suppressWarnings(lapply(ii, function(i) corFun(x=m[i[1],],
                                  y=m[i[2],],
                                  ...)))
}


#' Find nearby pairs of clusters and calculate pairwise correlations.
#'
#' Finds all links or pairs of clusters within a certain distance of each other
#' and then calculates the correlation between them. The links found can be
#' restricted to only be between two classes, for example TSSs to enhancers.
#'
#' @param object GRanges or RangedSummarizedExperiment: Clusters, possibly with
#'   expression for calculating correlations.
#' @param inputAssay character: Name of assay holding expression values (if
#'   object is a RangedSummarizedExperiment)
#' @param maxDist integer: Maximum distance between links.
#' @param directional character: Name of a column in object holding a grouping
#'   of the clusters. This must be a factor with two levels. The first level is
#'   used as the basis for calculating orientation (see below).
#' @param corFun function: Function for calculating pairwise correlations. See
#'   notes for supplying custom functions.
#' @param vals character: Statistics extracted from the results produced by
#'   corFun. See notes for supplying custom functions.
#' @param ... additional arguments passed to methods or ultimately corFun.
#'
#' @return A GInteractions holding the links, along with the distance between
#'   them and correlation estimate and p-value calculated from their expression.
#'   If a directional analysis was performed, the two anchors are always
#'   connecting members of the two classes and the orientation of the second
#'   anchor relative to the first is additionaly calculated (e.g. whether an
#'   enhancers is upstream or downstream of the TSS).
#' @family Spatial functions
#'
#' @details A custom function for calculation correlations can be supplied by
#'   the user. The output of this function must be a named list or vector of
#'   numeric values. The names of the vals to be extracted should be supplied to
#'   vals.
#' @export
#'
#' @examples
#' library(InteractionSet)
#'
#' # Subset to highly expressed unidirectional clusters
#' TCs <- subset(exampleUnidirectional, score > 10)
#'
#' # Find links within a certain distance
#' findLinks(TCs, inputAssay="counts", maxDist=10000L)
#'
#' # To find TSS-to-enhancer type links, first merge the clusters:
#' colData(exampleBidirectional) <- colData(TCs)
#' rowRanges(TCs)$clusterType <- "TSS"
#' rowRanges(exampleBidirectional)$clusterType <- "Enhancer"
#' SE <- combineClusters(TCs, exampleBidirectional, removeIfOverlapping="object1")
#' rowRanges(SE)$clusterType <- factor(rowRanges(SE)$clusterType, levels=c("TSS", "Enhancer"))
#'
#' # Calculate kendall correlations of TPM values:
#' SE <- calcTPM(SE, totalTags="totalTags")
#' findLinks(SE, inputAssay="TPM", maxDist=10000L, directional="clusterType", method="kendall")
setGeneric("findLinks", function(object, ...) standardGeneric("findLinks"))

#' @rdname findLinks
setMethod("findLinks", signature = "GRanges",
          definition = function(object, maxDist=10000L, directional=NULL){
              # Pre-checks
              assert_that(is.count(maxDist))

              # Check if supplied factor
              if(!is.null(directional)){
                  assert_that(is.string(directional),
                              directional %in% colnames(mcols(object)),
                              is.factor(mcols(object)[,directional]),
                              length(levels(mcols(object)[,directional])) == 2)

                  # Grouping factor
                  g <- mcols(object)[,directional]
                  message("Finding directional links from ",
                          levels(g)[1], " to ", levels(g)[2], "...")

                  # Remember original order then split
                  by_dir <- granges(object)
                  by_dir$.id <- seq_len(length(by_dir))
                  by_dir <- split(by_dir, g)

                  # Find hits only between pairs
                  hits <- findOverlaps(by_dir[[1]],
                                       by_dir[[2]],
                                       maxgap=maxDist)

                  # Rebuild with original order
                  links <- InteractionSet::GInteractions(
                      anchor1 = by_dir[[1]]$.id[from(hits)],
                      anchor2 = by_dir[[2]]$.id[to(hits)],
                      regions = object)

                  # Determine orientation
                  links$orientation <- determineOrientation(links)
              }else{
                  message("Finding links...")
                  # Find nearby clusters
                  hits <- findOverlaps(object,
                                       drop.self=TRUE,
                                       drop.redundant=TRUE,
                                       maxgap=maxDist)

                  # Coerce to GI
                  links <- InteractionSet::GInteractions(anchor1 = from(hits),
                                         anchor2 = to(hits),
                                         regions = object)
              }

              # Calculate distance
              links$distance <- InteractionSet::pairdist(links, type="gap")

              # Return
              links
          })

#' @rdname findLinks
setMethod("findLinks", signature = "RangedSummarizedExperiment",
          definition = function(object, inputAssay, maxDist=10000L,
                                directional=NULL, corFun=stats::cor.test,
                                vals=c("estimate", "p.value"), ...){
              # Pre-checks
              assert_that(is.string(inputAssay),
                          inputAssay %in% assayNames(object),
                          is.function(corFun),
                          is.character(vals),
                          length(vals) >= 1)

              # Find links
              links <- findLinks(rowRanges(object),
                                 maxDist=maxDist,
                                 directional=directional)

              message("Calculating ", length(links), " pairwise correlations...")
              # Get pair IDs and matrix
              iii <- Pairs(InteractionSet::anchorIds(links, "first"),
                           InteractionSet::anchorIds(links, "second"))
              iii <- zipup(iii)
              m <- assay(object, inputAssay)

              # Calculate
              # linkwise <- lapply(iii, function(i) corFun(x=m[i[1],],
              #                                           y=m[i[2],],
              #                                           ...))
              linkwise <- BiocParallel::bpvec(X=iii,
                                FUN=corHelper,
                                m=m,
                                corFun=corFun,
                                ...)
              rm(iii, m)

              message("Preparing output...")
              # Extract and transpose
              linkwise <- vapply(linkwise,
                                 FUN=function(o) unlist(o[vals], use.names=FALSE),
                                 FUN.VALUE = numeric(length(vals)))

              # Extract output, will be either vector or matrix
              if(length(vals) == 1){
                  if(is.vector(linkwise)){
                      linkwise <- as.matrix(linkwise)
                  }else{
                      stop("Could not properly format output:\n",
                           "Output was neither a vector nor matrix!")
                  }
              }else{
                  linkwise <- t(linkwise)
              }
              colnames(linkwise) <- vals

              # Append
              mcols(links) <- cbind(mcols(links), linkwise)
              rm(linkwise)

              message("# Link summary:")
              message("Number of links: ", length(links))
              message("Summary of pairwise distance: ")
              tmp <- utils::capture.output(summary(links$distance))
              message(tmp[1])
              message(tmp[2])

              # Return
              links
          })



