#### Range utilities ####

#' Utility: Split Genomic Ranges by strand
#'
#' Utility function that attemps to split genomic ranges by strand with split(object, strand(object))
#'
#' @param object Any object with a split and strand method, e.g. GRanges/GPos
#'
#' @return Object split by strand, e.g. GRangesList.
#' @family Utility functions
#' @export
#' @examples
#' gp <- GPos(seqnames=Rle(c("chr1", "chr2", "chr1"), c(10, 6, 4)),
#'             pos=c(44:53, 5:10, 2:5),
#'             strand=c(rep("+", 10), rep("-", 10)))
#' gr <- as(gp, "GRanges")
#' utilsDeStrand(gp)
#' utilsDeStrand(gr)
utilsDeStrand <- function(object) {
    split(object, strand(object))
}

# Kept here for backwards compatability
extendRanges <- function(object, upstream, downstream) {
    U <- flank(object, width = upstream, start = TRUE)
    D <- flank(object, width = downstream, start = FALSE)
    object <- punion(object, U)
    object <- punion(object, D)
    object
}

#### Counting Utilities ####

# Kept here for backwards compatability
countScoredOverlaps <- function(query, subject) {
    # Solution from bioconductor by M Morgan
    hits <- methods::as(findOverlaps(query = query, subject = subject), "List")
    weightedCount <- sum(extractList(score(subject), hits))

    # Return
    weightedCount
}

#' Utility: Counting overlaps taking into account scores
#'
#' Similar to countOverlaps, but takes the score column into account.
#'
#' @param query same as findOverlaps/countOverlaps
#' @param subject same as findOverlaps/countOverlaps
#' @param ... additional arguments passed to findOverlaps
#'
#' @return vector of number of overlaps weigthed by score column.
#' @seealso \url{https://support.bioconductor.org/p/87736/#87758}
#' @family Utility functions
#' @export
#' @examples
#' gr1 <- GRanges(seqnames="chr1",
#'                ranges=IRanges(start = c(4, 9, 10, 30),
#'                               end = c(4, 15, 20, 31)),
#'                strand="+")
#' gr2 <- GRanges(seqnames="chr1",
#'                ranges=IRanges(start = c(1, 4, 15, 25),
#'                               end = c(2, 4, 20, 26)),
#'                strand=c("+"),
#'                score=c(10, 20, 15, 5))
#' countOverlaps(gr1, gr2)
#' utilsScoreOverlaps(gr1, gr2)
utilsScoreOverlaps <- function(query, subject, ...) {
    # Solution from bioconductor by M Morgan
    o <- findOverlaps(query=query,  subject=subject, ...)
    o <- methods::as(o, "List")
    o <- sum(extractList(score(subject), o))

    # Return
    o
}

#' Utility: Aggregate rows
#'
#' Used by quantifyClusters and quantifyGenes. Wrapper around rowsum with a few
#' improvements: 1) Handles dgCMatrix 2) Suppresses warnings from and discards
#' NAs in grouping 3) Checks if output can be coerced to integer (useful when
#' aggregating a dgCMatrix), 4) For the dgCMatrix case, has the option to keep
#' unused levels and output a sparse matrix.
#'
#' @param x matrix or dgCMatrix: Matrix to be aggregated.
#' @param group factor: Grouping, can cannot NAs which will be discarded.
#' @param drop logical: Whether to drop unused levels (TRUE) or keep assign them
#'   0 (FALSE).
#' @param sparse logical: Whether output should be coerced to a dense matrix.
#'
#' @return matrix (or dgCMatrix if sparse=TRUE)
#'
#' @family Utility functions
#' @export
#' @examples
#' library(Matrix)
#' data("exampleCTSSs")
#' data("exampleUnidirectional")
#'
#' # Sparse and dense examples
#' sparse_matrix <- assay(exampleCTSSs)
#' dense_matrix <- as(sparse_matrix, "matrix")
#'
#' # Groupings
#' grp <- findOverlaps(query = exampleCTSSs,
#'                   subject = exampleUnidirectional,
#'                   select="arbitrary")
#'
#' # Aggregate rows and compare
#' sparse_res <- utilsAggregateRows(sparse_matrix, grp)
#' dense_res <- utilsAggregateRows(dense_matrix, grp)
#' all(sparse_res == dense_res)
#'
#' # Note that storage type was converted to integers!
#' storage.mode(sparse_res)
#' storage.mode(dense_res)
#'
#' # You can also elect to keep a sparse representation
#' utilsAggregateRows(sparse_matrix, grp, sparse = TRUE)
#'
#' #### Examples with unused levels ####
#'
#' # Silly example
#' dense_mat <- replicate(5, runif(10))
#' sparse_mat <- as(dense_mat, "dgCMatrix")
#' fct_unused <- factor(c(1, 1, NA, NA, 3, 3, NA, NA, 5, 5), levels=1:5)
#'
#' # The default is to drop unused levels
#' utilsAggregateRows(dense_mat, fct_unused, drop=TRUE)
#' utilsAggregateRows(sparse_mat, fct_unused, drop=TRUE)
#'
#' # For dgCMatrix, one can elect to retain these:
#' utilsAggregateRows(sparse_mat, fct_unused, drop=FALSE)
#'
#' # For matrix, a warning is produced if either drop or sparse is requested
#' utilsAggregateRows(dense_mat, fct_unused, drop=FALSE)
#' utilsAggregateRows(dense_mat, fct_unused, sparse=TRUE)
setGeneric("utilsAggregateRows", function(x, group, drop=TRUE, sparse=FALSE) {
    standardGeneric("utilsAggregateRows")
})

#' @rdname utilsAggregateRows
#' @export
setMethod("utilsAggregateRows", signature(x = "matrix"),
          function(x, group, drop=TRUE, sparse=FALSE) {
              # Pre-checks
              assert_that(nrow(x) == length(group),
                          is.flag(drop),
                          is.flag(sparse))

              # Drop sparse doesn't work on raw matrix
              if (isFALSE(drop) | isTRUE(sparse)) {
                  warning("Aggregating base::matrix always drops unused levels",
                          " and return a base::matrix!")
              }

              # Aggregate
              o <- suppressWarnings(rowsum(x, group))

              # Discard NAs
              if (anyNA(group)) {
                  o <- o[-nrow(o), ]
              }

              # Try and coerce output
              if(storage.mode(o)  != "integer"){
                  if (all(o == floor(o))){
                      storage.mode(o) <- "integer"
                  }
              }

              # Check output
              group <- factor(group)
              stopifnot(nrow(o) == length(levels(group)),
                        rownames(o) == levels(group))

              # Return
              o
          })

#' @rdname utilsAggregateRows
#' @importClassesFrom Matrix dgCMatrix
#' @export
setMethod("utilsAggregateRows", signature(x = "dgCMatrix"),
          function(x, group, drop=TRUE, sparse=FALSE) {
              # Pre-checks
              assert_that(nrow(x) == length(group),
                          is.flag(drop),
                          is.flag(sparse))

              # Use matrix multiplixation
              o <- Matrix::fac2sparse(group, drop.unused.levels = drop)
              o <- o %*% x

              # Coerce to dense integer matrix
              if(isFALSE(sparse)){
                  o <- methods::as(o, "matrix")
                  if(storage.mode(o)  != "integer"){
                      if (all(o == floor(o))){
                          storage.mode(o) <- "integer"
                      }
                  }
              }

              # Return
              o
          })

#### Annotation utilities ####

#' Utility: Extract annotation hierachy from a TxDb.
#'
#' Used by assignTxType. This function extracts the hierachical annotations used
#' by assignTxType from a TxDb object. If you are annotating many ranges, it can
#' be time saving to built the hierachy first, to avoid processing the TxDb for
#' every assignTxDb call.
#'
#' @param object TxDb: Transcript database
#' @param tssUpstream integer: Distance to extend annotated promoter upstream.
#' @param tssDownstream integer: Distance to extend annotated promoter
#'   downstream.
#' @param proximalUpstream integer: Maximum distance upstream of promoter to be
#'   considered proximal.
#' @param detailedAntisense logical: Wether to mirror all txType categories in
#'   the antisense direction (TRUE) or lump them all together (FALSE).
#'
#' @return GRangesList of annotation hierachy
#' @family Utility functions
#' @seealso assignTxType
#' @export
#'
#' @examples
#' \dontrun{
#' data(exampleUnidirectional)
#'
#' # Obtain transcript models from a TxDb-object:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'
#' # Simplify txdb
#' hierachy <- utilsSimplifyTxDb(txdb)
#'
#' # Standard way of calling
#' x <- assignTxType(exampleUnidirectional,
#'                   txModels=txdb)
#'
#' # Calling with premade hierachy
#' y <- assignTxType(exampleUnidirectional, txModels=hierachy)
#'
#' # These are identical
#' stopifnot(all(rowRanges(x)$txType == rowRanges(y)$txType))
#' }
utilsSimplifyTxDb <- function(object,
                              tssUpstream = 100,
                              tssDownstream = 100,
                              proximalUpstream = 1000,
                              detailedAntisense = FALSE){
    # Pre-checks
    assert_that(methods::is(object, "TxDb"),
                is.count(tssUpstream),
                is.count(tssDownstream),
                is.count(proximalUpstream),
                is.flag(detailedAntisense))

    # Build sense hierachy
    hierachy <- GRangesList(
        promoter = granges(trim(promoters(object, upstream = tssUpstream,
                                          downstream = tssDownstream))),
        proximal = granges(trim(promoters(object, upstream = proximalUpstream,
                                          downstream = 0))),
        fiveUTR = granges(unlist(fiveUTRsByTranscript(object))),
        threeUTR = granges(unlist(threeUTRsByTranscript(object))),
        CDS = granges(cds(object)),
        exon = granges(exons(object)),
        intron = granges(unlist(intronsByTranscript(object))))

    # Build antisense hierachy message('Adding antisense categories...')
    if (detailedAntisense) {
        antisense <- invertStrand(hierachy)
        names(antisense) <- paste0("antisense_", names(antisense))
        hierachy <- c(hierachy, antisense)
    } else if (!detailedAntisense) {
        antisense <- invertStrand(granges(transcripts(object)))
        hierachy$antisense <- antisense
    } else {
        stop("detailedAntisense must be either TRUE/FALSE!")
    }

    #  Return
    hierachy
}
