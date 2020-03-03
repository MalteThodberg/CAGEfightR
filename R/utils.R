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

extendRanges <- function(object, upstream, downstream) {
    U <- flank(object, width = upstream, start = TRUE)
    D <- flank(object, width = downstream, start = FALSE)
    object <- punion(object, U)
    object <- punion(object, D)
    object
}

#### Counting Utilities ####

countScoredOverlaps <- function(query, subject) {
    # Solution from bioconductor by M Morgan
    hits <- methods::as(findOverlaps(query = query, subject = subject), "List")
    weightedCount <- sum(extractList(score(subject), hits))

    # Return
    weightedCount
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
