#### Mappers ####

bg_mapper <- function(file, seqinfo, strand = "*"){
    # Import
    o <- import.bedGraph(file, genome=seqinfo)

    # Instead of conversion, stop if not proper CTSSs-like
    stopifnot(all(width(o) == 1))
    # o <- as(o, "StitchedGPos")
    # score(o) <- as.integer(score(o))

    # Set strand
    strand(o) <- strand

    # Return
    o
}

bw_mapper <- function(range, file, seqinfo, strand = "*"){
    # Import
    o <- suppressWarnings(import.bw(con=file,
                                    format="bigwig",
                                    which=range,
                                    as="GRanges"))

    # Instead of conversion, stop if not proper CTSSs-like
    stopifnot(all(width(o) == 1))
    # o <- as(o, "StitchedGPos")
    # score(o) <- as.integer(score(o)) # Not used by dgC anyway

    # Force seqlevels and strand
    seqlevels(o) <- seqlevels(seqinfo)
    seqinfo(o) <- seqinfo
    strand(o) <- strand

    # Return
    o
}

#### Reducers ####

gf_reducer2 <- function(mapped){
    # Make table of j = sample indicies
    o <- unlist(GRangesList(mapped), use.names=FALSE)
    o$j <- rep(seq_along(mapped), elementNROWS(mapped))

    # Find unique sites
    unique_pos <- granges(o)
    unique_pos <- unique(unique_pos)
    unique_pos <- sort(unique_pos)

    # Overlap to get i = position indices
    o$i <- match(o, unique_pos)

    # Determine the shape of the sparse matrix
    mat_dim <- c(length(unique_pos), length(mapped))
    mat_names <- list(NULL, names(mapped))

    # Build sparse matrix
    if (length(o) > 0) {
      o <- Matrix::sparseMatrix(i = o$i,
                                     j = o$j,
                                     x = score(o),
                                     dims = mat_dim,
                                     dimnames = mat_names)
    } else {
      # Simply create an empty matrix
      o <- Matrix::sparseMatrix(i = integer(0),
                                j = integer(0),
                                x = integer(0),
                                dims = c(0, mat_dim[2]),
                                dimnames = mat_names)
    }

    # # Simplify so sparse matrix encoding
    # indices <- as.matrix(mcols(indices)[,c("i", "j", "score")])
    # rle_dim <- c(length(unique_pos), length(mapped))
    # o <- indices[,1:2]
    # o <- SparseArraySeed(dim=rle_dim,
    #                      aind=o,
    #                      nzdata=indices[,3])
    #
    # # To DataFrame of Rles via RleMatrix
    # # o <- as(o, "RleMatrix")
    # o <- DelayedArray(o)
    # o <- as(o, "DataFrame")

    # Build RSE
    o <- SummarizedExperiment(assays=list(counts=o),
                              rowRanges=unique_pos)
    colnames(o) <- names(mapped)

    # Return
    o
}

gf_wrapper2 <- function(files, ranges, seqinfo, strand) {
    # Run GenomicFiles
    o <- GenomicFiles::reduceByRange(ranges = ranges,
                                     files = files,
                                     MAP = pryr::partial(bw_mapper,
                                                         seqinfo = seqinfo,
                                                         strand = strand),
                                     REDUCE = gf_reducer2,
                                     iterate = FALSE)

    # Merge output
    #message("Attempting to bind chunks...")
    o <- do.call(rbind, o)

    # Return
    o
}

#### Other helpers ####

setGeneric("cleanDesign", function(design, SampleIDs) {
    standardGeneric("cleanDesign")
})

setMethod("cleanDesign", signature(design = "NULL"), function(design, SampleIDs) {
    DataFrame(row.names = SampleIDs)
})

setMethod("cleanDesign", signature(design = "DataFrame"), function(design, SampleIDs) {
    assert_that(nrow(design) == length(SampleIDs),
                all(rownames(design) == SampleIDs))
    design
})

setMethod("cleanDesign", signature(design = "data.frame"), function(design, SampleIDs) {
    assert_that(nrow(design) == length(SampleIDs),
                all(rownames(design) == SampleIDs))
    methods::as(design, "DataFrame")
})

setMethod("cleanDesign", signature(design = "ANY"), function(design, SampleIDs) {
    stop("design must either NULL or a DataFrame-/data.frame-object!")
})

# pseudoSparsity <- function(x){
#     stopifnot(methods::is(x, "DataFrame"))
#
#     # Find non-zero
#     o <- vapply(X=x,
#                 FUN=function(x) sum(x == 0),
#                 FUN.VALUE=integer(1))
#
#     # Calc sparsity
#     o <- sum(o) / (ncol(x) * nrow(x))
#
#     # Return
#     o
# }

# clean_design <- function(ids, design=NULL){
#     # Decide what to
#     if (is.null(design)) {
#         design <- DataFrame(row.names = ids)
#     } else if (methods::is(design, "DataFrame")) {
#         assert_that(all(rownames(design) == ids))
#     } else if (is.data.frame(design)) {
#         assert_that(all(rownames(design) == ids))
#         design <- methods::as(design, "DataFrame")
#     } else {
#         stop("design must either NULL or a DataFrame-/data.frame-object!")
#     }
#
#     # Return
# }

#### Main function ####

#' Quantify CAGE Transcriptions Start Sites (CTSSs)
#'
#' This function reads in CTSS count data from a series of BigWig-files (or bedGraph-files) and
#' returns a CTSS-by-library count matrix. For efficient processing, the count
#' matrix is stored as a sparse matrix (dgCMatrix from the Matrix package), and CTSSs are compressed to a GPos object if possible.
#'
#' @param plusStrand BigWigFileList or character: BigWig/bedGraph files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList or character: BigWig/bedGraph files with minus-strand CTSS data.
#' @param design DataFrame or data.frame: Additional information on samples which will be added to the ouput
#' @param genome Seqinfo: Genome information. If NULL the smallest common genome
#'   will be found using bwCommonGenome when BigWig-files are analyzed.
#' @param nTiles integer: Number of genomic tiles to parallelize over.
#' @param ... additional arguments passed to methods.
#'
#' @return RangedSummarizedExperiment, where assay is a sparse matrix
#'   (dgCMatrix) of CTSS counts and design stored in colData.
#' @family Quantification functions
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' \dontrun{
#' # Load the example data
#' data('exampleDesign')
#' # Use the BigWig-files included with the package:
#' bw_plus <- system.file('extdata', exampleDesign$BigWigPlus,
#'                        package = 'CAGEfightR')
#' bw_minus <- system.file('extdata', exampleDesign$BigWigMinus,
#'                         package = 'CAGEfightR')
#'
#' # Create two named BigWigFileList-objects:
#' bw_plus <- BigWigFileList(bw_plus)
#' bw_minus <- BigWigFileList(bw_minus)
#' names(bw_plus) <- exampleDesign$Name
#' names(bw_minus) <- exampleDesign$Name
#'
#' # Quantify CTSSs, by default this will use the smallest common genome:
#' CTSSs <- quantifyCTSSs(plusStrand=bw_plus,
#'                        minusStrand=bw_minus,
#'                        design=exampleDesign)
#'
#' # Alternatively, a genome can be specified:
#' si <- seqinfo(bw_plus[[1]])
#' si <- si['chr18']
#' CTSSs_subset <- quantifyCTSSs(plusStrand=bw_plus,
#'                        minusStrand=bw_minus,
#'                        design=exampleDesign,
#'                        genome=si)
#'
#' # Quantification can be speed up by using multiple cores:
#' library(BiocParallel)
#' register(MulticoreParam(workers=3))
#' CTSSs_subset <- quantifyCTSSs(plusStrand=bw_plus,
#'                        minusStrand=bw_minus,
#'                        design=exampleDesign,
#'                        genome=si)
#'
#' # CAGEfightR also support bedGraph files, first BigWig is converted
#' bg_plus <- replicate(n=length(bw_plus), tempfile(fileext="_plus.bedGraph"))
#' bg_minus <- replicate(n=length(bw_minus), tempfile(fileext="_minus.bedGraph"))
#' names(bg_plus) <- names(bw_plus)
#' names(bg_minus) <- names(bw_minus)
#'
#' convertBigWig2BedGraph(input=sapply(bw_plus, resource), output=bg_plus)
#' convertBigWig2BedGraph(input=sapply(bw_minus, resource), output=bg_minus)
#'
#' # Then analyze: Note a genome MUST be supplied here!
#' si <- bwCommonGenome(bw_plus, bw_minus)
#' CTSSs_via_bg <- quantifyCTSSs(plusStrand=bg_plus,
#'                         minusStrand=bg_minus,
#'                         design=exampleDesign,
#'                         genome=si)
#'
#' # Confirm that the two approaches yield the same results
#' all(assay(CTSSs_via_bg) == assay(CTSSs))
#' }
setGeneric("quantifyCTSSs", function(plusStrand, minusStrand, design=NULL, genome=NULL, ...) {
    standardGeneric("quantifyCTSSs")
})

#' @rdname quantifyCTSSs
#' @export
setMethod("quantifyCTSSs",
          signature(plusStrand = "BigWigFileList",
                    minusStrand = "BigWigFileList"),
  function(plusStrand, minusStrand, design=NULL, genome=NULL, nTiles=1L) {
      # Pre-checks
      assert_that(length(plusStrand) == length(minusStrand),
                  bwValid(plusStrand),
                  bwValid(minusStrand),
                  !is.null(names(plusStrand)),
                  !is.null(names(minusStrand)),
                  all(names(plusStrand) == names(minusStrand)),
                  is.count(nTiles))

      # Set design
      message("Checking design...")
      design <- cleanDesign(design, SampleIDs=names(plusStrand))

      # Aquire seqinfo if missing.
      if (is.null(genome)) {
          message("Finding common genome...")
          genome <- bwCommonGenome(plusStrand, minusStrand, method = "intersect")
      } else if (methods::is(genome, "Seqinfo")) {
          message("Checking supplied genome compatibility...")
          bwGenomeCompatibility(plusStrand, minusStrand, genome)
      } else {
          stop("genome must either NULL or a Seqinfo-object!")
      }

      # Setup tiles and workers
      grl <- tileGenome(genome, ntile=nTiles)
      n_workers <- BiocParallel::bpworkers()
      n_workers <- ifelse(is.numeric(n_workers), n_workers, 1)

      message("Iterating over ", length(grl),
              " genomic tiles in ", length(plusStrand),
              " samples using ", n_workers, " worker(s)...")

      # Load data
      message("Importing CTSSs from plus strand...")
      plus_strand <- gf_wrapper2(files = plusStrand,
                                 ranges = grl,
                                 seqinfo = genome,
                                 strand = "+")

      message("Importing CTSSs from minus strand...")
      minus_strand <- gf_wrapper2(files = minusStrand,
                                  ranges = grl,
                                  seqinfo = genome,
                                  strand = "-")

      # Merge and format
      message("Merging strands...")
      o <- rbind(plus_strand, minus_strand)
      rm(plus_strand, minus_strand)

      # Attach design
      colData(o) <- design

      message("Formatting output...")
      # Calculate sparsity before seed
      zero_frac <- 1 - (Matrix::nnzero(assay(o)) / length(assay(o)))

      # Wrap in DelayedArray
      #assay(o) <- DelayedArray(assay(o))

      # Convert to GPos if possible:
      if(length(o) < .Machine$integer.max){
        rowRanges(o) <- methods::as(rowRanges(o), "StitchedGPos")
      }

      # Attach design

      # Sparsity
      #  Convert to appropriate seed
      # message("Formatting to specified seed...")
      # if(seed == "DataFrame"){
      #   assays(o) <- List(counts=DelayedArray(assay(o)))
      # }else if(seed == "dgCMatrix"){
      #   assays(o) <- List(counts=DelayedArray(assay(o)))
      #   assay(o) <- DelayedArray(as(assay(o), "dgCMatrix"))
      # }else if(seed == "RleMatrix"){
      #   assays(o) <- List(counts=RleArray(assay(o))) # Chunk by column
      # }else if(seed == "RleMatrix2"){
      #   assays(o) <- List(counts=as(assay(o), "RleMatrix")) # High compression
      # }else if(seed == "SparseArraySeed"){
      #   assays(o) <- List(counts=DelayedArray(assay(o)))
      #   assay(o) <- DelayedArray(as(assay(o), "SparseArraySeed"))
      # }else{
      #   warning("Did not wrap seed in DelayedArray!")
      # }

      # Calculate sparsity
      #pseudo_sparsity <- pseudoSparsity(assay(o))

      # Final conversion
      #assays(o) <- List(counts=assayCompression(assay(o, "counts")))

      # # Calculate sparsity
      # sparse_percent <- vapply(X=assay(o),
      #                          FUN=function(x) sum(x == 0),
      #                          FUN.VALUE=integer(1))
      # sparse_percent <- sum(sparse_percent) / (ncol(o) * nrow(o))

      # Final conversion
      # if(assayStorage == "RleMatrix"){
      #     assays(o) <- List(counts=RleArray(assay(o, "counts")))
      # }else if(assayStorage == "DelayedArray") {
      #     assays(o) <- List(counts=DelayedArray(assay(o, "counts")))
      # }else if(assayStorage == "DataFrame") {
      #     warning("Output is a simple DataFrame of Rle's: ",
      #     "This is mainly intended for debugging and will not work with remaining CAGEfightR functions!")
      # }else{
      #     stop("assayStorage must be in: RleMatrix, DelayedArray, DataFrame")
      # }

      # Post warnings
      if(!length(plusStrand) == ncol(o)){
        warning("Output has the wrong dimensions! ",
                "Check input files with checkCTSSs.")
      }

      if(!identical(seqinfo(o), genome)){
        warning("Seqinfo was not properly set on output! ",
                "Check input files with checkCTSSs")
      }

      # Basic Summary
      message("### CTSS summary ###")
      message("Number of samples: ", ncol(o))
      message("Number of CTSSs: ", format(nrow(o) / 1e6L, digits = 4), " millions")
      message("Sparsity: ", format(zero_frac * 100, digits=4), " %")
      #message("Type of DelayedArray seed: dgCMatrix")
      message("Type of rowRanges: ", class(rowRanges(o)))
      message("Final object size: ", utils::capture.output(pryr::object_size(o)))

      # Return
      o
})

#' @rdname quantifyCTSSs
#' @export
setMethod("quantifyCTSSs",
          signature(plusStrand = "character",
                    minusStrand = "character"),
          function(plusStrand, minusStrand, design=NULL, genome=NULL) {
              # Pre-checks
              assert_that(all(checkFileSeries(plusStrand)),
                          all(checkFileSeries(minusStrand)),
                          all(checkExtensions(plusStrand, "bedGraph")),
                          all(checkExtensions(minusStrand, "bedGraph")),
                          length(plusStrand) == length(minusStrand),
                          !is.null(names(plusStrand)),
                          !is.null(names(minusStrand)),
                          all(names(plusStrand) == names(minusStrand)),
                          methods::is(genome, "Seqinfo"))

              # Set design
              message("Checking design...")
              design <- cleanDesign(design, SampleIDs=names(plusStrand))

              # Setup tiles and workers
              n_workers <- BiocParallel::bpworkers()
              n_workers <- ifelse(is.numeric(n_workers), n_workers, 1)

              message("Iterating over ", length(plusStrand),
                      " samples using ", n_workers, " worker(s)...")

              # Load data
              message("Importing CTSSs from plus strand...")
              plus_strand <- BiocParallel::bplapply(plusStrand,
                                                    bg_mapper,
                                                    seqinfo=genome,
                                                    strand="+")
              plus_strand <- gf_reducer2(plus_strand)

              message("Importing CTSSs from minus strand...")
              minus_strand <- BiocParallel::bplapply(minusStrand,
                                                    bg_mapper,
                                                    seqinfo=genome,
                                                    strand="-")
              minus_strand <- gf_reducer2(minus_strand)

              # Merge and format
              message("Merging strands...")
              o <- rbind(plus_strand, minus_strand)
              rm(plus_strand, minus_strand)

              # Attach design
              colData(o) <- design

              message("Formatting output...")
              # Calculate sparsity before seed
              zero_frac <- 1 - (Matrix::nnzero(assay(o)) / length(assay(o)))

              # Wrap in DelayedArray
              #assay(o) <- DelayedArray(assay(o))

              # Convert to GPos if possible:
              if(length(o) < .Machine$integer.max){
                rowRanges(o) <- methods::as(rowRanges(o), "StitchedGPos")
              }

              # Basic Summary
              message("### CTSS summary ###")
              message("Number of samples: ", ncol(o))
              message("Number of CTSSs: ", format(nrow(o) / 1e6L, digits = 4), " millions")
              message("Sparsity: ", format(zero_frac * 100, digits=4), " %")
              #message("Type of DelayedArray seed: dgCMatrix")
              message("Type of rowRanges: ", class(rowRanges(o)))
              message("Final object size: ", utils::capture.output(pryr::object_size(o)))

              # Return
              o
})
