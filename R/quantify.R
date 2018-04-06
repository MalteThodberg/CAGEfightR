#### CTSSs ####

gf_mapper <- function(range, file, seqinfo, strand = "*") {
    # Load and force seqlevels
    o <- import(file, which = range)
    seqlevels(o) <- seqlevels(seqinfo)
    # isCircular(o) <- isCircular(seqinfo) genome(o) <- genome(seqinfo)
    seqinfo(o) <- seqinfo
    strand(o) <- strand
    
    # Compress and rename scores
    score(o) <- as.integer(score(o))
    
    # Return (GPos in comments) methods::as(o, 'GPos') #
    o
}

#' @importClassesFrom Matrix dgCMatrix
gf_reducer <- function(mapped) {
    # Save names
    j_names <- names(mapped)
    
    # Add j mapped <- mapply(function(gp, j){gp$j <- j; gp}, mapped,
    # seq_along(mapped))#
    mapped <- mapply(function(gp, j) {
        mcols(gp)[, "j"] <- j
        gp
    }, mapped, seq_along(mapped))
    
    # Coerce to mat mapped <- do.call(pryr::partial(c, x=GPos()), mapped)
    mapped <- do.call(pryr::partial(c, x = GRanges()), mapped)
    
    # Add i (erase mcols to save memory) dj <- unique(granges(mapped)) # Old version,
    # coerced to GRanges
    dj <- mapped
    mcols(dj) <- NULL
    dj <- unique(dj)
    dj <- sort(dj)
    
    mapped$i <- match(mapped, dj)
    
    # To matrix
    if (length(mapped) > 0) {
        mapped <- Matrix::sparseMatrix(i = mapped$i, j = mapped$j, x = score(mapped), 
            dimnames = list(NULL, j_names))
    } else {
        mapped <- methods::as(matrix(nrow = 0, ncol = length(j_names), dimnames = list(NULL, 
            j_names)), "dgCMatrix")
    }
    
    # Assemble into summarized experiment
    mapped <- SummarizedExperiment(assays = SimpleList(counts = mapped), rowRanges = dj)
    
    # Return
    mapped
}

gf_wrapper <- function(files, ranges, seqinfo, strand) {
    # Run GenomicFiles
    o <- GenomicFiles::reduceByRange(ranges = ranges, files = files, MAP = pryr::partial(gf_mapper, 
        seqinfo = seqinfo, strand = strand), REDUCE = gf_reducer, iterate = FALSE)
    
    # Merge output
    o <- do.call(rbind, o)
    
    # Return
    o
}

#' Quantify CAGE Transcriptions Start Sites (CTSSs)
#'
#' This function reads in CTSS count data from a series of BigWig-files and
#' returns a CTSS-by-library count matrix. For efficient processing, the count
#' matrix is stored as a sparse matrix (dgCMatrix).
#'
#' @param plusStrand BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param design DataFrame or data.frame: Additional information on samples.
#' @param genome Seqinfo: Genome information. If NULL the smallest common genome
#'   will be found using bwCommonGenome.
#' @param tileWidth integer: Size of tiles to parallelize over.
#'
#' @return RangedSummarizedExperiment, where assay is a sparse matrix
#'   (dgCMatrix) of CTSS counts..
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
#' CTSSs <- quantifyCTSSs(plusStrand=bw_plus,
#'                        minusStrand=bw_minus,
#'                        design=exampleDesign,
#'                        genome=si)
#'
#' # Quantification can be speed up by using multiple cores:
#' library(BiocParallel)
#' register(MulticoreParam(workers=3))
#' CTSSs <- quantifyCTSSs(plusStrand=bw_plus,
#'                        minusStrand=bw_minus,
#'                        design=exampleDesign,
#'                        genome=si)
#' }
quantifyCTSSs <- function(plusStrand, minusStrand, design = NULL, genome = NULL, 
    tileWidth = 100000000L) {
    # Pre-checks
    assert_that(class(plusStrand) == "BigWigFileList", class(minusStrand) == "BigWigFileList", 
        length(plusStrand) == length(minusStrand), bwValid(plusStrand), bwValid(minusStrand), 
        all(names(plusStrand) == names(minusStrand)), is.count(tileWidth))
    
    # Set design
    if (is.null(design)) {
        design <- DataFrame(row.names = names(plusStrand))
    } else if (class(design) == "DataFrame") {
        assert_that(all(rownames(design) == names(plusStrand)))
    } else if (is.data.frame(design)) {
        assert_that(all(rownames(design) == names(plusStrand)))
        design <- methods::as(design, "DataFrame")
    } else {
        stop("design must either NULL or a DataFrame-/data.frame-object!")
    }
    
    # Aquire seqinfo if missing.
    if (is.null(genome)) {
        message("Finding common genome...")
        genome <- bwCommonGenome(plusStrand, minusStrand, method = "intersect")
    } else if (class(genome) == "Seqinfo") {
        message("Checking supplied genome compatibility...")
        bwGenomeCompatibility(plusStrand, minusStrand, genome)
    } else {
        stop("genome must either NULL or a Seqinfo-object!")
    }
    
    # Setup tiles
    grl <- GenomicRanges::tileGenome(genome, tilewidth = tileWidth)
    message("Iterating over ", length(grl), " genomic tiles in ", length(plusStrand), 
        " samples using ", BiocParallel::bpworkers(), " worker(s)...")
    
    # Load data
    message("Importing CTSSs from plus strand...")
    plus_strand <- gf_wrapper(files = plusStrand, ranges = grl, seqinfo = genome, 
        strand = "+")
    
    message("Importing CTSSs from minus strand...")
    minus_strand <- gf_wrapper(files = minusStrand, ranges = grl, seqinfo = genome, 
        strand = "-")
    
    # Merge
    message("Merging strands...")
    o <- rbind(plus_strand, minus_strand)
    rm(plus_strand, minus_strand)
    
    # Attach design
    colData(o) <- design
    
    # Post-checks (GPos commented out)
    stopifnot(class(rowRanges(o)) == "GRanges", length(plusStrand) == ncol(o), identical(seqinfo(o), 
        genome))
    
    message("### CTSS summary ###")
    message("Number of samples: ", ncol(o))
    message("Number of CTSSs: ", format(nrow(o)/1000000L, digits = 4), " millions")
    message("Sparsity: ", format((1 - (Matrix::nnzero(assay(o))/length(assay(o)))) * 
        100, digits = 4), " %")
    message("Final object size: ", utils::capture.output(pryr::object_size(o)))
    
    # Return
    o
}

#### TCs and Genes ####

#' Quantify expression of clusters (TSSs or enhancers) by summing CTSSs within
#' clusters.
#'
#' @param object RangedSummarizedExperiment: CTSSs.
#' @param clusters GRanges: Clusters to be quantified.
#' @param inputAssay character: Name of assay holding expression values to be
#'   quantified (usually counts).
#' @param sparse logical: If the input is a sparse matrix, TRUE will keep the
#'   output matrix sparse while FALSE will coerce it into a normal matrix.
#'
#' @return RangedSummarizedExperiment with row corresponding to clusters.
#'   seqinfo and colData is copied over from object.
#'
#' @family Quantification functions
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' # CTSSs stored in a RangedSummarizedExperiment:
#' data(exampleCTSS)
#'
#' # Clusters to be quantified as a GRanges:
#' data(exampleUnidirectional)
#' clusters <- rowRanges(exampleUnidirectional)
#'
#' # Quantify clusters:
#' quantifyClusters(exampleCTSSs, clusters)
#'
#' # For exceptionally large datasets,
#' # the resulting count matrix can be left sparse:
#' quantifyClusters(exampleCTSSs, rowRanges(exampleUnidirectional), sparse=TRUE)
quantifyClusters <- function(object, clusters, inputAssay = "counts", sparse = FALSE) {
    # Pre-checks
    assert_that(class(object) == "RangedSummarizedExperiment", not_empty(object), 
        isDisjoint(object), methods::is(clusters, "GRanges"), not_empty(clusters), 
        isDisjoint(clusters), identical(seqinfo(object), seqinfo(clusters)), is.string(inputAssay), 
        inputAssay %in% assayNames(object), is.flag(sparse))
    
    # Find overlaps
    message("Finding overlaps...")
    hits <- findOverlaps(query = object, subject = clusters, select = "arbitrary")
    hits <- factor(hits, levels = seq_along(clusters))
    missing_hits <- is.na(hits)
    
    # Warn if there is no overlap
    if (all(missing_hits)) {
        warning("The supplied clusters had no overlapping CTSSs!")
    }
    
    # Summarize
    message("Aggregating within clusters...")
    mat <- rowsum2(x = assay(object, inputAssay), group = hits, drop = FALSE, sparse = sparse)
    
    # Check output is the right format and assign names.
    stopifnot(nrow(mat) == length(clusters))
    rownames(mat) <- names(clusters)
    
    # Coerce to RSE
    o <- SummarizedExperiment(assays = SimpleList(mat), rowRanges = clusters, colData = colData(object))
    assayNames(o) <- inputAssay
    
    # Return
    o
}

#' Quantify expression of genes
#'
#' Obtain gene-level expression estimates by summing clusters annotated to the
#' same gene. Unannotated transcripts (NAs) are discarded.
#'
#' @param object RangedSummarizedExperiment: Cluster-level expression values.
#' @param genes character: Name of column in rowData holding gene IDs (NAs will
#'   be discarded).
#' @param inputAssay character: Name of assay holding values to be quantified,
#'   (usually counts).
#' @param sparse logical: If the input is a sparse matrix, TRUE will keep the
#'   output matrix sparse while FALSE will coerce it into a normal matrix.
#'
#' @return RangedSummarizedExperiment with rows corresponding to genes. Location
#'   of clusters within genes is stored as a GRangesList in rowRanges. seqinfo
#'   and colData is copied over from object.
#' @family Quantification functions
#' @export
#' @examples
#' data(exampleUnidirectional)
#'
#' # Annotate clusters with geneIDs:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#' exampleUnidirectional <- assignGeneID(exampleUnidirectional,
#'                                       geneModels=txdb,
#'                                       outputColumn='geneID')
#'
#' # Quantify counts within genes:
#' quantifyGenes(exampleUnidirectional, genes='geneID', inputAssay='counts')
#'
#' # For exceptionally large datasets,
#' # the resulting count matrix can be left sparse:
#' quantifyGenes(exampleUnidirectional,
#'               genes='geneID',
#'               inputAssay='counts',
#'               sparse=TRUE)
quantifyGenes <- function(object, genes, inputAssay = "counts", sparse = FALSE) {
    # Pre-checks
    assert_that(methods::is(object, "RangedSummarizedExperiment"), isDisjoint(object), 
        is.string(genes), is.element(genes, colnames(rowData(object))), is.character(rowData(object)[, 
            genes]), is.string(inputAssay), inputAssay %in% assayNames(object), is.flag(sparse))
    
    # Factor genes
    genes <- factor(rowData(object)[, genes])
    
    # Split ranges
    new_gr <- splitAsList(rowRanges(object), f = genes, drop = TRUE)
    
    # Sum matrix
    new_m <- rowsum2(assay(object, inputAssay), group = genes, drop = TRUE, sparse = sparse)
    
    # Check that names match
    stopifnot(setequal(rownames(new_m), names(new_gr)))
    if (!all(rownames(new_m) == names(new_gr))) {
        new_m <- new_m[names(new_gr), ]
    }
    
    # Reassemble and copy over
    o <- SummarizedExperiment(assays = list(new_m), rowRanges = new_gr, colData = colData(object), 
        metadata = metadata(object))
    assayNames(o) <- inputAssay
    
    # Calculate some stats
    rowData(o)[, "nClusters"] <- elementNROWS(new_gr)
    
    # Return
    o
}
