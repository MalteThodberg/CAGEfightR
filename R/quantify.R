### CTSS

#' @import rtracklayer GenomeInfoDb GenomicRanges
gf_mapper <- function(range, file, seqinfo, strand="*"){
	# Load and force seqlevels
	o <- import(file, which=range)
	seqlevels(o) <- seqlevels(seqinfo)
	# isCircular(o) <- isCircular(seqinfo)
	# genome(o) <- genome(seqinfo)
	seqinfo(o) <- seqinfo
	strand(o) <- strand

	# Compress and rename scores
	score(o) <- as.integer(score(o))

	# Return as GenomicPositions
	#methods::as(o, "GPos")
	o
}

#' @importClassesFrom Matrix dgCMatrix
#' @import GenomicRanges SummarizedExperiment
gf_reducer <- function(mapped){
	# Save names
	j_names <- names(mapped)

	# Add j
	#mapped <- mapply(function(gp, j){gp$j <- j; gp}, mapped, seq_along(mapped))#
	mapped <- mapply(function(gp, j){mcols(gp)[,"j"] <- j; gp}, mapped, seq_along(mapped))

	# Coerce to mat
	#mapped <- do.call(pryr::partial(c, x=GPos()), mapped)
	mapped <- do.call(pryr::partial(c, x=GRanges()), mapped)

	# Add i (erase mcols to save memory)
	#dj <- unique(granges(mapped))  # Old  version,  coerced to  GRanges
	dj <- mapped
	mcols(dj) <- NULL
	dj <- unique(dj)
	dj <- sort(dj)

	mapped$i <- match(mapped, dj)

	# To matrix
	if(length(mapped) > 0){
		mapped <- Matrix::sparseMatrix(i=mapped$i, j=mapped$j, x=score(mapped), dimnames = list(NULL, j_names))
	}else{
		mapped <- methods::as(matrix(nrow=0, ncol=length(j_names), dimnames=list(NULL, j_names)), "dgCMatrix")
	}

	# Assemble into summarized experiment
	mapped <- SummarizedExperiment(assays=SimpleList(counts=mapped),
														rowRanges=dj)

	# Return
	mapped
}

#' @import GenomeInfoDb GenomicRanges SummarizedExperiment
gf_wrapper <- function(files, ranges, seqinfo, strand){
	# Run GenomicFiles
	o <- GenomicFiles::reduceByRange(ranges=ranges, files=files,
										 MAP=pryr::partial(gf_mapper, seqinfo=seqinfo, strand=strand),
										 REDUCE=gf_reducer,
										 iterate=FALSE)

	# Merge output
	o <- do.call(rbind, o)

	# Return
	o
}

#' Quantify CAGE Transcriptions Start Sites (CTSSs)
#'
#' This function reads in CTSS data from a series of BigWig-files and returns a basepair-by-sample count matrix. To save memory, the count matrix is stored as a sparse matrix (dgCMatrix), and basepair positions as GPos-object.
#'
#' @param plusStrand BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param design DataFrame or data.frame: Additional information on samples.
#' @param genome Seqinfo: Genome information. If NULL the smallest common genome will be found using bwCommonGenome.
#' @param tileWidth integer: Size of tiles to parallelize over.
#'
#' @return RangedSummarizedExperiment, where assay is a sparse matrix (dgCMatrix) of CTSS counts and rowRanges is GPos-object.
#'
#' @family Quantification functions
#' @importClassesFrom Matrix dgCMatrix
#' @import S4Vectors rtracklayer SummarizedExperiment assertthat
#' @export
quantifyCTSSs <- function(plusStrand, minusStrand, design=NULL, genome=NULL, tileWidth=1e8L){
	# Pre-checks
	assert_that(class(plusStrand) == "BigWigFileList",
							class(minusStrand) == "BigWigFileList",
							length(plusStrand) == length(minusStrand),
							bwValid(plusStrand),
							bwValid(minusStrand),
							all(names(plusStrand) == names(minusStrand)),
							is.count(tileWidth))

	# Set design
	if(is.null(design)){
		design <- DataFrame(row.names=names(plusStrand))
	}else if(class(design) == "DataFrame"){
		assert_that(all(rownames(design) == names(plusStrand)))
	}else if(is.data.frame(design)){
		assert_that(all(rownames(design) == names(plusStrand)))
		design <- methods::as(design, "DataFrame")
	}else{
		stop("design must either NULL or a DataFrame-/data.frame-object!")
	}

	# Aquire seqinfo if missing.
	if(is.null(genome)){
		message("Finding common genome...")
		genome <- bwCommonGenome(plusStrand, minusStrand, method="intersect")
	}else if(class(genome) == "Seqinfo"){
		message("Checking supplied genome compatibility...")
		bwGenomeCompatibility(plusStrand, minusStrand, genome)
	}else{
		stop("genome must either NULL or a Seqinfo-object!")
	}

	# Setup tiles
	grl <- GenomicRanges::tileGenome(genome, tilewidth=tileWidth)
	message("Iterating over ", length(grl), " genomic tiles in ", length(plusStrand), " samples using ", BiocParallel::bpworkers(), " worker(s)...")

	# Load data
	message("Importing CTSSs from plus strand...")
	plus_strand <- gf_wrapper(files=plusStrand, ranges=grl, seqinfo=genome, strand="+")

	message("Importing CTSSs from minus strand...")
	minus_strand <- gf_wrapper(files=minusStrand, ranges=grl, seqinfo=genome, strand="-")

	# Merge and realize
	message("Merging strands...")
	o <- rbind(plus_strand, minus_strand)
	rm(plus_strand, minus_strand)

	# Post-checks
	stopifnot(class(rowRanges(o)) == "GRanges",
						#class(rowRanges(o)) == "GPos", Comment  in for GPos
						length(plusStrand) == ncol(o),
						identical(seqinfo(o), genome))

	message("### CTSS summary ###")
	message("Number of samples: ", ncol(o))
	message("Number of CTSSs: ", format(nrow(o) / 1e6L, digits=4), " millions")
	message("Sparsity: ", format((1 - (Matrix::nnzero(assay(o)) / length(assay(o)))) * 100, digits=4), " %")
	message("Final object size: ", utils::capture.output(pryr::object_size(o)))

	# Return
	o
}

### TCs and Genes

#' Quantify Tag Clusters (TCs)
#'
#' @param object RangedSummarizedExperiment: CTSSs.
#' @param clusters GRanges: Clusters ro be quantified.
#' @param inputAssay character: Name of holding values to be quantified
#'
#' @return RangedSummarizedExperiment with quantified CTSSs for each region in cluster, with seqinfo and colData is copied over.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @import S4Vectors GenomicRanges rtracklayer SummarizedExperiment assertthat
#' @export
quantifyClusters <- function(object, clusters, inputAssay="counts"){
	# Pre-checks
	assert_that(class(object) == "RangedSummarizedExperiment",
							not_empty(object),
							isDisjoint(object),
							class(clusters) == "GRanges",
							not_empty(clusters),
							isDisjoint(clusters),
							identical(seqinfo(object), seqinfo(clusters)),
							is.string(inputAssay),
							inputAssay %in% assayNames(object))

	# Find overlaps
	message("Finding overlaps...")
	hits <- findOverlaps(query=object, subject=clusters, select="arbitrary")
	missing_hits <- is.na(hits)

	# Warn if there is no overlap
	if(all(missing_hits)){
		warning("The supplied clusters had no overlapping CTSSs!")
	}

	# Summarize
	message("Aggregating within clusters...")
	mat <- Matrix.utils::aggregate.Matrix(x=assay(object, inputAssay), groupings=hits, fun="sum")

	# Coerce to integer matrix
	message("Preparing output...")
	# Coerce to basic matrix
	mat <- Matrix::as.matrix(mat)

	# Remove NA row if present
	if(any(missing_hits)){
		mat <- mat[-nrow(mat),]
	}

	# Simplify if possible
	if(all(mat == floor(mat))){
		storage.mode(mat) <- "integer"
	}

	# Check output is the right format and assign names.
	stopifnot(nrow(mat) == length(clusters))
	rownames(mat) <- names(clusters)

	# Coerce to RSE
	o <- SummarizedExperiment(assays=SimpleList(mat),
														rowRanges=clusters,
														colData=colData(object))
	assayNames(o) <- inputAssay

	# Return
	o
}

#' Quantify genes
#'
#' Obtain a gene-level EM by summing transcripts within genes. Unannotated transcripts (NAs) are discarded.
#'
#' @param object RangedSummarizedExperiment: Transcript-level counts.
#' @param genes character: Name of column in rowData holding genes (NAs will be discarded).
#' @param inputAssay character: Name of assay holding values to be quantified, usually raw counts.
#'
#' @return RangedSummarizedExperiment with rows corresponding to genes. All other information is copied over from object.
#' @import assertthat S4Vectors GenomicRanges SummarizedExperiment
#' @export
quantifyGenes <- function(object, genes, inputAssay="counts"){
	# Pre-checks
	assert_that(methods::is(object, "RangedSummarizedExperiment"),
							isDisjoint(object),
							is.string(genes),
							is.element(genes, colnames(rowData(object))),
							is.character(rowData(object)[,genes]),
							is.string(inputAssay),
							inputAssay %in% assayNames(object))

	# Split ranges
	new_gr <- splitAsList(rowRanges(object), rowData(object)[,genes], drop=TRUE)

	# Sum matrix
	new_m <- suppressWarnings(rowsum(assay(object, inputAssay), group=rowData(object)[,genes]))

	# Discard NA
	new_m <- new_m[!is.na(rownames(new_m)),]

	# Check that names match
	stopifnot(setequal(rownames(new_m), names(new_gr)))
	if(!all(rownames(new_m) == names(new_gr))){
		new_m <- new_m[names(new_gr),]
	}

	# Reassemble and copy over
	o <- SummarizedExperiment(assays=list(new_m),
														rowRanges=new_gr,
														colData=colData(object),
														metadata=metadata(object))
	assayNames(o) <- inputAssay

	# Calculate some stats
	rowData(o)[,"nFeatures"] <- elementNROWS(new_gr)

	# Return
	o
}
