#' Swap ranges
#'
#' Swap out the range of a GRanges-object with another IRanges-object stored inside the same object. I.e., swapping TC widths with TC peaks.
#'
#' @param gr GRanges: GRanges
#' @param column character: Name of column in GRanges containing an IRanges.
#'
#' @return GRanges with swapped ranges. Old ranges are in the column "swap"
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
swapRanges <- function(gr, column="thick"){
	# TO DO
	# Check all columns are actually ranges
	# Allow other name than swap in function

	# Copy ranges to a separate column
	gr$swap <- ranges(gr)

	# Switch in the new column
	ranges(gr) <- mcols(gr)[,column]

	gr
}

#' Annotate GRanges with transcript type
#'
#' Annotate a set of ranges in a GRanges object with transcript type (TxType) based on their genic context. Transcripts are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation.
#' @param tssUpstream integer: Distance to extend annotated promoter upstream.
#' @param tssDownstream integer: Distance to extend annotated promoter downstream.
#' @param proximalUpstream integer: Maximum distance upstream of promoter to be considered proximal.
#' @param rangeColumn character: Name of column in gr with an alternative set of ranges (For example Tag Cluster (TC) peaks)
#'
#' @return character vector of same length as gr containing TxType for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignTxType <- function(gr, txdb, tssUpstream=100, tssDownstream=100, proximalUpstream=1000, rangeColumn=NULL){
	# TO DO
	# Add arguments for promoter and upstream coordinates.
	# Option to not swap ranges
	# Pre/post checks

	# Reference categories
	message("Extracting reference categories...")
	Promoters <- promoters(txdb, upstream=tssUpstream, downstream=tssDownstream)
	Proximal <- promoters(txdb, upstream=proximalUpstream, downstream=0)
	FivePrimeUTRs <- fiveUTRsByTranscript(txdb)
	ThreePrimeUTRs <- threeUTRsByTranscript(txdb)
	Introns <- intronsByTranscript(txdb)
	Exons <- exons(txdb)
	CDSs <- cds(txdb)
	Antis <- transcripts(txdb)
	strand(Antis) <- ifelse(strand(Antis) == "+", "-", "+")

	# Swap in different range for overlaps
	if(!is.null(rangeColumn)){
		message("Swapping peaks...")
		gr <- swapRanges(gr, column="thick")
	}

	# Overlap sequentially
	message("Overlapping sequentially...")
	featureType <- rep("intergenic", length(gr))
	featureType <- ifelse(overlapsAny(query=gr, subject=Antis), "antisense", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=Introns), "intron", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=Exons), "exon", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=CDSs), "CDS", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=ThreePrimeUTRs), "threeUTR", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=FivePrimeUTRs), "fiveUTR", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=Proximal), "proximal", featureType)
	featureType <- ifelse(overlapsAny(query=gr, subject=Promoters), "promoter", featureType)

	# Return
	featureType
}

#' Annotate GRanges with gene ID
#'
#' Annotate a set of ranges in a GRanges object with Entrez Gene Identifier (ID) based on their genic context. Genes are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation
#' @param proximalUpstream integer: Maximum distance upstream of annotated promoter to be considered part of genes.
#' @param rangeColumn character: Name of column in gr with an alternative set of ranges (For example Tag Cluster (TC) peaks)
#' @param multiHits character: How to choose when multiple genes overlap.
#' @return character vector of same length as gr containing Entrez gene IDs for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignGene <- function(gr, txdb, proximalUpstream=1000, rangeColumn=NULL, multiHits="arbitrary"){
	# TO DO
	# Option for reporting more than one gene.
	# Option to discard multi strand genes.
	# Option to not swap ranges
	# Pre/post checks

	# Reference categories
	message("Extracting reference categories...")
	#Genes <- genes(txdb, single.strand.genes.only=TRUE)
	Genes <- transcriptsBy(x=txdb, by="gene")
	Promoters <- promoters(Genes, upstream=proximalUpstream, downstream=0)

	# Swap in different range for overlaps
	if(!is.null(rangeColumn)){
		message("Swapping peaks...")
		gr <- swapRanges(gr, column="thick")
	}

	# Find overlaps
	message("Overlapping sequentially...")
	geneOverlaps <- findOverlaps(gr, Genes, select=multiHits)
	promoterOverlaps <- findOverlaps(gr, Promoters, select=multiHits)

	# Overlap sequentially
	geneName <- rep("novel", length(gr))
	geneName <- ifelse(!is.na(promoterOverlaps), promoterOverlaps, geneName)
	geneName <- ifelse(!is.na(geneOverlaps), geneOverlaps, geneName)

	# Rename novel genes to unique names
	novelGenes <- geneName == "novel"
	geneName[novelGenes] <- paste0("Novel", seq_len(sum(novelGenes)))

	# Return
	geneName
}
