#' Annotate GRanges with transcript type
#'
#' Annotate a set of ranges in a GRanges object with transcript type (TxType) based on their genic context. Transcripts are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation.
#' @param tssUpstream integer: Distance to extend annotated promoter upstream.
#' @param tssDownstream integer: Distance to extend annotated promoter downstream.
#' @param proximalUpstream integer: Maximum distance upstream of promoter to be considered proximal.
#'
#' @return character vector of same length as gr containing TxType for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignTxType <- function(gr, txdb, tssUpstream=100, tssDownstream=100, proximalUpstream=1000){
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

#' Annotate GRanges with transcript ID
#'
#' Annotate a set of ranges in a GRanges object with Entrez Transcript Identifier (ID) based on their genic context. Genes are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation
#' @param proximalUpstream integer: Maximum distance upstream of annotated promoter to be considered part of genes.
#' @param proximalDownstream integer: Maximum distance downstream of annotated promoter to be considered part of genes.
#' @param multiHits character: How to choose when multiple genes overlap.
#' @return character vector of same length as gr containing Entrez transcript IDs for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignTxID <- function(gr, txdb, proximalUpstream=1000, proximalDownstream=100, multiHits="arbitrary"){
	# TO DO
	# Option for reporting more than one gene.
	# Option to discard multi strand genes.
	# Option to not swap ranges
	# Pre/post checks

	# Reference categories
	message("Extracting transcripts...")
	Transcripts <- transcripts(txdb)
	Promoters <- promoters(txdb, upstream=proximalUpstream, downstream=proximalDownstream)

	# Find overlaps
	message("Overlapping sequentially...")
	transcriptOverlaps <- findOverlaps(gr, Transcripts, select=multiHits)
	promoterOverlaps <- findOverlaps(gr, Promoters, select=multiHits)

	# Overlap sequentially
	transcriptName <- rep(NA, length(gr))
	transcriptName <- ifelse(!is.na(promoterOverlaps), Promoters$tx_name[promoterOverlaps], transcriptName)
	transcriptName <- ifelse(!is.na(transcriptOverlaps), Transcripts$tx_name[transcriptOverlaps], transcriptName)

	# Return
	transcriptName
}

#' Annotate GRanges with gene ID
#'
#' Annotate a set of ranges in a GRanges object with Entrez Gene Identifier (ID) based on their genic context. Genes are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation
#' @param proximalUpstream integer: Maximum distance upstream of annotated promoter to be considered part of genes.
#' @param proximalDownstream integer: Maximum distance downstream of annotated promoter to be considered part of genes.
#' @param multiHits character: How to choose when multiple genes overlap.
#' @return character vector of same length as gr containing Entrez gene IDs for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignGeneID <- function(gr, txdb, proximalUpstream=1000, proximalDownstream=100, multiHits="arbitrary"){
	# TO DO
	# Option for reporting more than one gene.
	# Option to discard multi strand genes.
	# Option to not swap ranges
	# Pre/post checks

	# Reference categories
	message("Extracting genes...")
	#Genes <- genes(txdb, single.strand.genes.only=TRUE)
	Genes <- transcriptsBy(x=txdb, by="gene")
	Promoters <- promoters(Genes, upstream=proximalUpstream, downstream=proximalDownstream)

	# Find overlaps
	message("Overlapping sequentially...")
	geneOverlaps <- findOverlaps(gr, Genes, select=multiHits)
	promoterOverlaps <- findOverlaps(gr, Promoters, select=multiHits)

	# Overlap sequentially
	geneName <- rep(NA, length(gr))
	geneName <- ifelse(!is.na(promoterOverlaps), names(Promoters)[promoterOverlaps], geneName)
	geneName <- ifelse(!is.na(geneOverlaps), names(Genes)[geneOverlaps], geneName)

	# Return
	geneName
}

#' Collapse an Expression Matrix by summin over genes
#'
#' Summarise a transcript-level EM to a gene-level EM, i.e. for GO-term analysis.
#'
#' @param em matrix: Expression matrix.
#' @param genes factor or character: Vector of gene names to group by.
#' @param keepUnannotated logical: Whether to discard (FALSE) or keep (TRUE) un annotated transcripts.
#' @param prefix character: If keepUnannotated is TRUE, this prefix is added to all unannotated transcripts along with counter.
#'
#' @return Expression matrix summed by genes.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import data.table
#' @export
sumByGene <- function(em, genes, keepUnannotated=FALSE, prefix="Novel"){
	# Build data.table
	d <- data.table(gene=genes, em)

	# Either remove or relabel NAs
	if(keepUnannotated==FALSE){
		d <- d[!is.na(gene)]
	}else if(keepUnannotated==TRUE){
		d[,gene := ifelse(is.na(gene),
											paste0(prefix, seq_len(sum(is.na(gene)))),
											gene)]
	}

	# Sum by gene
	d <- d[, lapply(.SD, sum), by="gene"]

	# Data.table to matrix
	rowgene <- d$gene
	d[,gene := NULL]
	m <- as.matrix(d)
	rownames(m) <- rowgene

	# Return
	m
}
