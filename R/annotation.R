#' Annotate GRanges with transcript type
#'
#' Annotate a set of ranges in a GRanges object with transcript type (TxType) based on their genic context. Transcripts are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation.
#' @param tssUpstream integer: Distance to extend annotated promoter upstream.
#' @param tssDownstream integer: Distance to extend annotated promoter downstream.
#' @param proximalUpstream integer: Maximum distance upstream of promoter to be considered proximal.
#' @param asFactor logical: Whether to output a character (FALSE) or a factor (TRUE).
#'
#' @return character vector of same length as gr containing TxType for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignTxType <- function(gr, txdb, tssUpstream=100, tssDownstream=100, proximalUpstream=1000, asFactor=FALSE){
	# TO DO
	# Add arguments for promoter and upstream coordinates.
	# Option to not swap ranges
	# Pre/post checks

	# Reference categories
	message("Extracting reference categories...")
	Promoters <- trim(promoters(txdb, upstream=tssUpstream, downstream=tssDownstream))
	Proximal <- trim(promoters(txdb, upstream=proximalUpstream, downstream=0))
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

	if(asFactor){
		featureType <- factor(featureType,
													levels=c("intergenic",
																	 "proximal", "promoter", "fiveUTR",
																	 "CDS", "intron", "exon", "threeUTR",
																	 "antisense"))
	}

	# Return
	featureType
}

#' Simplify transcript types
#'
#' Simplify annotation of transcript types (txTypes).
#'
#' @param txTypes factor: txTypes returned by asignTxTypes-function.
#' @param scheme character or list: Scheme for simplification: Currently supports "genic" and "promoter proximity. Advanced users may directly supply a named list of levels.
#'
#' @return factor with simplified txTypes
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
simplifyTxTypes <- function(txTypes, scheme="genic"){
	stopifnot(is.factor(txTypes))

	if(scheme == "genic"){
		levels(txTypes) <- list(intergenic=c("proximal", "intergenic"),
													 intragenic=c("promoter", "fiveUTR", "CDS", "intron", "exon", "threeUTR"),
													 antisense=c("antisense"))
	}else if(scheme == "promoter proximity"){
		levels(txTypes) <- list(promoter=c("promoter"),
													 adjacent=c("proximal", "fiveUTR"),
													 intragenic=c("CDS", "intron", "exon", "threeUTR"),
													 intergenic=c("intergenic"),
													 antisense=c("antisense"))
	}else{
		message("Using supplied conversion scheme...")
		levels(txTypes) <- scheme
	}

	if(anyNA(txTypes)){
		warning("Output factor has missing values!")
	}

	# Return
	txTypes
}

#' Relabel unannotated genes
#'
#' Relabels unannotated genes, i.e. as "Novel1". Usesful when testing for Alternative Transcript Usage (ATU).
#'
#' @param geneIDs character: Gene names, with NA indicating unannotated.
#' @param prefix character: Prefix to add before counter.
#'
#' @return character vector of the same length as geneID
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
relabelUnannotated <- function(geneIDs, prefix="Novel"){
	stopifnot(is.character(geneIDs))

	# Find unnannotated
	unannotated <- is.na(geneIDs)

	# Make new names
	new_names <- paste0(prefix, 1:sum(unannotated))

	# Replace
	geneIDs[unannotated] <- new_names

	# Return
	geneIDs
}

#' Annotate GRanges with transcript ID
#'
#' Annotate a set of ranges in a GRanges object with Transcript Identifier (TxID) based on their genic context. Transcripts are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation
#' @param upstream integer: Maximum distance upstream of annotated TSS to be considered part of genes.
#' @param downstream integer: Maximum distance downstream of annotated TSS to be considered part of genes.
#' @return character vector of same length as gr containing transcript names.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignTxID <- function(gr, txdb, upstream=1000, downstream=100){
	# TO DO
	# Option for reporting more than one gene.
	# Option to discard multi strand genes.
	# Option to not swap ranges
	# Pre/post checks

	# TO DO
	# Add arguments for promoter and upstream coordinates.
	# Option to not swap ranges
	# Pre/post checks

	# Reference categories
	message("Extracting transcripts...")
	Tx <- transcripts(txdb)
	TxExpanded <- punion(promoters(Tx, upstream=1000, downstream=100), Tx)
	TxExpanded <- trim(TxExpanded)
	names(TxExpanded) <- Tx$tx_name

	# Overlaps
	message("Overlapping...")
	hits <- findOverlaps(query=gr, subject=TxExpanded, select="all")

	# Extract names
	hits <- methods::as(hits, "List")
	TxNames <- extractList(names(TxExpanded), hits)
	TxNames <- paste(TxNames, collapse=";")
	TxNames <- ifelse(TxNames == "", NA, TxNames)

	# Checks
	stopifnot(length(TxNames) == length(gr),
						is.character(TxNames))

	# Return
	TxNames
}

#' Annotate GRanges with gene ID
#'
#' Annotate a set of ranges in a GRanges object with gene IDs (i.e. Entrez Gene Identifiers) based on their genic context. Genes are obtained from a TxDb object.
#'
#' @param gr GRanges: Ranges to annotated.
#' @param txdb TxDb: Transcript database to use for annotation.
#' @param upstream integer: Region upstream of gene to be considered part of the gene.
#' @param downstream integer: Region downstream of gene to be considered part of the gene.
#' @param resolveRule character: How to resolve multiple overlaps: "first", "last", "arbitrary", "shortest" or "longest".
#'
#' @return character vector of same length as gr containing a single gene ID for each range.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' @export
assignGeneID <- function(gr, txdb, upstream=1000, downstream=0, resolveRule="shortest"){
	# Extract and expand
	message("Extracting genes...")
	Genes <- genes(txdb)
	Genes <- extendRanges(Genes, upstream=upstream, downstream=downstream)
	Genes <- trim(Genes)

	# Overlap
	message("Overlapping...")
	o <- resolvedOverlaps(query=gr, subject=Genes, select=resolveRule)

	# Translate to gene names
	o <- ifelse(is.na(o), NA, names(Genes)[o])

	# Check
	stopifnot(length(o) == length(gr),
						is.character(o))

	# Return
	o
}

#' Collapse Expression Matrix and Ranges by summing over genes
#'
#' Summarize a transcript-level EM to a gene-level EM,  i.e. for GO-term analysis. Transcripts are grouped within genes. Unannotated transcripts (NAs) are discarded.
#'
#' @param RSE RangedSummarizedExperiment: Transcript-level data to be summarized.
#' @param inputAssay character: Name of assay holding expression values to be summed, i.e. "counts".
#' @param geneID character: Name of column in rowData holding geneIDs.
#' @param calcStats logical: Whether to calculate extra statistics for the merged transcripts (i.e. number of transcripts in each gene).
#'
#' @return RangedSummarizedExperiment with inputAssay summed over genes and rowRanges as a GRangesList. If calcStats=TRUE, rowData contains gene-wise statistics. colData and metaData are unchanged from the original RangedSummarizedExperiment.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges SummarizedExperiment
#' @export
sumOverGenes <- function(RSE, inputAssay, geneID, calcStats=TRUE){
	# Checks
	stopifnot(class(RSE) == "RangedSummarizedExperiment")
	stopifnot(geneID %in% colnames(rowData(RSE)))
	stopifnot(inputAssay %in% assayNames(RSE))

	# Split gr
	old_gr <- rowRanges(RSE)
	geneID <- mcols(old_gr)[,geneID]
	new_gr <- split(old_gr, geneID)

	# Sum matrix
	old_m <- assay(RSE, inputAssay)
	new_m <- suppressWarnings(rowsum(old_m, group=geneID, na.rm=FALSE))

	# Discard NA and check names match
	new_m <- new_m[!is.na(rownames(new_m)),]
	stopifnot(rownames(new_m) == names(new_gr))

	# Reassemble

	o <- SummarizedExperiment(assays=list(new_m),
														rowRanges=new_gr,
														colData=colData(RSE),
														metadata=metadata(RSE))
	assayNames(o) <- inputAssay

	# Calculate some stats
	if(calcStats){
		rowData(o)[,"nFeatures"] <- elementNROWS(new_gr)
	}

	# Return
	o
}



### LEGACY VERSIONS

#' #' Annotate GRanges with gene ID
#' #'
#' #' Annotate a set of ranges in a GRanges object with Entrez Gene Identifier (ID) based on their genic context. Genes are obtained from a TxDb object.
#' #'
#' #' @param gr GRanges: Ranges to annotated.
#' #' @param txdb TxDb: Transcript database to use for annotation
#' #' @param proximalUpstream integer: Maximum distance upstream of annotated promoter to be considered part of genes.
#' #' @param proximalDownstream integer: Maximum distance downstream of annotated promoter to be considered part of genes.
#' #' @param multiHits character: How to choose when multiple genes overlap.
#' #' @return character vector of same length as gr containing Entrez gene IDs for each range.
#' #' @examples
#' #' # ADD_EXAMPLES_HERE
#' #' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' #' @export
#' assignGeneID <- function(gr, txdb, proximalUpstream=1000, proximalDownstream=100, multiHits="arbitrary"){
#' 	# TO DO
#' 	# Option for reporting more than one gene.
#' 	# Option to discard multi strand genes.
#' 	# Option to not swap ranges
#' 	# Pre/post checks
#'
#' 	# Reference categories
#' 	message("Extracting genes...")
#' 	#Genes <- genes(txdb, single.strand.genes.only=TRUE)
#' 	Genes <- transcriptsBy(x=txdb, by="gene")
#' 	Promoters <- promoters(Genes, upstream=proximalUpstream, downstream=proximalDownstream)
#'
#' 	# Find overlaps
#' 	message("Overlapping sequentially...")
#' 	geneOverlaps <- findOverlaps(gr, Genes, select=multiHits)
#' 	promoterOverlaps <- findOverlaps(gr, Promoters, select=multiHits)
#'
#' 	# Overlap sequentially
#' 	geneName <- rep(NA, length(gr))
#' 	geneName <- ifelse(!is.na(promoterOverlaps), names(Promoters)[promoterOverlaps], geneName)
#' 	geneName <- ifelse(!is.na(geneOverlaps), names(Genes)[geneOverlaps], geneName)
#'
#' 	# Return
#' 	geneName
#' }

#' #' Annotate GRanges with transcript ID
#' #'
#' #' Annotate a set of ranges in a GRanges object with Entrez Transcript Identifier (ID) based on their genic context. Genes are obtained from a TxDb object.
#' #'
#' #' @param gr GRanges: Ranges to annotated.
#' #' @param txdb TxDb: Transcript database to use for annotation
#' #' @param proximalUpstream integer: Maximum distance upstream of annotated promoter to be considered part of genes.
#' #' @param proximalDownstream integer: Maximum distance downstream of annotated promoter to be considered part of genes.
#' #' @param multiHits character: How to choose when multiple genes overlap.
#' #' @return character vector of same length as gr containing Entrez transcript IDs for each range.
#' #' @examples
#' #' # ADD_EXAMPLES_HERE
#' #' @import S4Vectors IRanges GenomicRanges GenomicFeatures
#' #' @export
#' assignTxID <- function(gr, txdb, proximalUpstream=1000, proximalDownstream=100, multiHits="arbitrary"){
#' 	# TO DO
#' 	# Option for reporting more than one gene.
#' 	# Option to discard multi strand genes.
#' 	# Option to not swap ranges
#' 	# Pre/post checks
#'
#' 	# Reference categories
#' 	message("Extracting transcripts...")
#' 	Transcripts <- transcripts(txdb)
#' 	Promoters <- promoters(txdb, upstream=proximalUpstream, downstream=proximalDownstream)
#'
#' 	# Find overlaps
#' 	message("Overlapping sequentially...")
#' 	transcriptOverlaps <- findOverlaps(gr, Transcripts, select=multiHits)
#' 	promoterOverlaps <- findOverlaps(gr, Promoters, select=multiHits)
#'
#' 	# Overlap sequentially
#' 	transcriptName <- rep(NA, length(gr))
#' 	transcriptName <- ifelse(!is.na(promoterOverlaps), Promoters$tx_name[promoterOverlaps], transcriptName)
#' 	transcriptName <- ifelse(!is.na(transcriptOverlaps), Transcripts$tx_name[transcriptOverlaps], transcriptName)
#'
#' 	# Return
#' 	transcriptName
#' }

#' #' Collapse an Expression Matrix by summing over genes
#' #'
#' #' Summarise a transcript-level EM to a gene-level EM, i.e. for GO-term analysis.
#' #'
#' #' @param em matrix: Expression matrix.
#' #' @param genes factor or character: Vector of gene names to group by.
#' #' @param keepUnannotated logical: Whether to discard (FALSE) or keep (TRUE) un annotated transcripts.
#' #' @param prefix character: If keepUnannotated is TRUE, this prefix is added to all unannotated transcripts along with a counter.
#' #'
#' #' @return Expression matrix summed by genes.
#' #' @examples
#' #' # ADD_EXAMPLES_HERE
#' #' @export
#' sumByGeneLegacy <- function(em, genes, keepUnannotated=FALSE, prefix="Novel"){
#' 	# Check that input is matrix and factor
#' 	stopifnot(is.matrix(em))
#' 	stopifnot(is.character(genes))
#'
#' 	# Relabel unannotated if they should be kept
#' 	nNovel <- sum(is.na(genes))
#'
#' 	if(keepUnannotated){
#' 		message(paste0("Relabelling ", nNovel, " unannotated features..."))
#' 		genes <- ifelse(is.na(genes), paste0(prefix, seq_len(sum(is.na(genes)))), genes)
#' 	}else{
#' 		message(paste0("Discarding ", nNovel, " unannotated features..."))
#' 	}
#'
#' 	# Sum by gene (this automatically skips NAs)
#' 	o <- by(data=em, INDICES=factor(genes), FUN=colSums)
#'
#' 	# Coerce back to matrix
#' 	o <- do.call(rbind, o)
#'
#' 	# Return
#' 	o
#' }

# #' Collapse an Expression Matrix by summin over genes
# #'
# #' Summarise a transcript-level EM to a gene-level EM, i.e. for GO-term analysis.
# #'
# #' @param em matrix: Expression matrix.
# #' @param genes factor or character: Vector of gene names to group by.
# #' @param keepUnannotated logical: Whether to discard (FALSE) or keep (TRUE) un annotated transcripts.
# #' @param prefix character: If keepUnannotated is TRUE, this prefix is added to all unannotated transcripts along with counter.
# #'
# #' @return Expression matrix summed by genes.
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import data.table
# #' @export
# sumByGene <- function(em, genes, keepUnannotated=FALSE, prefix="Novel"){
# 	# Build data.table
# 	d <- data.table(gene=genes, em)
#
# 	# Either remove or relabel NAs
# 	if(keepUnannotated==FALSE){
# 		d <- d[!is.na(gene)]
# 	}else if(keepUnannotated==TRUE){
# 		d[,gene := ifelse(is.na(gene),
# 											paste0(prefix, seq_len(sum(is.na(gene)))),
# 											gene)]
# 	}
#
# 	# Sum by gene
# 	d <- d[, lapply(.SD, sum), by="gene"]
#
# 	# Data.table to matrix
# 	rowgene <- d$gene
# 	d[,gene := NULL]
# 	m <- as.matrix(d)
# 	rownames(m) <- rowgene
#
# 	# Return
# 	m
# }
