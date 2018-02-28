#### txType ####

#' Annotate ranges with transcript type.
#'
#' Annotate a set of ranges in a GRanges object with transcript type (TxType)
#' based on their genic context. Transcripts are obtained from a TxDb object,
#' but can alternatively be specified manually as a GRangesList.
#'
#' @param object GRanges or RangedSummarizedExperiment: Ranges to be annotated.
#' @param txModels TxDb or GRangesList: Transcript models via a TxDb, or
#'   manually specified as a GRangesList.
#' @param outputColumn character: Name of column to hold txType values.
#' @param swap character or NULL: If not NULL, use another set of ranges
#'   contained in object to calculate overlaps, for example peaks in the thick
#'   column.
#' @param tssUpstream integer: Distance to extend annotated promoter upstream.
#' @param tssDownstream integer: Distance to extend annotated promoter
#'   downstream.
#' @param proximalUpstream integer: Maximum distance upstream of promoter to be
#'   considered proximal.
#' @param detailedAntisense logical: Wether to mirror all txType categories in
#'   the antisense direction (TRUE) or lump them all together (FALSE).
#' @param noOverlap character: In case categories are manually supplied with as
#'   a GRangesList, what to call regions with no overlap.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with txType added as factor in metadata columns (inside
#'   rowRanges in case object is a RangedSummarizedExperiment)
#'
#' @family Annotation functions
#' @export
#' @examples
#' \dontrun{
#' data(exampleUnidirectional)
#'
#' # Obtain transcript models from a TxDb-object:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'
#' # Assign txIDs
#' assignTxType(exampleUnidirectional,
#'              txModels=txdb)
#'
#' # Assign txIDs using only TC peaks:
#' exampleUnidirectional <- assignTxType(exampleUnidirectional,
#'                                       txModels=txdb,
#'                                       swap="thick")
#'
#' # The "promoter" and "proximal" category boundaries can be changed:
#' assignTxType(exampleUnidirectional,
#'              txModels=txdb,
#'              swap="thick",
#'              tssUpstream=50,
#'              tssDownstream=50,
#'              proximalUpstream=100)
#'
#' # Annotation using complete antisense categories:
#' exampleUnidirectional <- assignTxType(exampleUnidirectional,
#'                                     txModels=txdb,
#'                                     outputColumn="txTypeExtended",
#'                                     swap="thick",
#'                                     detailedAntisense=TRUE)
#'
#' # The output is always a factor added as a column:
#' summary(rowRanges(exampleUnidirectional)$txType)
#' summary(rowRanges(exampleUnidirectional)$txTypeExtended)
#'
#' # To avoid using a TxDb-object, a GRangesList can be supplied:
#' custom_hierarchy <- GRangesList(promoters=granges(promoters(txdb)),
#'                                 exons=granges(exons(txdb)))
#' assignTxType(exampleUnidirectional,
#'              txModels=custom_hierarchy,
#'              outputColumn="customType",
#'              swap="thick",
#'              noOverlap = "intergenic")
#' }
setGeneric("assignTxType", function(object, txModels, ...) {
	standardGeneric("assignTxType")
})

#' @rdname assignTxType
setMethod("assignTxType", signature(object="GenomicRanges", txModels="GenomicRangesList"), function(object, txModels, outputColumn="txType", swap=NULL, noOverlap="intergenic"){
	# Pre-checks
	assert_that(!is.null(names(txModels)),
							is.string(outputColumn),
							is.string(noOverlap),
							identical(seqlengths(object), seqlengths(txModels)))

	# Warnings
	if(outputColumn %in% colnames(mcols(object))){
		warning("object already has a column named ", outputColumn," in mcols: It will be overwritten!")
	}

	# Find overlaps
	if(is.null(swap)){
		message("Finding hierachical overlaps...")
		hits <- findOverlaps(object, txModels)
	}else{
		message("Finding hierachical overlaps with swapped ranges...")
		hits <- findOverlaps(swapRanges(object, inputColumn=swap), txModels)
	}

	# Choose highest in hierachy
	hits <- breakTies(hits, method="first")

	# Extract corresponding categories
	hits <- methods::as(hits, "List")
	hits <- extractList(names(txModels), hits)
	hits <- as.character(hits)

	# Replace NAs and format as factor
	hits <- ifelse(is.na(hits), noOverlap, hits)
	hits <- factor(hits, levels=c(names(txModels), noOverlap))

	# Append to object
	 mcols(object)[,outputColumn] <- hits

	 # Summarise annotation
	 s <- as.data.frame(table(hits))
	 colnames(s) <- c("txType", "count")
	 s$percentage <- round((s$count / sum(s$count)) * 100, digits=1)
	 s <- paste(utils::capture.output(print(s)), collapse = "\n")
	 message("### Overlap summary: ###")
	 message(s)

	# Return
	object
})

#' @rdname assignTxType
setMethod("assignTxType", signature(object="RangedSummarizedExperiment", txModels="GenomicRangesList"), function(object, txModels, ...){
	rowRanges(object) <- assignTxType(rowRanges(object), txModels=txModels, ...)

	# Return
	object
})

#' @rdname assignTxType
setMethod("assignTxType", signature(object="GenomicRanges", txModels="TxDb"), function(object, txModels, outputColumn="txType", swap=NULL, tssUpstream=100, tssDownstream=100, proximalUpstream=1000, detailedAntisense=FALSE){
	# Pre-checks
	assert_that(is.count(tssUpstream),
							is.count(tssDownstream),
							is.count(proximalUpstream),
							is.flag(detailedAntisense),
							identical(seqlengths(object), seqlengths(txModels)))

	# Build sense hierachy
	hierachy <- GRangesList(promoter=granges(trim(promoters(txModels, upstream=tssUpstream, downstream=tssDownstream))),
									 proximal=granges(trim(promoters(txModels, upstream=proximalUpstream, downstream=0))),
									 fiveUTR=granges(unlist(fiveUTRsByTranscript(txModels))),
									 threeUTR=granges(unlist(threeUTRsByTranscript(txModels))),
									 CDS=granges(cds(txModels)),
									 exon=granges(exons(txModels)),
									 intron=granges(unlist(intronsByTranscript(txModels))))

	# Build sense hierachy
	# message("Extracting txType categories...")
	# hierachy <- List(promoter=trim(promoters(txModels, upstream=tssUpstream, downstream=tssDownstream)),
	# 								 proximal=trim(promoters(txModels, upstream=proximalUpstream, downstream=0)),
	# 								 fiveUTR=fiveUTRsByTranscript(txModels),
	# 								 threeUTR=threeUTRsByTranscript(txModels),
	# 								 CDS=cds(txModels),
	# 								 exon=exons(txModels),
	# 								 intron=intronsByTranscript(txModels))
	#
	# # Coerce to GRangesList
	# #message("Coercing to GRangesList...")
	# hierachy <- lapply(hierachy, unlist)
	# hierachy <- lapply(hierachy, granges)
	# hierachy <- GRangesList(hierachy)

	# Build antisense hierachy
	#message("Adding antisense categories...")
	if(detailedAntisense){
		antisense <- invertStrand(hierachy)
		names(antisense) <- paste0("antisense_", names(antisense))
		hierachy <- c(hierachy, antisense)
	}else if(!detailedAntisense){
		antisense <- invertStrand(granges(transcripts(txModels)))
		hierachy$antisense <- antisense
	}else{
		stop("detailedAntisense must be either TRUE/FALSE!")
	}
	rm(antisense)

	# Overlap
	object <- assignTxType(object=object, txModels=hierachy, outputColumn=outputColumn, swap=swap, noOverlap="intergenic")

	# Return
	object
})

#' @rdname assignTxType
setMethod("assignTxType", signature(object="RangedSummarizedExperiment", txModels="TxDb"), function(object, txModels, ...){
	rowRanges(object) <- assignTxType(rowRanges(object), txModels=txModels, ...)

	# Return
	object
})

#### geneID ####

#' Annotate ranges with gene ID.
#'
#' Annotate a set of ranges in a GRanges object with gene IDs (i.e. Entrez Gene
#' Identifiers) based on their genic context. Features overlapping multiple
#' genes are resolved by distance to the nearest TSS. Genes are obtained from a
#' TxDb object, or can manually supplied as a GRanges.
#'
#' @param object GRanges or RangedSummarizedExperiment: Ranges to be annotated.
#' @param geneModels TxDb or GRanges: Gene models via a TxDb, or manually
#'   specified as a GRangesList.
#' @param outputColumn character: Name of column to hold geneID values.
#' @param swap character or NULL: If not NULL, use another set of ranges
#'   contained in object to calculate overlaps, for example peaks in the thick
#'   column.
#' @param upstream integer: Distance to extend annotated promoter upstream.
#' @param downstream integer: Distance to extend annotated promoter downstream.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with geneID added in columns outputColumn (inside rowRanges
#'   in case object is a RangedSummarizedExperiment)
#'
#' @family Annotation functions
#' @export
#' @examples
#' data(exampleUnidirectional)
#'
#' # Obtain gene models from a TxDb-object:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'
#' # Assign geneIDs
#' assignGeneID(exampleUnidirectional,
#'              geneModels=txdb,
#'              outputColumn="geneID")
#'
#' # Assign geneIDs using only TC peaks:
#' assignGeneID(exampleUnidirectional,
#'              geneModels=txdb,
#'              outputColumn="geneID",
#'              swap="thick")
setGeneric("assignGeneID", function(object, geneModels, ...) {
	standardGeneric("assignGeneID")
})

#' @rdname assignGeneID
setMethod("assignGeneID", signature(object="GenomicRanges", geneModels="GenomicRanges"), function(object, geneModels, outputColumn="geneID", swap=NULL, upstream=1000, downstream=100){
	# Pre-checks
	assert_that(!is.null(names(geneModels)),
							is.string(outputColumn),
							is.number(upstream),
							is.number(downstream),
							identical(seqlengths(object), seqlengths(geneModels)))

	# Warnings
	if(outputColumn %in% colnames(mcols(object))){
		warning("object already has a column named ", outputColumn," in mcols: It will be overwritten!")
	}

	# Extract anchor points
	message("Overlapping while taking distance to nearest TSS into account...")
	extendedGeneModels <- punion(geneModels, promoters(geneModels, upstream=upstream, downstream=downstream))
	extendedGeneModels <- trim(extendedGeneModels)
	TSSs <- resize(geneModels, width=1, fix="start")

	# Find overlaps
	if(is.null(swap)){
		message("Finding hierachical overlaps...")
		hits <- findOverlaps(object, extendedGeneModels)
	}else{
		message("Finding hierachical overlaps with swapped ranges...")
		hits <- findOverlaps(swapRanges(object, inputColumn=swap), extendedGeneModels)
	}

	# Calculate distances to TSSs
	mcols(hits)$distance <- distance(Pairs(object, TSSs, hits=hits))
	rm(TSSs, extendedGeneModels)

	# Resolve by distance to nearest TSS by split apply, THIS CAN BE UPDATED IN NEXT VERSION OF S4Vectors!
	hits <- hits[which.min(splitAsList(mcols(hits)$distance, queryHits(hits)), global=TRUE)]

	# Extract ids
	hits <- methods::as(hits, "List")
	hits <- extractList(names(geneModels), hits)
	hits <- as.character(hits)
	stopifnot(length(hits) == length(object))

	# Append to object
	mcols(object)[,outputColumn] <- hits

	# Return
	message("### Overlap Summary: ###")
	message("Features overlapping genes: ", round(mean(!is.na(hits)) * 100, digits=2), " %")
	message("Number of unique genes: ", length(unique(hits)))

	# Return
	object
})

#' @rdname assignGeneID
setMethod("assignGeneID", signature(object="RangedSummarizedExperiment", geneModels="GenomicRanges"), function(object, geneModels, ...){
	rowRanges(object) <- assignGeneID(rowRanges(object), geneModels=geneModels, ...)

	# Return
	object
})

#' @rdname assignGeneID
setMethod("assignGeneID", signature(object="GenomicRanges", geneModels="TxDb"), function(object, geneModels, outputColumn="geneID", swap=NULL, upstream=1000, downstream=100){
	# Pre-checks
	assert_that(is.string(outputColumn),
							is.number(upstream),
							is.number(downstream),
							identical(seqlengths(object), seqlengths(geneModels)))

	# Extract anchor points
	message("Extracting genes...")
	as_gr <- genes(geneModels)

	# Overlap
	object <- assignGeneID(object, geneModels=as_gr, outputColumn=outputColumn, swap=swap, upstream=upstream, downstream=downstream)

	# Return
	object
})

#' @rdname assignGeneID
setMethod("assignGeneID", signature(object="RangedSummarizedExperiment", geneModels="TxDb"), function(object, geneModels, ...){
	rowRanges(object) <- assignGeneID(rowRanges(object), geneModels=geneModels, ...)

	# Return
	object
})

#### txID ####

#' Annotate ranges with transcript ID.
#'
#' Annotate a set of ranges in a GRanges object with transcript IDs based on
#' their genic context. All overlapping transcripts are returned. Transcripts
#' are obtained from a TxDb object, or can manually supplied as a GRanges.
#'
#' @param object GRanges or RangedSummarizedExperiment: Ranges to be annotated.
#' @param txModels TxDb or GRanges: Transcript models via a TxDb, or manually
#'   specified as a GRanges.
#' @param outputColumn character: Name of column to hold txID values.
#' @param swap character or NULL: If not NULL, use another set of ranges
#'   contained in object to calculate overlaps, for example peaks in the thick
#'   column.
#' @param upstream integer: Distance to extend annotated promoter upstream.
#' @param downstream integer: Distance to extend annotated promoter downstream.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with txID added in column outputColumn (inside rowRanges in
#'   case object is a RangedSummarizedExperiment)
#'
#' @family Annotation functions
#' @export
#' @examples
#' data(exampleUnidirectional)
#'
#' # Obtain transcript models from a TxDb-object:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'
#' # Assign txIDs
#' assignTxID(exampleUnidirectional,
#'            txModels=txdb,
#'            outputColumn="geneID")
#'
#' # Assign txIDs using only TC peaks:
#' assignTxID(exampleUnidirectional,
#'              txModels=txdb,
#'              outputColumn="geneID",
#'              swap="thick")
setGeneric("assignTxID", function(object, txModels, ...) {
	standardGeneric("assignTxID")
})

#' @rdname assignTxID
setMethod("assignTxID", signature(object="GenomicRanges", txModels="GenomicRanges"), function(object, txModels, outputColumn="txID", swap=NULL){
	# Pre-checks
	assert_that(!is.null(names(txModels)),
							is.string(outputColumn),
							identical(seqlengths(object), seqlengths(txModels)))

	# Warnings
	if(outputColumn %in% colnames(mcols(object))){
		warning("object already has a column named ", outputColumn," in mcols: It will be overwritten!")
	}

	# Find overlaps
	if(is.null(swap)){
		message("Finding hierachical overlaps...")
		hits <- findOverlaps(object, txModels)
	}else{
		message("Finding hierachical overlaps with swapped ranges...")
		hits <- findOverlaps(swapRanges(object, inputColumn=swap), txModels)
	}

	# Count unique transcrips for later
	nTxs <- length(unique(to(hits)))

	# Extract names
	hits <- methods::as(hits, "List")
	hits <- extractList(names(txModels), hits)
	hits <- paste(hits, collapse=";")
	hits <- ifelse(hits == "", NA, hits)

	# Checks
	stopifnot(length(hits) == length(object),
						is.character(hits))

	# Append to object
	mcols(object)[,outputColumn] <- hits

	# Return
	message("### Overlap Summary: ###")
	message("Features overlapping transcripts: ", round(mean(!is.na(hits)) * 100, digits=2), " %")
	message("Number of unique transcripts: ", nTxs)

	# Return
	object
})

#' @rdname assignTxID
setMethod("assignTxID", signature(object="RangedSummarizedExperiment", txModels="GenomicRanges"), function(object, txModels, ...){
	rowRanges(object) <- assignTxID(rowRanges(object), txModels=txModels, ...)

	# Return
	object
})

#' @rdname assignTxID
setMethod("assignTxID", signature(object="GenomicRanges", txModels="TxDb"), function(object, txModels, outputColumn="txID", swap=NULL, upstream=1000, downstream=0){
	# Pre-checks
	assert_that(is.string(outputColumn),
							is.number(upstream),
							upstream >= 0,
							is.number(downstream),
							downstream >= 0,
							identical(seqlengths(object), seqlengths(txModels)))

	# Extract anchor points
	message("Extracting transcripts...")
	txs <- transcripts(txModels, columns="tx_name")
	names(txs) <- txs$tx_name
	txs <- punion(txs, promoters(txs, upstream=upstream, downstream=downstream))
	txs <- trim(txs)

	# Overlap
	object <- assignTxID(object, txModels=txs, outputColumn=outputColumn)

	# Return
	object
})

#' @rdname assignTxID
setMethod("assignTxID", signature(object="RangedSummarizedExperiment", txModels="TxDb"), function(object, txModels, ...){
	rowRanges(object) <- assignTxID(rowRanges(object), txModels=txModels, ...)

	# Return
	object
})

#### NAs ####

#' Annotate ranges with missing IDs.
#'
#' This function can relabel ranges with missing IDs (i.e. returned by
#' assignTxID and assignGeneID), in case they need to be retained for further
#' analysis.
#'
#' @param object character, GRanges or RangedSummarizedExperiment: IDs to have
#'   NAs replaces with new IDs.
#' @param outputColumn character: Name of column to hold txID values.
#' @param prefix character: New name to assign to ranges with missing IDs, in
#'   the style prefix1, prefix2, etc.
#' @param ... additional arguments passed to methods.
#'
#' @return GRanges with NAs replaced in outputColumn (inside rowRanges in case
#'   object is a RangedSummarizedExperiment)
#'
#' @family Annotation functions
#' @export
#' @examples
#' data(exampleUnidirectional)
#'
#' # Obtain gene models from a TxDb-object:
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'
#' # Assign geneIDs using only TC peaks:
#' exampleUnidirectional <- assignGeneID(exampleUnidirectional,
#'                                       geneModels=txdb,
#'                                       outputColumn="geneID",
#'                                       swap="thick")
#'
#' # Replace NAs with "Novel"
#' assignMissingID(exampleUnidirectional)
#'
#' # Replace NAs with "NovelTSS"
#' assignMissingID(exampleUnidirectional, prefix = "NovelTSS")
setGeneric("assignMissingID", function(object, ...) {
	standardGeneric("assignMissingID")
})

#' @rdname assignMissingID
setMethod("assignMissingID", signature(object="character"), function(object, prefix="Novel"){
	# Pre-checks
	assert_that(is.string(prefix))

	missingIDs <- is.na(object)
	totalMissing <- sum(missingIDs)
	object[missingIDs] <- paste0(prefix, seq_len(totalMissing))

	message("Assigned ", totalMissing, " missing IDs")
	object
})

#' @rdname assignMissingID
setMethod("assignMissingID", signature(object="GenomicRanges"), function(object, outputColumn="geneID", prefix="Novel"){
	# Pre-checks
	assert_that(is.string(outputColumn),
							outputColumn %in% colnames(mcols(object)))

	# Replace
	mcols(object)[,outputColumn] <- assignMissingID(mcols(object)[,outputColumn], prefix=prefix)

	# Return
	object
})

#' @rdname assignMissingID
setMethod("assignMissingID", signature(object="RangedSummarizedExperiment"), function(object, outputColumn="geneID", prefix="Novel"){
	rowRanges(object) <- assignMissingID(rowRanges(object), outputColumn=outputColumn, prefix=prefix)

	# Return
	object
})
