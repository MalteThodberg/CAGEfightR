#' Combine two CAGE experiments.
#'
#' This function can safely combine two CAGE experiments, for example TCs and
#' enhancers, for later analysis, by making sure no ranges in the final object
#' are overlapping.
#'
#' @param object1 RangedSummarizedExperiment: First experiment to be combined.
#' @param object2 RangedSummarizedExperiment: First experiment to be combined.
#' @param removeIfOverlapping character: Whether to keep overlapping ranges
#'   ("none") or discard from either the first ("object1") or second ("object2")
#'   experiment.
#' @param ... arguments passed to methods.
#'
#' @return RangedSummarizedExperiment with merged and sorted ranges (colData and
#'   metadata are carried over unchanged).
#' @export
#' @examples
#' data(exampleUnidirectional)
#' data(exampleBidirectional)
#'
#' # Clusters must have identical colData to be combined:
#' exampleUnidirectional$totalTags <- NULL
#'
#' # Combine, keeping potential overlaps
#' combineClusters(object1=exampleUnidirectional, object2=exampleBidirectional)
#'
#' # If features overlap, keep only from object1
#' combineClusters(object1=exampleUnidirectional, object2=exampleBidirectional,
#'    removeIfOverlapping="object2")
#'
#' # If features overlap, keep only from object2
#' combineClusters(object1=exampleUnidirectional, object2=exampleBidirectional,
#'    removeIfOverlapping="object1")
setGeneric("combineClusters", function(object1, object2, ...) {
	standardGeneric("combineClusters")
})

#' @rdname combineClusters
#' @export
setMethod("combineClusters",
					signature(object1="RangedSummarizedExperiment",
										"RangedSummarizedExperiment"),
					function(object1, object2, removeIfOverlapping="none"){
	# Pre-checks
	assert_that(!is.null(rownames(object1)),
							!is.null(rownames(object2)),
							identical(seqinfo(object1), seqinfo(object2)),
							identical(colData(object1), colData(object2)),
							removeIfOverlapping %in% c("object1", "object2", "none"))

	# Remove overlapping features
	if(removeIfOverlapping == "object1"){
		# Find overlaps
		os <- overlapsAny(object1, object2)

		# Remove
		message("Removing overlapping features from object1: ", sum(os))
		object1 <- subset(object1, !os)
	}else if(removeIfOverlapping == "object2"){
		# Find overlaps
		os <- overlapsAny(object2, object1)

		# Remove
		message("Removing overlapping features from object2: ", sum(os))
		object2 <- subset(object2, !os)
	}else{
		os <- intersect(rowRanges(object1), rowRanges(object2))
		message("Overlapping features were not removed: ", length(os))
	}
	rm(os)

	# Only retain shared assays
	shared_assays <- intersect(assayNames(object1), assayNames(object2))
	message("Keeping assays: ", paste(shared_assays, collapse=", "))
	assays(object1) <- assays(object1)[shared_assays]
	assays(object2) <- assays(object2)[shared_assays]
	rm(shared_assays)

	# Only retain shared rowcolumns
	shared_cols <- intersect(colnames(mcols(object1)), colnames(mcols(object2)))
	message("Keeping columns: ", paste(shared_cols, collapse=", "))
	mcols(object1) <- mcols(object1)[,shared_cols]
	mcols(object2) <- mcols(object2)[,shared_cols]
	rm(shared_cols)

	# Attempt to merge metadata
	message("Merging metadata...")
	new_meta <- c(metadata(object1), metadata(object2))
	metadata(object1) <- list()
	metadata(object2) <- list()
	new_meta <- new_meta[!duplicated(new_meta)]

	# Bind, add meta and re-sort
	message("Stacking and re-sorting...")
	o <- rbind(object1, object2)
	metadata(o) <- new_meta
	rm(new_meta)
	o <- sort(o)

	# Return
	o
})
