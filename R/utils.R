#' Split by strand.
#'
#' Splits an object into a list by strand
#'
#' @param object object with strand information (GRanges, Gpos, RangedSummarizedExperiment, etc.).
#'
#' @return List with three elements: +, - and *.
#' @family Helper functions
#' @import GenomicRanges SummarizedExperiment
#' @export
splitByStrand <- function(object){
	split(object, strand(object))
}

#' @import GenomicRanges SummarizedExperiment
extendRanges <- function(object, upstream, downstream){
	U <- flank(object, width=upstream, start=TRUE)
	D <- flank(object, width=downstream, start=FALSE)
	object <- punion(object, U)
	object <- punion(object, D)
	object
}
