#' Stranded Coverage Genome Browser Track
#'
#' Create a Gviz-track with stranded coverage (Plus/minus coverage as positive/negative). This representation makes it easy to identify bidirectional peaks.
#'
#' @param gr GRanges: GRanges with coverage in the score column.
#' @param plusColor character: Color for plus-strand coverage.
#' @param minusColor character: Color for minus-strand coverage.
#' @param ... additional names passed on to DataTrack.
#'
#' @return DataTrack-object.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Genome Browser functions
#' @import S4Vectors IRanges GenomicRanges Gviz
#' @export
strandedCoverageTrack <- function(gr, plusColor="cornflowerblue", minusColor="tomato", ...){
	# Vector by strand
	message("Splitting coverage by strand...")
	by_strand <- splitByStrand(gr)
	plus_coverage <- coverage(by_strand$`+`, weight="score")
	minus_coverage <- 0 - coverage(by_strand$`-`, weight="score")

	# Back to GRanges
	message("Preparing track...")
	names(minus_coverage) <- names(plus_coverage)
	o <- bindAsGRanges(plus=plus_coverage, minus=minus_coverage)

	# Build track
	o <- DataTrack(o, type="histogram", groups=c("plus", "minus"), col=c(minusColor, plusColor), ...)

	# Return
	o
}

#' Tag Cluster Genome Browser Track
#'
#' Create a Gviz-track of Tag Clusters (TCs), indicating TC peak with thickness and strand with color.
#'
#' @param gr GRanges: GRanges with peaks in the thick-column.
#' @param plusColor character: Color for plus-strand coverage.
#' @param minusColor character: Color for minus-strand coverage.
#' @param ... additional names passed on to GeneRegionTrack.
#'
#' @return GeneRegionTrack-object.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Genome Browser functions
#' @import S4Vectors IRanges GenomicRanges Gviz
#' @export
tagClusterTrack <- function(gr, plusColor="cornflowerblue", minusColor="tomato", ...){
	# Extract peaks
	message("Setting thick and thin features...")
	insideThick <- swapRanges(gr)

	# Remove mcols and add features for thin feature
	names(insideThick) <- NULL
	mcols(insideThick) <- NULL
	insideThick$feature <- ifelse(strand(insideThick) == "+", "thickPlus", "thickMinus")

	# Remove peaks from TCs
	outsideThick <- setdiff(gr, insideThick)

	# Remove mcols and add features
	mcols(outsideThick) <- NULL
	outsideThick$feature <- ifelse(strand(outsideThick) == "+", "thinPlus", "thinMinus")

	# Temporary to GRangesList for easy sorting
	message("Merging and sorting...")
	o <- sort(c(insideThick, outsideThick))
	fo <- findOverlaps(o, gr, select="arbitrary")
	o <- split(o, fo)
	names(o) <- names(gr)
	o <- unlist(o)

	# Add necessary columns for track
	o$transcript <- names(o)
	o$gene <- o$transcript
	o$symbol <- o$transcript

	# Build track
	message("Preparing track...")
	o <- GeneRegionTrack(o, thinBoxFeature=c("thinPlus", "thinMinus"),
											 min.distance=0, collapse=FALSE,
											 thinPlus=plusColor, thickPlus=plusColor,
											 thinMinus=minusColor, thickMinus=minusColor,
											 ...)

	# Return
	o
}
