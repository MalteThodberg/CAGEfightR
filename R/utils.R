splitByStrand <- function(object){
	split(object, strand(object))
}

extendRanges <- function(object, upstream, downstream){
	U <- flank(object, width=upstream, start=TRUE)
	D <- flank(object, width=downstream, start=FALSE)
	object <- punion(object, U)
	object <- punion(object, D)
	object
}

countScoredOverlaps <- function(query, subject){
	# Solution from bioconductor by M Morgan
	hits <- methods::as(findOverlaps(query=query, subject=subject), "List")
	weightedCount <- sum(extractList(score(subject), hits))

	# Return
	weightedCount
}
