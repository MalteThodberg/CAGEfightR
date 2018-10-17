splitByStrand <- function(object) {
	split(object, strand(object))
}

splitPooled <- function(object){
	# Checks never done inside...

	# Split by strand
	o <- splitByStrand(object)

	# Calculate coverage
	o <- lapply(o, coverage, weight="score")

	# Round to handle floating point errors
	o <- lapply(o, round, digits=getOption("CAGEfightR.round"))

	# Return
	o
}
