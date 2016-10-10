quantifyFeatures <- function(ctss, features) UseMethod("quantifyFeatures")

# Default method
quantifyFeatures.default <- function(ctss, features){
	stop("ctss must be either a character or GRangesList")
}

# Count overlaps with scores
countScoredOverlaps <- function(gr, features){
		# Solution from bioconductor by M Morgan
		hits <- as(findOverlaps(query=features, subject=gr), "List")
		weightedCount <- sum(extractList(score(gr), hits))

		# Return
		weightedCount
}

# High memory BED
quantifyFeatures.GRangesList <- function(ctss, features){
	# Quantify each GRange in GRanges list
	x <- pblapply(ctss, countScoredOverlaps, features=features)

	# Merge into matrix
	x <- do.call(cbind, x)
	rownames(x) <- names(features)

	# Return
	x
}

# High memory BED
quantifyFeatures.character <- function(ctss, features){
	# Quantify each GRange in GRanges list
	tmp_fun <- function(x, f){
		# Load file
		d <- fastReadBED(x)

		# Turn into GRanges
		d <- makeGRangesFromDataFrame(df=d,
																	keep.extra.columns=TRUE,
																	starts.in.df.are.0based=TRUE)

		# Overlap
		d <- countScoredOverlaps(gr=d, features=f)

		# Return
		d
	}

	x <- pblapply(ctss, tmp_fun, f=features)

	# Merge into matrix
	x <- do.call(cbind, x)
	rownames(x) <- names(features)

	# Return
	x
}
