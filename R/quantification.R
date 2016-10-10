quantifyFeatures <- function(ctss, features, cores) UseMethod("quantifyFeatures")

# Default method
quantifyFeatures.default <- function(ctss, features, cores){
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
quantifyFeatures.GRangesList <- function(ctss, features, cores=NULL){
	# Get some basic info for printing
	nFeatures <- length(features)
	nSamples <- length(ctss)
	m <- sprintf("Quantifying expression of %d features in %d samples",
							 nFeatures, nSamples)
	message(m)

	# Quantify each GRange in GRanges list
	if(is.null(cores)){
		x <- pblapply(ctss, countScoredOverlaps, features=features)
	}else{
		message(paste0("Using cores: ", cores))
		x <- mclapply(ctss, countScoredOverlaps, features=features,
									mc.preschedule=TRUE, mc.cores=cores)
	}

	# Merge into matrix
	message("Building Expression Matrix")
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
