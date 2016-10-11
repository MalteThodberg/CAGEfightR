#' Calculate CTSS coverage
#'
#' Calculates TPM-coverage for a series of CTSS-files.
#'
#' @param ctss GRangesList or character: CTSS data stores as a GRangesList or a vector of paths to CTSS-files.
#' @param ctssCutoff integer: Trim away positions with this or less CAGE-tags.
#'
#' @return GRanges with bp-wise TPM-coverage.
#'
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import data.table
#' @export
coverageOfCTSS <- function(ctss, ctssCutoff=0) UseMethod("coverageOfCTSS")

#' @describeIn coverageOfCTSS default method (always throws an error)
#' @export
coverageOfCTSS.default <- function(ctss, ctssCutoff=0){
	stop("ctss must be either a character or GRangesList")
}

#' @describeIn coverageOfCTSS High memory method
#' @export
coverageOfCTSS.GRangesList <- function(ctss, ctssCutoff=0){
	message("Filtering CTSS & Calculating TPM")
	d <- as.data.table(ctss)
	d <- d[score > ctssCutoff,]
	d[, TPM := score / (sum(score)/1e6), by=group]

	message("Calculating global coverage")
	d <- d[,.(score=sum(TPM)), by = c("seqnames", "start", "end", "strand")]

	message("Building GRanges")
	d <- GenomicRanges::makeGRangesFromDataFrame(df=d,
																 keep.extra.columns=TRUE,
																 starts.in.df.are.0based=FALSE)
	# Return
	d
}

#' @describeIn coverageOfCTSS Low memory method
#' @export
coverageOfCTSS.character <- function(ctss, ctssCutoff=0){
	# Setup progress bar
	message("Reading files")
	pb <- utils::txtProgressBar(min=0, max=length(ctss), style=3)
	pbCounter <- 0
	utils::setTxtProgressBar(pb, pbCounter)

	# Read first file
	x <- fastReadBED(ctss[1])

	# Remove low CTSS
	x <- x[score > ctssCutoff,]

	# Calculate TPM
	x[, score := score / (sum(score)/1e6)]

	pbCounter <- pbCounter + 1
	utils::setTxtProgressBar(pb, pbCounter)

	for(i in ctss[-1]){
		# Read first file
		y <- fastReadBED(i)

		# Remove low CTSS
		y <- y[score > ctssCutoff,]

		# Calculate TPM
		y[, score := score / (sum(score)/1e6)]

		# Concatenate
		x <- rbind(x, y)

		# Calculate Coverage
		x <- x[, .(score=sum(score)), by = c("seqnames", "start", "end", "strand")]

		pbCounter <- pbCounter + 1
		utils::setTxtProgressBar(pb, pbCounter)
}

	# To GRanges
	message("Building GRanges")
	makeGRangesFromDataFrame(df=x,
													 keep.extra.columns=TRUE,
													 starts.in.df.are.0based=TRUE )
}
