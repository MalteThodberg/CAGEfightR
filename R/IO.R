#' Quickly read simple BED-file
#'
#' Uses fread from data.table to quickly reading a simple bed file. Only scores are kept, any names are discarded
#'
#' @param fname character: Path to BED file.
#'
#' @return data.table
#' @examples
#' # ADD_EXAMPLES_HERE
fastReadBED <- function(fname){
	# TO DO:
	# Handle compressed files by calling gunzip.

	data.table::fread(input=fname,
										data.table=TRUE,
										sep="\t",
										col.names=c("seqnames", "start", "end", "score", "strand"),
										select=c(1,2,3,5,6),
										verbose=FALSE,
										showProgress=FALSE,
										na.strings=NULL)
}

#' Import CTSS data into a GRangesList
#'
#' Import a series of files containing CTSS data into a GRangesList. The file extensions are used to determine how the file is loaded, currently supportet is .bed, .bedGraph and .bw.
#'
#' @param ctss character: Vector of paths to CTSS-files.
#'
#' @return GRangesList of the same length as ctss, containing CTSS scores for each file.
#' @import S4Vectors IRanges GenomicRanges
#' @export
importAsGRangesList <- function(ctss) UseMethod("importAsGRangesList")

#' @describeIn importAsGRangesList Default method (always throws an error)
#' @export
importAsGRangesList.default <- function(ctss){
	stop("Supported classes: character, BigWigList or BAM")
}

#' @describeIn importAsGRangesList Character method
#' @export
importAsGRangesList.character <- function(ctss){
	# Import as data.table
	message("Reading CTSS-BED files")
	d <- pbapply::pblapply(ctss, fastReadBED)

	# Concatenate
	message("Concatenating")
	d <- data.table::rbindlist(d, idcol="sample")

	# Turn into GRanges
	message("Building GRangesList")
	d <- makeGRangesFromDataFrame(df=d,
																	keep.extra.columns=TRUE,
																	starts.in.df.are.0based=TRUE)

	# Run length encode scores
	score(d) <- Rle(score(d))
	d$sample <- Rle(d$sample)

	# Split into GRangesList
	d <- split(d, d$sample)
	names(d) <- basename(ctss)

	# Return
	d
}
