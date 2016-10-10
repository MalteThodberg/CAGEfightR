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

# Generic for importing as GRanges list
importAsGRangesList <- function(ctss) UseMethod("importAsGRangesList")

importAsGRangesList.default <- function(ctss){
	stop("Supported classes: character, BigWigList or BAM")
}

importAsGRangesList.character <- function(ctss){
	# Import as data.table
	message("Reading CTSS-BED files")
	d <- pblapply(ctss, fastReadBED)

	# Concatenate
	message("Concatenating")
	d <- rbindlist(d, idcol="sample")

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
