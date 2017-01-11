# #' Internal function: Quickly read simple BED-file
# #'
# #' Uses fread from data.table to quickly reading a simple bed file. Only scores are kept, any names are discarded
# #'
# #' @param fname character: Path to BED file.
# #'
# #' @return data.table
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @export
# fastReadBED <- function(fname){
# 	# TO DO:
# 	# Handle compressed files by calling gunzip.
#
# 	data.table::fread(input=fname,
# 										data.table=TRUE,
# 										sep="\t",
# 										col.names=c("seqnames", "start", "end", "score", "strand"),
# 										select=c(1,2,3,5,6),
# 										verbose=FALSE,
# 										showProgress=FALSE,
# 										na.strings=NULL)
# }
#
# #' Import CTSS data into a GRangesList
# #'
# #' Import a series of files containing CTSS data into a GRangesList. The file extensions are used to determine how the file is loaded, currently supportet is .bed, .bedGraph and .bw.
# #'
# #' @param ctss character: Vector of paths to CTSS-files.
# #' @param cores integer: Number of cores to use on import.
# #'
# #' @return GRangesList of the same length as ctss, containing CTSS scores for each file.
# #' @import S4Vectors IRanges GenomicRanges
# #' @export
# importAsGRangesList <- function(ctss, cores=NULL) UseMethod("importAsGRangesList")
#
# #' @describeIn importAsGRangesList Default method (always throws an error)
# #' @export
# importAsGRangesList.default <- function(ctss, cores=NULL){
# 	stop("Supported classes: character, BigWigList or BAM")
# }
#
# #' @describeIn importAsGRangesList Character method
# #' @export
# importAsGRangesList.character <- function(ctss, cores=NULL){
# 	# Import as data.table
# 	message("Reading CTSS-BED files")
# 	if(is.null(cores)){
# 		d <- pbapply::pblapply(ctss, fastReadBED)
# 	}else{
# 		message(paste0("Using cores: ", cores))
# 		d <- parallel::mclapply(ctss, fastReadBED,
# 											 mc.preschedule=FALSE, mc.cores=cores)
# 	}
#
# 	# Concatenate
# 	message("Concatenating")
# 	d <- data.table::rbindlist(d, idcol="sample")
#
# 	# Turn into GRanges
# 	message("Building GRangesList")
# 	# d <- makeGRangesFromDataFrame(df=d,
# 	# 																keep.extra.columns=TRUE,
# 	# 																starts.in.df.are.0based=TRUE)
# 	#
# 	# # Run length encode scores
# 	# score(d) <- Rle(score(d))
# 	# d$sample <- Rle(d$sample)
#
# 	d <- GRanges(seqnames=Rle(d$seqnames),
# 							 ranges=IRanges(d$start+1, width=1),
# 							 strand=Rle(d$strand),
# 							 score=Rle(d$score),
# 							 sample=Rle(d$sample))
#
# 	# Split into GRangesList
# 	d <- split(d, d$sample)
# 	names(d) <- basename(ctss)
#
# 	# Return
# 	d
# }
