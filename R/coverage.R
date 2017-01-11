# #' Calculate CTSS coverage
# #'
# #' Calculates TPM-coverage for a series of CTSS-files.
# #'
# #' @param ctss GRangesList or character: CTSS data stores as a GRangesList or a vector of paths to CTSS-files.
# #' @param ctssCutoff integer: Trim away positions with this or less CAGE-tags.
# #' @param cores integer: Number of cores to use for calculations.
# #'
# #' @return GRanges with bp-wise TPM-coverage.
# #'
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import S4Vectors IRanges GenomicRanges data.table
# #' @export
# coverageOfCTSS <- function(ctss, ctssCutoff=0, cores=NULL) UseMethod("coverageOfCTSS")
#
# #' @describeIn coverageOfCTSS default method (always throws an error)
# #' @export
# coverageOfCTSS.default <- function(ctss, ctssCutoff=0, cores=NULL){
# 	stop("ctss must be either a character or GRangesList")
# }
#
# #' @describeIn coverageOfCTSS High memory method
# #' @export
# coverageOfCTSS.GRangesList <- function(ctss, ctssCutoff=0, cores=NULL){
# 	# Don't use multiple cores
# 	if(!is.null(cores)){
# 		stop("Multiple cores not supported!")
# 	}
#
# 	message("Filtering CTSS & Calculating TPM")
# 	d <- as.data.table(ctss)
# 	d <- d[score > ctssCutoff,]
# 	d[, TPM := score / (sum(score)/1e6), by=group]
#
# 	message("Calculating global coverage")
# 	d <- d[,.(score=sum(TPM)), by = c("seqnames", "start", "end", "strand")]
#
# 	message("Building GRanges")
# 	d <- GRanges(seqnames=Rle(d$seqnames),
# 							 ranges=IRanges(d$start, width=1),
# 							 strand=Rle(d$strand),
# 							 score=d$score)
# 	# Return
# 	d
# }
#
# #' Internal function: Reduce-style coverage calculation
# #'
# #' Uses a reduce approach to calculate coverage. Each file is TPM-transformed prior to being reduced.
# #'
# #' @param ctss character: CTSS-BED files
# #' @param ctssCutoff integer: Trim away positions with this or less CAGE-tags.
# #' @param progressBar logical: Whether to track progress.
# #'
# #' @return data.table
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import data.table
# #' @export
# reduceCoverage <- function(ctss, ctssCutoff=0, progressBar=FALSE){
# 	### This should be made general!
#
# 	# Setup progress bar
# 	if(progressBar == TRUE){
# 		pb <- utils::txtProgressBar(min=0, max=length(ctss), style=3)
# 		pbCounter <- 0
# 		utils::setTxtProgressBar(pb, pbCounter)
# 	}
#
# 	# Read first file
# 	x <- fastReadBED(ctss[1])
#
# 	# Remove low CTSS
# 	x <- x[score > ctssCutoff,]
#
# 	# Calculate TPM
# 	x[, score := score / (sum(score)/1e6)]
#
# 	if(progressBar == TRUE){
# 		pbCounter <- pbCounter + 1
# 		utils::setTxtProgressBar(pb, pbCounter)
# 	}
#
# 	# Loop
# 	for(i in ctss[-1]){
# 		# Read first file
# 		y <- fastReadBED(i)
#
# 		# Remove low CTSS
# 		y <- y[score > ctssCutoff,]
#
# 		# Calculate TPM
# 		y[, score := score / (sum(score)/1e6)]
#
# 		# Concatenate
# 		x <- rbind(x, y)
#
# 		# Calculate Coverage
# 		x <- x[, .(score=sum(score)), by = c("seqnames", "start", "end", "strand")]
#
# 		if(progressBar == TRUE){
# 			pbCounter <- pbCounter + 1
# 			utils::setTxtProgressBar(pb, pbCounter)
# 		}
# 	}
#
# 	# Return
# 	x
# }
#
# #' Internal function: Merge data.tables with coverage
# #'
# #' To be used in a Reduced-loop, merge data.tables containing coverage information.
# #'
# #' @param dt1 data.table: Coverage in a data.table
# #' @param dt2 data.table: Coverage in a data.table
# #'
# #' @return Coverage in a data.table
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @import data.table
# #' @export
# mergeCoverage <- function(dt1, dt2){
# 	# Bind all data.frame
# 	dt <- rbindlist(list(dt1, dt2))
#
# 	# Merge without preprocess
# 	dt <- dt[,.(score=sum(score)), by=c("seqnames", "start", "end", "strand")]
#
# 	# Return
# 	dt
# }
#
# #' @describeIn coverageOfCTSS Low memory method
# #' @export
# coverageOfCTSS.character <- function(ctss, ctssCutoff=0, cores=NULL){
# 	if(is.null(cores)){
# 		#	Setup progress bar
# 		message("Processing files")
# 		d <- reduceCoverage(ctss=ctss, ctssCutoff=ctssCutoff, progressBar=TRUE)
#
# 	}else{
# 		# Split files into chunks, terminate if ill formed
# 		message(paste0("Processing files in chunks: ", cores))
# 		chunks <- split(ctss, factor(rep_len(seq_len(cores), length.out=length(ctss))))
# 		stopifnot(min(sapply(chunks, length)) > 1)
#
# 		# Process each chunk
# 		d <- parallel::mclapply(chunks, reduceCoverage, ctssCutoff=ctssCutoff, progressBar=FALSE,
# 												 mc.preschedule=FALSE, mc.cores=cores)
#
# 		# Merge results
# 		message("Merging chunks")
# 		d <- Reduce(mergeCoverage, d)
# 	}
#
# 	# To GRanges
# 	message("Building GRanges")
# 	d <- GRanges(seqnames=Rle(d$seqnames),
# 							 ranges=IRanges(d$start+1, width=1),
# 							 strand=Rle(d$strand),
# 							 score=d$score)
# 	# Return
# 	d
# }
#
# # #' @describeIn coverageOfCTSS Low memory method
# # #' @export
# # coverageOfCTSS.character <- function(ctss, ctssCutoff=0, cores=NULL){
# # 	if(is.null(cores)){
# # 		#	Setup progress bar
# # 		message("Reading files")
# # 		pb <- utils::txtProgressBar(min=0, max=length(ctss), style=3)
# # 		pbCounter <- 0
# # 		utils::setTxtProgressBar(pb, pbCounter)
# #
# # 		# Read first file
# # 		x <- fastReadBED(ctss[1])
# #
# # 		# Remove low CTSS
# # 		x <- x[score > ctssCutoff,]
# #
# # 		# Calculate TPM
# # 		x[, score := score / (sum(score)/1e6)]
# #
# # 		pbCounter <- pbCounter + 1
# # 		utils::setTxtProgressBar(pb, pbCounter)
# #
# # 		for(i in ctss[-1]){
# # 			# Read first file
# # 			y <- fastReadBED(i)
# #
# # 			# Remove low CTSS
# # 			y <- y[score > ctssCutoff,]
# #
# # 			# Calculate TPM
# # 			y[, score := score / (sum(score)/1e6)]
# #
# # 			# Concatenate
# # 			x <- rbind(x, y)
# #
# # 			# Calculate Coverage
# # 			x <- x[, .(score=sum(score)), by = c("seqnames", "start", "end", "strand")]
# #
# # 			pbCounter <- pbCounter + 1
# # 			utils::setTxtProgressBar(pb, pbCounter)
# # 		}
# # 	}else{
# # 		# Split files into chunks
# # 		message(paste0("Splitting files into chunks: ", cores))
# # 		chunks <- splitFiles(ctss, cores)
# #
# # 		# Process each chunk
# # 		chunked_res <- mclapply(chunks, reduceCoverage, ctssCutoff=ctssCutoff,
# # 												 mc.preschedule=FALSE, mc.cores=cores)
# #
# # 		# Merge results
# # 		x <- Reduce(mergeCoverage, innerRes)
# # }
# #
# # 	# To GRanges
# # 	message("Building GRanges")
# # 	# makeGRangesFromDataFrame(df=x,
# # 	# 												 keep.extra.columns=TRUE,
# # 	# 												 starts.in.df.are.0based=TRUE)
# # 	GRanges(seqnames=Rle(x$seqnames),
# # 							 ranges=IRanges(x$start+1, width=1),
# # 							 strand=Rle(x$strand),
# # 							 score=x$score)
# #
# # }
#
