#' @import S4Vectors IRanges GenomicRanges
TCstats <- function(coverage_plus, coverage_minus, tcs_plus, tcs_minus){
	# Check classes
	stopifnot(class(coverage_plus)=="SimpleRleList",
						class(coverage_minus)=="SimpleRleList",
						class(tcs_plus)=="CompressedIRangesList",
						class(tcs_minus)=="CompressedIRangesList")

	# Check seqlevels
	stopifnot(length(unique(list(names(coverage_plus),
															 names(tcs_plus),
															 names(coverage_minus),
															 names(tcs_minus)))) == 1)

	# Obtain views
	views_plus <- Views(coverage_plus, tcs_plus)
	views_minus <- Views(coverage_minus, tcs_minus)

	# Calculate Sums
	sum_plus <- unlist(viewSums(views_plus))
	sum_minus <- unlist(viewSums(views_minus))

	# Find peaks
	ranges_plus <- viewRangeMaxs(views_plus)
	ranges_minus <- viewRangeMaxs(views_minus)
	ranges_plus <- resize(unlist(ranges_plus), width=1, fix="center")
	ranges_minus <- resize(unlist(ranges_minus), width=1, fix="center")

	# Merge into GRanges
	TCs <- c(GRanges(tcs_plus, strand="+", thick=ranges_plus, score=sum_plus),
					 GRanges(tcs_minus, strand="-", thick=ranges_minus, score=sum_minus))

	# Names as IDs for both ranges and peaks
	TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
	names(TCs) <- TC_ids
	names(TCs$thick) <- TC_ids

	# Return
	TCs
}

#' @import S4Vectors IRanges GenomicRanges
summarizeWidths <- function(gr){
	# Checks
	stopifnot(class(gr) == "GRanges")

	# Cut up widths
	x <- cut(width(gr),
					 breaks=c(1, 10, 100, 1000, Inf),
					 labels=c(">=1", ">=10", ">=100", ">=1000"),
					 include.lowest=TRUE)

	# Get freqs and props
	y <- table(Width=x)
	z <- prop.table(y)

	# Format to data.frame
	w <- merge(as.data.frame(y, responseName="Count"),
						 as.data.frame(z, responseName="Percent"))

	# Add Total row
	w <- rbind(data.frame(Width="Total",
												Count=sum(w$Count),
												Percent=sum(w$Percent)),
						 w)

	# Reformat to percent
	w$Percent <- paste0(format(w$Percent * 100, digits=1), " %")

	# To string and message
	s <- paste(utils::capture.output(print(w, row.names=FALSE)), collapse="\n")
	message(s)
}

#' Tag Clustering of genome-wide coverage
#'
#' Finds Tag Cluster (TCs) with a TPM above a certain threshold. Addtionally calculates the sum and peak position of the TCs.
#'
#' @param ctssCoverage GRanges: CTSS coverage.
#' @param tpmCutoff numeric: Minium TPM coverage to be considered as TC.
#' @param mergeDist integer: Merge TCs within this distance.
#'
#' @return GRanges with TPM sum as the score column, and TC peak as the thick column. Genome info is carried over
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Tag-clustering functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
tagClustering <- function(ctssCoverage, tpmCutoff=0, mergeDist=20){
	# Split coverage by strand
	message("Splitting by strand...")
	coverage_stranded <- splitByStrand(ctssCoverage)
	coverage_plus <- coverage(coverage_stranded$`+`, weight="score")
	coverage_minus <- coverage(coverage_stranded$`-`, weight="score")
	rm(coverage_stranded)

	# Peak calling: Find peaks
	message("Finding tag clusters...")
	peaks_plus <- slice(coverage_plus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
	peaks_minus <- slice(coverage_minus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)

	# Peak calling: Merge nearby peaks
	reduced_plus <- reduce(peaks_plus, min.gapwidth=mergeDist)
	reduced_minus <- reduce(peaks_minus, min.gapwidth=mergeDist)

	### Coverage and peaks
	message("Calculating tag cluster statistics...")
	TCs <- TCstats(coverage_plus=coverage_plus, coverage_minus=coverage_minus,
												tcs_plus=reduced_plus, tcs_minus=reduced_minus)

	# Carry over seqinfo and sort
	message("Preparing output...")
	seqinfo(TCs) <- seqinfo(ctssCoverage)
	TCs <- sort(TCs)

	# Print some basic stats
	message("Tag clustering summary:")
	summarizeWidths(TCs)

	# Return
	TCs
}

#' Simple Tag Clustering - DEPRECATED
#'
#' Finds Tag Cluster (TCs) with a TPM above a certain threshold. Addtionally calculates the sum and peak position of the TCs.
#'
#' @param ctssCoverage GRanges: CTSS coverage.
#' @param tpmCutoff numeric: Minium TPM coverage to be considered as TC.
#' @param mergeDist integer: Merge TCs within this distance.
#'
#' @return GRanges with TPM sum as the score column, and TC peak as the thick column.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Tag-clustering functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
simpleTagClustering <- function(ctssCoverage, tpmCutoff=0, mergeDist=20){
	# Mark as deprecated
	.Deprecated("tagClustering")

	# Split coverage by strand
	message("Calculating coverage by strand...")
	coverage_plus <- coverage(subset(ctssCoverage, strand == "+"), weight="score")
	coverage_minus <- coverage(subset(ctssCoverage, strand == "-"), weight="score")

	# Peak calling: Find peaks
	message("Finding peaks...")
	peaks_plus <- slice(coverage_plus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
	peaks_minus <- slice(coverage_minus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)

	# Peak calling: Merge nearby peaks
	reduced_plus <- reduce(peaks_plus, min.gapwidth=mergeDist)
	reduced_minus <- reduce(peaks_minus, min.gapwidth=mergeDist)

	### Coverage and peaks
	message("Calculating global peak characteristics...")
	# Peaks stats: summed TPM
	views_plus <- Views(coverage_plus, reduced_plus)
	views_minus <- Views(coverage_minus, reduced_minus)
	sum_plus <- unlist(viewSums(views_plus))
	sum_minus <- unlist(viewSums(views_minus))

	# Peaks stats: Dominant TSS
	peaks_plus <- unlist(viewWhichMaxs(views_plus))
	peaks_minus <- unlist(viewWhichMaxs(views_minus))
	ranges_plus <- IRanges(start=peaks_plus, width=1)
	ranges_minus <- IRanges(start=peaks_minus, width=1)

	### Merge into final TC set
	message("Merging Tag Cluster info...")
	TCs <- c(GRanges(reduced_plus, strand="+", thick=ranges_plus, score=sum_plus),
					 GRanges(reduced_minus, strand="-", thick=ranges_minus, score=sum_minus))

	# Copy seqinfo
	seqinfo(TCs) <- seqinfo(ctssCoverage)

	# Add ids
	TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
	names(TCs) <- TC_ids

	# Return
	TCs
}

#' @import S4Vectors IRanges GenomicRanges
countClusters <- function(thresholds, cv, mergeDist=20){
	# Slice
	o <- lapply(thresholds, slice, x=cv, includeLower=FALSE, upper=Inf, rangesOnly=TRUE)

	# Reduce
	o <- lapply(o, reduce, min.gapwidth=mergeDist)

	# Count
	o <- lapply(o, elementNROWS)

	# Sum
	o <- vapply(o, sum, integer(1))

	# Return
	o
}

#' Determine optimal expression threshold for Tag Clustering.
#'
#' This function counts the number of Tag Clusters (TCs) for an exponentially increasing series of expression cutoffs.
#'
#' @param ctssCoverage GRanges: CTSS coverage in score column.
#' @param stepSize integer: Size of steps for increasing threshold
#' @param mergeDist integer: Merge TCs within this distance.
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#'
#' @return data.frame with two columns: threshold and nTCs (number of Tag Clusters)
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family Tag-clustering functions
#' @import S4Vectors IRanges GenomicRanges BiocParallel
#' @export
tuneTagClustering <- function(ctssCoverage, stepSize=10, mergeDist=20, biocParallel=bpparam()){
	# Split coverage by strand
	message("Splitting coverage by strand...")
	by_strand <- splitByStrand(ctssCoverage)
	coverage_plus <- coverage(by_strand$`+`, weight="score")
	coverage_minus <- coverage(by_strand$`-`, weight="score")
	rm(by_strand)

	# Highest and lowest coverage
	min_val <- min(score(ctssCoverage))
	max_val <- max(score(ctssCoverage))
	message("Min / max CTSS coverage: ", min_val, " / ", max_val)

	# Exponential series from min to mx
	exp_series <- min_val * 2^seq(0, 100, stepSize)
	tpms <- c(0, exp_series[exp_series < max_val])
	rm(min_val, max_val, exp_series)
	message("Number of thresholds scored: ", length(tpms))

	# Count
	message("Plus strand...")
	n_plus <- bpvec(tpms, countClusters, cv=coverage_plus, mergeDist=mergeDist, BPPARAM=biocParallel)
	message("Minus strand...")
	n_minus <- bpvec(tpms, countClusters, cv=coverage_minus, mergeDist=mergeDist, BPPARAM=biocParallel)

	# Assemble results
	message("Preparing output...")
	o <- data.frame(threshold=tpms, nTCs=n_plus+n_minus)

	# Return
	o
}

# OLD CODE #################################

#' #' Tag Clustering of genome-wide coverage
#' #'
#' #' Finds Tag Cluster (TCs) with a TPM above a certain threshold. Addtionally calculates the sum and peak position of the TCs.
#' #'
#' #' @param ctssCoverage GRanges: CTSS coverage.
#' #' @param tpmCutoff numeric: Minium TPM coverage to be considered as TC.
#' #' @param mergeDist integer: Merge TCs within this distance.
#' #'
#' #' @return GRanges with TPM sum as the score column, and TC peak as the thick column. Genome info is carried over
#' #' @examples
#' #' # ADD_EXAMPLES_HERE
#' #' @family Tag-clustering functions
#' #' @import S4Vectors IRanges GenomicRanges
#' #' @export
#' tagClustering <- function(ctssCoverage, tpmCutoff=0, mergeDist=20){
#' 	# Split coverage by strand
#' 	message("Splitting coverage by strand...")
#' 	coverage_stranded <- splitByStrand(ctssCoverage)
#' 	coverage_plus <- coverage(coverage_stranded$`+`, weight="score")
#' 	coverage_minus <- coverage(coverage_stranded$`-`, weight="score")
#' 	rm(coverage_stranded)
#'
#' 	# Peak calling: Find peaks
#' 	message("Finding peaks...")
#' 	peaks_plus <- slice(coverage_plus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
#' 	peaks_minus <- slice(coverage_minus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
#'
#' 	# Peak calling: Merge nearby peaks
#' 	reduced_plus <- reduce(peaks_plus, min.gapwidth=mergeDist)
#' 	reduced_minus <- reduce(peaks_minus, min.gapwidth=mergeDist)
#'
#' 	### Coverage and peaks
#' 	message("Calculating global peak characteristics...")
#' 	# Peaks stats: summed TPM
#' 	views_plus <- Views(coverage_plus, reduced_plus)
#' 	views_minus <- Views(coverage_minus, reduced_minus)
#' 	sum_plus <- unlist(viewSums(views_plus))
#' 	sum_minus <- unlist(viewSums(views_minus))
#'
#' 	# Peaks stats: Dominant TSS
#' 	ranges_plus <- viewRangeMaxs(views_plus)
#' 	ranges_minus <- viewRangeMaxs(views_minus)
#' 	ranges_plus <- resize(unlist(ranges_plus), width=1, fix="center")
#' 	ranges_minus <- resize(unlist(ranges_minus), width=1, fix="center")
#'
#' 	### Merge into final TC set
#' 	message("Merging Tag Cluster info...")
#' 	TCs <- c(GRanges(reduced_plus, strand="+", thick=ranges_plus, score=sum_plus),
#' 					 GRanges(reduced_minus, strand="-", thick=ranges_minus, score=sum_minus))
#'
#' 	# Copy seqinfo
#' 	seqinfo(TCs) <- seqinfo(ctssCoverage)
#'
#' 	# Add ids
#' 	TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
#' 	names(TCs) <- TC_ids
#'
#' 	# Return
#' 	TCs
#' }
