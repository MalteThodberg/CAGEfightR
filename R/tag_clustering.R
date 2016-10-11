#' Simple Tag Clustering
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
#' @export
simpleTagClustering <- function(ctssCoverage, tpmCutoff=0, mergeDist=20){
	# Split coverage by strand
	message("Calculating coverage by strand")
	coverage_plus <- coverage(subset(ctssCoverage, strand == "+"), weight="score")
	coverage_minus <- coverage(subset(ctssCoverage, strand == "-"), weight="score")

	# Peak calling: Find peaks
	message("Finding peaks")
	peaks_plus <- slice(coverage_plus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)
	peaks_minus <- slice(coverage_minus, lower=tpmCutoff, upper=Inf, includeLower=FALSE, rangesOnly=TRUE)

	# Peak calling: Merge nearby peaks
	reduced_plus <- reduce(peaks_plus, min.gapwidth=mergeDist)
	reduced_minus <- reduce(peaks_minus, min.gapwidth=mergeDist)

	### Coverage and peaks
	message("Calculating global peak characteristics")
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
	message("Merging Tag Cluster info")
	TCs <- c(GRanges(reduced_plus, strand="+", thick=ranges_plus, score=sum_plus),
					 GRanges(reduced_minus, strand="-", thick=ranges_minus, score=sum_minus))

	# Add ids
	TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
	names(TCs) <- TC_ids

	# Return
	TCs
}
