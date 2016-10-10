#' Swap ranges
#'
#' Swap out the range of a GRanges-object with another IRanges-object stored inside the same object. I.e., swapping TC widths with TC peaks.
#'
#' @param gr GRanges: GRanges
#' @param column character: Name of column in GRanges containing an IRanges.
#'
#' @return GRanges with swapped ranges. Old ranges are in the column "swap"
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
swapRanges <- function(gr, column="thick"){
	# TO DO
	# Check all columns are actually ranges.
	# Check that the swap column does not exist.
	# Allow other name than swap in function.

	# Copy ranges to a separate column
	gr$swap <- ranges(gr)

	# Switch in the new column
	ranges(gr) <- mcols(gr)[,column]

	# Use original names
	names(gr) <- names(gr$swap)

	# Return
	gr
}

### Not yet implemented peak trimming

# # Calculate peak stats based on coverage
# calculatePeakStats <- function(gr, ctssCoverage){
# 	# Split coverage by strand
# 	message("Calculating coverage by strand")
# 	coverage_plus <- coverage(subset(ctssCoverage, strand == "+"), weight="score")
# 	coverage_minus <- coverage(subset(ctssCoverage, strand == "-"), weight="score")
#
# 	# Peak calling: Find peaks
# 	message("Extracting peaks")
# 	gr_plus <- as(subset(gr, strand == "+"), "RangesList")
# 	gr_minus <- as(subset(gr, strand == "-"), "RangesList")
#
# 	### Coverage and peaks
# 	message("Calculating global peak characteristics")
# 	# Peaks stats: summed TPM
# 	views_plus <- Views(coverage_plus, gr_plus)
# 	views_minus <- Views(coverage_minus, gr_minus)
# 	sum_plus <- unlist(viewSums(views_plus))
# 	sum_minus <- unlist(viewSums(views_minus))
#
# 	# Peaks stats: Dominant TSS
# 	peaks_plus <- unlist(viewWhichMaxs(views_plus))
# 	peaks_minus <- unlist(viewWhichMaxs(views_minus))
# 	ranges_plus <- IRanges(start=peaks_plus, width=1)
# 	ranges_minus <- IRanges(start=peaks_minus, width=1)
#
# 	### Merge into final TC set
# 	message("Merging Tag Cluster info")
# 	TCs <- c(GRanges(gr_plus, strand="+", thick=ranges_plus, score=sum_plus),
# 					 GRanges(gr_minus, strand="-", thick=ranges_minus, score=sum_minus))
#
# 	# Add ids
# 	TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
# 	names(TCs) <- TC_ids
#
# 	# Return
# 	TCs
# }
#
# # Function for scaling peaks around
# trimAroundPeaks <- function(gr, ctssCoverage, size, peaks="thick"){
# 	# Extract peaks
# 	message("Trimming features around peaks")
# 	peaks <- swapRanges(gr)
#
# 	# Expand peaks
# 	expanded <- promoters(peaks, size, size)
#
# 	# Reset original range
# 	start(peaks) <- ifelse(start(expanded) >= start(gr), start(expanded), start(gr))
# 	end(peaks) <- ifelse(end(expanded) <= end(gr), end(expanded), end(gr))
#
# 	# Recalculate stats
# 	peaks <- calculatePeakStats(gr=peaks, ctssCoverage=ctssCoverage)
#
# 	# Return
# 	peaks
# }
#
# trimToQuantiles <- function(TCs, ctssCoverage, lower=0.05, upper=0.05){
# 	stop("Not yet implemented!")
# }
#
# trimToAssymetricQuantiles <- function(TCs, ctssCoverage, mass=0.9){
# 	stop("Not yet implemented!")
# }
