#' Swap ranges in a GRanges.
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

#' Extend GRanges upstrean and/or downstream
#'
#' Extends the ranges of GRanges upstream and/or downstream, while retaining associated information in mcols.
#'
#' @param gr GRanges: Ranges to be extended, any seqinfo will be used to trim out-of-bounds ranges.
#' @param upstream integer: Number of basepairs to extend upstream.
#' @param downstream integer: Number of basepairs to extend downstream.
#'
#' @return GRanges with ranges extended.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @author Hervé Pagès
#' @import S4Vectors IRanges GenomicRanges
#' @export
extendRanges <- function(gr, upstream=0, downstream=0){
	# Check strand
	if (any(strand(gr) == "*"))
		warning("'*' ranges were treated as '+'")
	on_plus <- strand(gr) == "+" | strand(gr) == "*"

	# Extend ranges
	new_start <- start(gr) - ifelse(on_plus, upstream, downstream)
	new_end <- end(gr) + ifelse(on_plus, downstream, upstream)

	# New ranges
	ranges(gr) <- IRanges(new_start, new_end, names=names(gr))

	# Copy names
	gr
}

#' Finding overlapping genomic ranges resolved for multiple overlaps.
#'
#' Similar to findOverlaps, but provides to additional methods for resolving the case of a query overlapping multiple subjects using subject lenghts.
#'
#' @param query GRanges or GRangesList: Passed directly to findOverlaps.
#' @param subject GRanges or GRangesList: Passed directly to findOverlaps.
#' @param select character: How to resolve multiple overlaps: findOverlap's "first", "last", "arbitrary" or additional "shortest" or "longest".
#' @param ... additional arguments passed to findOverlaps.
#'
#' @return Similar to findOverlaps with select="first", "last" or "arbitrary".
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors IRanges GenomicRanges
#' @export
resolvedOverlaps <- function(query, subject, select, ...){
	stopifnot(select %in% c("first", "last", "arbitrary", "longest", "shortest"))

	if(select %in% c("first", "last", "arbitrary")){
		o <- findOverlaps(query=query, subject=subject, select=select, ...)
	}else if(select %in% c("longest", "shortest")){
		# Hits and widths as a list
		hits <- findOverlaps(query=query, subject=subject, select="all", ...)
		hits_ids <- methods::as(hits, "List")
		hits_widths <- extractList(width(subject), hits_ids)

		# Either longest or shortests
		if(select != "all"){
			# Find genes
			o <- switch(select,
									longest=hits_ids[hits_widths == max(hits_widths)],
									shortest=hits_ids[hits_widths == min(hits_widths)])

			# Clean in case some genes have the same lenght
			el <- elementNROWS(o)
			if(max(el) > 1){
				warning("Some genes have the same length! Resolving by using first-rule")
				o[el > 1] <- lapply(o[el > 1], function(x) x[1])
			}

			# Coerce to vector
			o <- as.integer(o)
		}
	}else{
		stop("select-argument not supported!")
	}

	# Checks
	stopifnot(length(o) == length(query),
						is.integer(o))

	# Return
	o
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
