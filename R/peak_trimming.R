#' @import S4Vectors
symmetricPercentiles <- function(r, prop){
	# Convert to normal vector
	r <- as.vector(r)

	# Calculate proportions
	r_sum <- sum(r)
	r_prop <- r / r_sum

	# Cumulate sum from either side
	f_cum <- cumsum(r_prop)
	r_cum <- rev(cumsum(rev(r_prop)))

	# Trim both sides
	half_prop <- prop / 2
	left <- Position(function(x) x >= half_prop, f_cum, right=FALSE, nomatch=1)
	right <- Position(function(x) x >= half_prop, r_cum, right = TRUE, nomatch=length(r))

	# Return range
	c(left, right)
}

#' @import S4Vectors
asymmetricPercentiles <- function(r, prop){
	# Convert to normal vector
	r <- as.vector(r)

	# Calculate proportions
	r_sum <- sum(r)
	r_prop <- r / r_sum

	# Sort be minimum cumsum from either side
	f_cum <- cumsum(r_prop)
	r_cum <- rev(cumsum(rev(r_prop)))
	min_cum <- pmin(f_cum, r_cum)
	o_cum <- order(min_cum)

	# Order original proportions and cut
	r_ord <- cumsum(r_prop[o_cum])
	left <- Position(function(x) x > prop, r_ord, right=FALSE, nomatch=1)
	o_out <- o_cum[left:length(r_ord)]

	# Return range
	range(o_out)
}

#' Trim width of TCs to expression percentiles
#'
#' Given a set of TCs and genome-wide CTSS coverage, reduce the width of TC until a certain amount of expression has been removed.
#'
#' @param TCs GRanges: TCs to be trimmed.
#' @param ctssCoverage GRanges: CTSS coverage.
#' @param percentile numeric: Fraction of expression to remove from TCs.
#' @param symmetric logical: Whether to trim the same amount from both edges of the TC (TRUE) or always trim from the least expressed end (FALSE).
#'
#' @return GRanges with trimmed TCs, including recalculated peaks and scores.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family TC trimming functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
trimToPercentiles <- function(TCs, ctssCoverage, percentile=0.1, symmetric=FALSE){
	# Checks
	stopifnot(class(TCs) == "GRanges",
						class(ctssCoverage) == "GRanges")
	stopifnot(identical(seqinfo(ctssCoverage), seqinfo(TCs)))

	# Split  by strand
	message("Splitting by strand...")
	coverage_stranded <- splitByStrand(ctssCoverage)
	coverage_plus <- coverage(coverage_stranded$`+`, weight="score")
	coverage_minus <- coverage(coverage_stranded$`-`, weight="score")
	tcs_stranded <- splitByStrand(TCs)
	rm(coverage_stranded)

	# Convert to IRangesList
	irl_plus <- methods::as(tcs_stranded$`+`,"RangesList")
	irl_minus <- methods::as(tcs_stranded$`-`,"RangesList")
	rm(tcs_stranded)

	# Obtain views
	views_plus <- Views(coverage_plus, irl_plus)
	views_minus <- Views(coverage_minus, irl_minus)

	# Extract Adjustments
	if(symmetric){
		message(paste0("Symmetric trimming to percentile: ", percentile * 100, "%"))
		adjust_plus <- viewApply(views_plus, symmetricPercentiles, prop=percentile, simplify=FALSE)
		adjust_minus <- viewApply(views_minus, symmetricPercentiles, prop=percentile, simplify=FALSE)
	}else if(!symmetric){
		message(paste0("Asymmetric trimming to percentile: ", percentile * 100, "%"))
		adjust_plus <- viewApply(views_plus, asymmetricPercentiles, prop=percentile, simplify=FALSE)
		adjust_minus <- viewApply(views_minus, asymmetricPercentiles, prop=percentile, simplify=FALSE)
	}else{
		stop("Additional percentile functions not yet implemented!")
	}
	rm(views_plus, views_minus)

	# Adjustments as IntegerLists
	message("Adjusting ranges...")
	left_plus <- methods::as(lapply(adjust_plus, function(x) lapply(x, function(x) x[1])), "IntegerList")
	right_plus <- methods::as(lapply(adjust_plus, function(x) lapply(x, function(x) x[2])), "IntegerList")
	left_minus <- methods::as(lapply(adjust_minus, function(x) lapply(x, function(x) x[1])), "IntegerList")
	right_minus <- methods::as(lapply(adjust_minus, function(x) lapply(x, function(x) x[2])), "IntegerList")

	# Narrow ranges
	irl_plus <- narrow(irl_plus, start=left_plus, end=right_plus)
	irl_minus <- narrow(irl_minus, start=left_minus, end=right_minus)
	rm(left_plus, right_plus, left_minus, right_minus)

	# Calculate new stats
	message("Calculating new stats...")
	trimmedTCs <- TCstats(coverage_plus=coverage_plus, coverage_minus=coverage_minus,
												tcs_plus=irl_plus, tcs_minus=irl_minus)
	rm(coverage_plus, coverage_minus, irl_plus, irl_minus)

	# Carry over seqinfo and sort
	message("Preparing output...")
	seqinfo(trimmedTCs) <- seqinfo(TCs)
	trimmedTCs <- sort(trimmedTCs)

	# Print some basic stats
	message("Tag clustering summary:")
	summarizeWidths(trimmedTCs)

	# Return
	trimmedTCs
}

#' Trim width of TCs by distance from TC peak
#'
#' Trim the width of TCs by distance from the TC peaks.
#'
#' @param TCs GRanges: TCs to be trimmed.
#' @param ctssCoverage GRanges: CTSS coverage.
#' @param upstream integer: Maximum upstream distance from TC peak.
#' @param downstream integer: Maximum downstream distance from TC peak.
#' @param peaks character: Name of column in TCs holding TC peaks as an IRanges.
#'
#' @return GRanges with trimmed TCs, including recalculated peaks and scores.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family TC trimming functions
#' @import S4Vectors IRanges GenomicRanges
#' @export
trimAroundPeaks <- function(TCs, ctssCoverage, upstream, downstream, peaks="thick"){
	# Checks
	stopifnot(class(TCs) == "GRanges",
						class(ctssCoverage) == "GRanges")
	stopifnot(identical(seqinfo(ctssCoverage), seqinfo(TCs)))

	# Extract peaks
	message("Trimming TCs around peaks...")
	peaks <- swapRanges(TCs)

	# Expand peaks
	expanded <- promoters(peaks, upstream=upstream, downstream=downstream)

	# Reset original range
	start(peaks) <- ifelse(start(expanded) >= start(TCs), start(expanded), start(TCs))
	end(peaks) <- ifelse(end(expanded) <= end(TCs), end(expanded), end(TCs))

	# Split  by strand
	message("Splitting by strand...")
	coverage_stranded <- splitByStrand(ctssCoverage)
	coverage_plus <- coverage(coverage_stranded$`+`, weight="score")
	coverage_minus <- coverage(coverage_stranded$`-`, weight="score")
	tcs_stranded <- splitByStrand(peaks)
	rm(coverage_stranded, peaks)

	# Convert to IRangesList
	irl_plus <- methods::as(tcs_stranded$`+`,"RangesList")
	irl_minus <- methods::as(tcs_stranded$`-`,"RangesList")
	rm(tcs_stranded)

	# Calculate new stats
	message("Calculating new stats...")
	trimmedTCs <- TCstats(coverage_plus=coverage_plus, coverage_minus=coverage_minus,
												tcs_plus=irl_plus, tcs_minus=irl_minus)
	rm(coverage_plus, coverage_minus, irl_plus, irl_minus)

	# Carry over seqinfo and sort
	message("Preparing output...")
	seqinfo(trimmedTCs) <- seqinfo(TCs)
	trimmedTCs <- sort(trimmedTCs)

	# Print some basic stats
	message("Tag clustering summary:")
	summarizeWidths(trimmedTCs)

	# Return
	trimmedTCs
}


### OLD CODE

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
# trimToQuantiles <- function(TCs, ctssCoverage, lower=0.05, upper=0.05){
# 	stop("Not yet implemented!")
# }
#
# trimToAssymetricQuantiles <- function(TCs, ctssCoverage, mass=0.9){
# 	stop("Not yet implemented!")
# }
