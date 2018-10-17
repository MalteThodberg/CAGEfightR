#' Helper for checking pooled signal
#'
#' Checks whether a supplied pooled signal is valid: Single bp disjoint with
#' signal in the score column with supplied genome information.
#'
#' @param object GRanges or GPos: Pooled signal to be checked
#'
#' @return TRUE if object is correct format, otherwise an error is thrown
#' @family Checking functions
#' @export
#'
#' @examples
#' data(exampleCTSSs)
#' checkPooled(rowRanges(exampleCTSSs))
checkPooled <- function(object){
	# Checks
	x <- c(methods::is(object, "GenomicRanges"),
				 isDisjoint(object),
				 all(width(object) == 1),
				 !any(is.na(seqlengths(object))),
				 !is.null(score(object)),
				 is.numeric(score(object)))

	# Message if something is wrong:
	if(!all(x)){
	    stop("The supplied pooled signal is not valid!\n",
	         "A valid pooled signal must obey:\n",
             "- Stored as a GenomicRanges (eg. GRanges or GPos)\n",
             "- All ranges must have width=1 (single bp)\n",
             "- Ranges must be disjoint (non-overlapping)\n",
             "- Signal must be stored in the score metadata colum and be numeric\n",
             "- All seqlevels must have lengths (see seqlengths-function)")
	}

	# Return nothing
	TRUE
}

#' Helper for checking cluster with peaks
#'
#' Checks whether a supplied set of cluster have valid peaks: Whether the thick
#' column contains IRanges all contained within the main ranges.
#'
#' @param object GRanges or GPos: Clusters with peaks to be checked.
#'
#' @return TRUE if object is correct format, otherwise an error is thrown
#' @family Checking functions
#' @export
#'
#' @examples
#' data(exampleUnidirectional)
#' checkPeaked(rowRanges(exampleUnidirectional))
checkPeaked <- function(object){
	# Checks
    x <- c(methods::is(object, "GenomicRanges"),
           !is.null(object$thick),
           methods::is(object$thick, "IRanges"),
           all(poverlaps(mcols(object)$thick,
                         ranges(object),
                         type = "within")))

    # Message if something is wrong:
    if(!all(x)){
        stop("The the supplied regions does not have valid peaks!\n",
             "Valid regions with peaks must obey:\n",
             "- Stored as a GenomicRanges (eg. GRanges or GPos)\n",
             "- Have a metadata column named thick containing an IRanges-object\n",
             "- Main ranges must all contain the inner ranges in thick")
    }

    # Return nothing
    TRUE
}
