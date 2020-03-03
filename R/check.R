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

#' Helper for checking files containing CTSSs
#'
#' Checks whether a file (or GRanges/GPos) contains data formatted in the same
#' manner as CAGE Transcription Start Sites (CTSSs): Each basepair of the genome
#' is associated with a single integer count.
#'
#' @param object BigWigFile, character, GRanges or GPos: Path to the file
#'   storing CTSSs, or an already improted GRanges/GPos.
#'
#' @return TRUE if CTSSs are correctly  formatted, otherwise a (hopefully)
#'   informative error is thrown.
#' @export
#'
#' @note In the case that a character is supplied pointing to a file, checkCTSSs
#'   will not check any extensions, but simply try to read it using
#'   rtracklayer::import. This means that checkCTSSs can technically analyze
#'   BED-files, although CAGEfightR can only import CTSSs from BigWig or
#'   bedGraph files.
#'
#' @examples
#' # Load example data
#' data('exampleDesign')
#' bw_plus <- system.file('extdata',
#'                        exampleDesign$BigWigPlus,
#'                        package = 'CAGEfightR')
#' bw_plus <- BigWigFileList(bw_plus)
#'
#' # Check raw file
#' checkCTSSs(bw_plus[[1]])
#'
#' # Import first, then check
#' gr  <- import(bw_plus[[1]])
#' checkCTSSs(gr)
setGeneric("checkCTSSs", function(object) {
    standardGeneric("checkCTSSs")
})

#' @rdname checkCTSSs
setMethod("checkCTSSs", signature(object = "ANY"), function(object) {
    stop("CTSSs can be checked as: a path to a BigWig/bedGraph file or a GRanges/GPos")
})

#' @rdname checkCTSSs
setMethod("checkCTSSs", signature(object = "GRanges"), function(object) {
    # Warning  for help
    if(length(object) == 0){
        warning("No CTSSs found: There are no ranges in the object!")
    }

    # Give informative error messages
    if(!isDisjoint(object)){
        stop("CTSSs are not disjoint: The ranges should not be overlapping!")
    }

    if(is.null(score(object))){
        stop("CTSSs do not have scores: Missing score column holding integer counts!")
    }

    if(any(score(object) %% 1 != 0)){
        stop("CTSSs do not have integer scores: The score column should contain the discrete number of 5'-ends of reads mapping to each bp!")
    }

    if(any(width(object) != 1)){
        stop("CTSSs are not bp-wise: Each range should be exactly 1 bp wide!\nNote: You can use `convertGRanges2GPos` function to correct for this.")
    }

    if(table(strand(object))["*"] > 0){
        warning("CTSSs are unstranded: CTSSs must be either on the + or - strand. ",
                "This warning can be ignored if a BigWig/bedGraph file was supplied, in which case each strand is stored in a separate file.")
    }

    # Return TRUE if check succeds
    TRUE
})

#' @rdname checkCTSSs
setMethod("checkCTSSs", signature(object = "character"), function(object) {
    # Make sure file is readable
    assert_that(is.string(object),
                file.exists(object),
                is.readable(object))

    # Attempt import
    message("Attempting import of supplied filename...")
    object <- import(object)

    # Check returned GRanges
    checkCTSSs(object)
})

#' @rdname checkCTSSs
setMethod("checkCTSSs", signature(object = "GPos"), function(object) {
    checkCTSSs(methods::as(object, "GRanges"))
})

#' @rdname checkCTSSs
setMethod("checkCTSSs", signature(object = "BigWigFile"), function(object) {
    # Make sure file is readable
    bwValid(object)

    # Attempt import
    object <- import(object)

    # Check returned GRanges
    checkCTSSs(object)
})
