#' Check if BigWig-files are valid.
#'
#' Checks if a BigWigFile or BigWigFileList is composed of readable files with the proper .bw extension.
#'
#' @param object BigWigFile or BigWigFileList
#'
#' @return TRUE, if any tests fails an error is raised.
#'
#' @family BigWig functions
#' @import rtracklayer assertthat
#' @export
setGeneric("bwValid", function(object) {
	standardGeneric("bwValid")
})

#' @import rtracklayer assertthat
#' @rdname bwValid
setMethod("bwValid", signature(object="BigWigFile"), function(object){
	# Checks, maybe wrap all of this in it's own function
	assert_that(file.exists(object@resource),
							has_extension(object@resource, "bw"),
							is.readable(object@resource))

	# Only returned if all tests pass
	TRUE
})

#' @import rtracklayer assertthat
#' @rdname bwValid
setMethod("bwValid", signature(object="BigWigFileList"), function(object){
	# Checks, maybe wrap all of this in it's own function
	assert_that(all(vapply(object, bwValid, logical(1))),
							!is.null(names(object)))

	# Only returned if all tests pass
	TRUE
})

#' Find a common genome for a series of BigWig files.
#'
#' Finds a common genome for a series of BigWig-files, either using only levels present in all files (intersect) or in any file (union).
#'
#' @param plusStrand BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param method character: Either "intersect" or "union".
#'
#' @return Sorted Seqinfo-object.
#' @family BigWig functions
#' @import rtracklayer GenomeInfoDb assertthat
#' @export
bwCommonGenome <- function(plusStrand, minusStrand, method="intersect"){
	assert_that(class(plusStrand) == "BigWigFileList",
							class(minusStrand) == "BigWigFileList",
							is.string(method),
							method %in% c("intersect", "union"))

	# Get seqinfo
	seqInfoPlus <- lapply(plusStrand, seqinfo)
	seqInfoMinus <- lapply(minusStrand, seqinfo)

	# Unique seqinfos
	seqInfo <- unique(c(seqInfoPlus, seqInfoMinus))

	# Sort every element
	seqInfo <- seqInfo[order(vapply(seqInfo, length, numeric(1)))]

	# Merge using either function
	if(method == "intersect"){
		o <- suppressWarnings(Reduce(f=intersect, seqInfo))
	}else if(method == "union"){
		o <- suppressWarnings(Reduce(f=merge, seqInfo))
	}

	# Sort
	o <- sortSeqlevels(o)

	# Return
	o
}

#' Check if BigWig-files are compatible with a given genome.
#'
#' Given a genome, checks whether a series of BigWig-files are compatible by checking if all common seqlevels have equal seqlengths.
#'
#' @param plusStrand BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param genome Seqinfo: Genome information.
#'
#' @return TRUE, raises an error if the supplied genome is incompabtible.
#' @import rtracklayer GenomeInfoDb assertthat
#' @family BigWig functions
#' @export
bwGenomeCompatibility <- function(plusStrand, minusStrand, genome){
	assert_that(class(plusStrand) == "BigWigFileList",
							class(minusStrand) == "BigWigFileList",
							class(genome) == "Seqinfo")

	# Get seqinfo
	seqInfoPlus <- lapply(plusStrand, seqinfo)
	seqInfoMinus <- lapply(minusStrand, seqinfo)

	# Unique seqinfos
	seqInfos <- unique(c(seqInfoPlus, seqInfoMinus))

	# Check if error is produces
	lapply(seqInfos, merge, y=genome)

	# Return
	TRUE
}
