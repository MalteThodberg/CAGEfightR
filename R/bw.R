#' Check if BigWig-files are valid.
#'
#' Checks if a BigWigFile or BigWigFileList is composed of readable files with
#' the proper .bw extension.
#'
#' @param object BigWigFile or BigWigFileList
#'
#' @return TRUE, if any tests fails an error is raised.
#'
#' @family BigWig functions
#' @export
#' @examples
#' # Use the BigWig-files included with the package:
#' data('exampleDesign')
#' bw_plus <- system.file('extdata', exampleDesign$BigWigPlus,
#'                        package = 'CAGEfightR')
#'
#' # Create a named BigWigFileList-object with names
#' bw_plus <- BigWigFileList(bw_plus)
#' names(bw_plus) <- exampleDesign$Name
#'
#' # Check a single BigWigFile:
#' bwValid(bw_plus[[1]])
#'
#' # Check the entire BigWigFileList:
#' bwValid(bw_plus)
setGeneric("bwValid", function(object) {
    standardGeneric("bwValid")
})

#' @rdname bwValid
setMethod("bwValid", signature(object = "BigWigFile"), function(object) {
    # Checks, maybe wrap all of this in it's own function
    assert_that(file.exists(resource(object)),
                has_extension(resource(object), "bw"),
                is.readable(resource(object)))

    # Only returned if all tests pass
    TRUE
})

#' @rdname bwValid
setMethod("bwValid", signature(object = "BigWigFileList"), function(object) {
    # Checks, maybe wrap all of this in it's own function
    assert_that(all(vapply(object, bwValid, logical(1))), !is.null(names(object)))

    # Only returned if all tests pass
    TRUE
})

#' Find a common genome for a series of BigWig files.
#'
#' Finds a common genome for a series of BigWig-files, either using only levels
#' present in all files (intersect) or in any file (union).
#'
#' @param plusStrand BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param method character: Either 'intersect' or 'union'.
#'
#' @return Sorted Seqinfo-object.
#' @family BigWig functions
#' @export
#' @examples
#' if (.Platform$OS.type != "windows") {
#' # Use the BigWig-files included with the package:
#' data('exampleDesign')
#' bw_plus <- system.file('extdata', exampleDesign$BigWigPlus,
#'                        package = 'CAGEfightR')
#' bw_minus <- system.file('extdata', exampleDesign$BigWigMinus,
#'                         package = 'CAGEfightR')
#'
#' # Create two named BigWigFileList-objects:
#' bw_plus <- BigWigFileList(bw_plus)
#' bw_minus <- BigWigFileList(bw_minus)
#' names(bw_plus) <- exampleDesign$Name
#' names(bw_minus) <- exampleDesign$Name
#'
#' # Find the smallest common genome (intersect) across the BigWigList-objects:
#' bwCommonGenome(plusStrand=bw_plus, minusStrand=bw_minus, method='intersect')
#'
#' # Find the most inclusive genome (union) across the BigWigList-objects:
#' bwCommonGenome(plusStrand=bw_plus, minusStrand=bw_minus, method='union')
#' }
bwCommonGenome <- function(plusStrand, minusStrand, method = "intersect") {
    assert_that(methods::is(plusStrand, "BigWigFileList"),
                methods::is(minusStrand, "BigWigFileList"),
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
    if (method == "intersect") {
        o <- suppressWarnings(Reduce(f = intersect, seqInfo))
    } else if (method == "union") {
        o <- suppressWarnings(Reduce(f = merge, seqInfo))
    }

    # Sort
    o <- sortSeqlevels(o)

    # Return
    o
}

#' Check if BigWig-files are compatible with a given genome.
#'
#' Given a genome, checks whether a series of BigWig-files are compatible by
#' checking if all common seqlevels have equal seqlengths.
#'
#' @param plusStrand BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param minusStrand BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param genome Seqinfo: Genome information.
#'
#' @return TRUE, raises an error if the supplied genome is incompabtible.
#' @family BigWig functions
#' @export
#' @examples
#' if (.Platform$OS.type != "windows") {
#' # Use the BigWig-files included with the package:
#' data('exampleDesign')
#' bw_plus <- system.file('extdata', exampleDesign$BigWigPlus,
#'                        package = 'CAGEfightR')
#' bw_minus <- system.file('extdata', exampleDesign$BigWigMinus,
#'                         package = 'CAGEfightR')
#'
#' # Create two named BigWigFileList-objects:
#' bw_plus <- BigWigFileList(bw_plus)
#' bw_minus <- BigWigFileList(bw_minus)
#' names(bw_plus) <- exampleDesign$Name
#' names(bw_minus) <- exampleDesign$Name
#'
#' # Make a smaller genome:
#' si <- seqinfo(bw_plus[[1]])
#' si <- si['chr18']
#'
#' # Check if it is still compatible:
#' bwGenomeCompatibility(plusStrand=bw_plus, minusStrand=bw_minus, genome=si)
#' }
bwGenomeCompatibility <- function(plusStrand, minusStrand, genome) {
    assert_that(methods::is(plusStrand, "BigWigFileList"),
                methods::is(minusStrand, "BigWigFileList"),
                methods::is(genome, "Seqinfo"))

    # Get seqinfo
    seqInfoPlus <- lapply(plusStrand, seqinfo)
    seqInfoMinus <- lapply(minusStrand, seqinfo)

    # Unique seqinfos
    seqInfos <- unique(c(seqInfoPlus, seqInfoMinus))

    # Check if error is produces
    lapply(seqInfos, merge, y = genome)

    # Return
    TRUE
}
