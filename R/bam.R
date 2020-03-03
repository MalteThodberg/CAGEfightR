bam2bw <- function(i, o_plus, o_minus, minLength,...){
    # Attempt to read into R
    message("Reading BAM-file into memory: ", basename(i))
    o <- import(i, ...)

    # Cut reads and resize
    message("Finding 5'-ends...")
    o <- methods::as(o, "GRanges")
    o <- subset(o, width >= minLength)
    o <- resize(o, width = 1, fix = "start")

    # Count
    message("Locating CTSSs and counting tags...")
    o <- GenomicRanges::reduce(o, min.gapwidth = 0, with.revmap = TRUE)
    o <- methods::as(o, "StitchedGPos")
    score(o) <- elementNROWS(o$revmap)
    o$revmap <- NULL

    # Make sure output is ok
    message("Sorting and checking output...")
    o <- sort(o)
    n_tags <- sum(score(o))
    n_CTSSs <- length(o)

    # Split by strand
    message("Writing and writing BigWig-files...")
    by_strand <- utilsDeStrand(methods::as(o, "GRanges"))
    rm(o)

    # Write to file
    export.bw(object=by_strand$`+`, con=o_plus, format="bw")
    export.bw(object=by_strand$`-`, con=o_minus, format="bw")

    # Return
    c(CTSSs=n_CTSSs, Tags=n_tags)
}

#' Extract CTSSs from BAM-files (EXPERIMENTAL)
#'
#' Function for converting mapped reads in BAM-files to CAGE Transcription Start
#' Sites (CTSSs) in BigWig-files. Currently, this function will simply load a
#' (single-end) BAM-file (respecting a supplied ScanBamParam), optionally remove
#' short tags, and count the number of 5'-ends at each bp. Note, the BAM-file is
#' loaded as a single object, so you must be able to keep at least one complete
#' BAM-file in RAM.
#'
#' @note WARNING: This function is experimental, has not been thoroughly tested,
#'   and will most likely significantly change in upcoming CAGEfightR version.
#'   For comments/question please go to  the CAGEfightR github page.
#'
#' @param input character: Path to input BAM-file
#' @param outputPlus character: Path to output BigWig-file holding CTSSs on the
#'   plus strand.
#' @param outputMinus character: Path to output BigWig-file holding CTSSs on the
#'   minus strand.
#' @param minLength integer: Minimum length of mapped reads.
#' @param ... Additional arguments passed to rtracklayer::import. This will
#'   often include a ScanBamParam
#'
#' @return Number of CTSSs/Tags returned invisibly.
#' @importClassesFrom GenomicAlignments GAlignments
#' @export
#'
#' @examples
#' # TBA
convertBAM2BigWig <- function(input, outputPlus, outputMinus, minLength=1L, ...){
    assert_that(all(checkFileSeries(input)),
                is.character(outputPlus),
                is.character(outputMinus),
                length(input) == length(outputPlus),
                length(input) == length(outputMinus),
                is.count(minLength))

    # Only give warn about improper  extensions
    input_exts <- all(checkExtensions(input, ext="bam"))
    outputPlus_exts <- all(checkExtensions(outputPlus, ext="bw"))
    outputMinus_exts <- all(checkExtensions(outputMinus, ext="bw"))

    if(!input_exts){
        warning("Input BAM files do not all have the proper .bam extension! ",
                "Attempting import anyway...")
    }

    if(!outputPlus_exts | !outputMinus_exts){
        warning("Output BigWig files do not all have the proper .bw extension! ",
                "Attempting import anyway...")
    }

    # Simply loop
    tmp <- BiocParallel::bpmapply(bam2bw,
                                  input,
                                  outputPlus,
                                  outputMinus,
                                  MoreArgs = list(minLength=minLength, ...))

    # Return invisibly
    invisible(tmp)
}

