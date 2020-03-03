#### Helpers ####

# Check files series
checkFileSeries <- Vectorize(function(x){
    see_if(is.string(x),
              file.exists(x),
              is.readable(x))
})

# Quick extension checker
checkExtensions <- Vectorize(has_extension, vectorize.args="path")

# Bed to BigWig worker
bed2bw <- function(i, o_plus, o_minus, genome){
    # Checks done outside!

    # Attempt import
    message("Attempting to convert file: ", basename(i))
    gr <- import(i, genome=genome)

    # Check file
    tmp <- checkCTSSs(gr)

    # Split
    by_strand <- splitByStrand(gr)

    # Write to file
    export.bw(object=by_strand$`+`, con=o_plus, format="bw")
    export.bw(object=by_strand$`-`, con=o_minus, format="bw")

    # Return
    tmp
}

# Bed to bedgraph worker
bed2bg <- function(i, o_plus, o_minus){
    # Checks done outside!

    # Attempt import
    message("Attempting to convert file: ", basename(i))
    gr <- import(i)

    # Check file
    tmp <- checkCTSSs(gr)

    # Split
    by_strand <- splitByStrand(gr)

    # Write to file
    export.bedGraph(object=by_strand$`+`, con=o_plus, format="bedGraph")
    export.bedGraph(object=by_strand$`-`, con=o_minus, format="bedGraph")

    # Return
    tmp
}

# Bedgraph to BigWig worker
bg2bw <- function(i, o, genome){
    # Attempt import
    message("Attempting to convert file: ", basename(i))
    gr <- import(i, genome=genome)

    # Set dummy strand and check
    strand(gr) <- "+"
    tmp <- checkCTSSs(gr)

    # Write to file (won't use strand)
    export.bw(object=gr, con=o, format="bw")

    # Return
    tmp
}

# Bedgraph to BigWig worker
bw2bg <- function(i, o){
    # Attempt import
    message("Attempting to convert file: ", basename(i))
    gr <- import(i)

    # Set dummy strand and check
    strand(gr) <- "+"
    tmp <- checkCTSSs(gr)

    # Write to file (won't use strand)
    export.bedGraph(object=gr, con=o, format="bedGraph")

    # Return
    tmp
}

# BigWig/bedgraph to BED worker
bwg2bed <- function(i_plus, i_minus, o){
    # Attempt import of both strands
    message("Attempting to convert files: ",
            basename(i_plus), " & ", basename(i_minus))
    gr_plus <- import(i_plus)
    gr_minus <- import(i_minus)

    # Set strands
    strand(gr_plus) <- "+"
    strand(gr_minus) <- "-"

    # Merge and check
    gr <- c(gr_plus, gr_minus)
    tmp <- checkCTSSs(gr)

    # Write to file (won't use strand)
    export.bed(object=gr, con=o, format="bed")

    # Return
    tmp
}

#### Main ####

#' Convert GRanges with scores to GPos
#'
#' Converts a GRanges to a GPos, correctly expanding the score column. This is
#' useful is nearby CTSSs with the same count are grouped in the same range (see
#' example).
#'
#' @param object GRanges object with a score column
#'
#' @return GPos with score column
#' @export
#'
#' @examples
#' # Example GRanges
#' gr <- GRanges(Rle(c("chr2", "chr2", "chr3", "chr4")),
#'               IRanges(start=c(1, 10, 5, 3),
#'               end=c(5L, 10L, 5L, 4L)),
#'               strand="+",
#'               score=c(2, 1, 3, 11))
#'
#' # Expand to proper GPos / CTSS format:
#' gp <- convertGRanges2GPos(gr)
#'
#' # Double check that the total number of counts remains the same
#' stopifnot(sum(score(gr) * width(gr)) == sum(score(gp)))
convertGRanges2GPos <- function(object){
    # Pre-checks
    assert_that(methods::is(object, "GRanges"),
                isDisjoint(object),
                !is.null(score(object)),
                is.numeric(score(object)))

    # Convert to GPos and expand score
    o <- GPos(object)
    score(o) <- rep(score(object), width(object))

    # Return
    o
}

#' Convert CTSSs stored in different file formats.
#'
#' Collection of functions for converting CTSSs/CTSSs-like data stored in
#' BigWig, bedGraph or BED file formats. BigWig and bedGraph files use a file
#' for each strand, while BED-files stores both strands in a single file. As
#' BigWig files stores info about the chromosome lenghts, conversion from
#' bedGraph/BED to BigWig requires a genome. Note that CAGEfightR will only
#' import BigWig or bedGraph files!
#'
#' @note These functions will warn if input files do not have the correct
#'   extensions (.bw, .bedGraph, .bed), but otherwise simply pass input to
#'   rtracklayer::import. This makes them able to handle compressed files (like
#'   .gz). The same applies to the genome argument, which can also be the name
#'   of a UCSC genome.
#'
#' @param input charater: Path to input files holding CTSSs on both strands.
#' @param inputPlus character: Path to input files holding CTSSs on plus strand.
#' @param inputMinus character: Path to input files holding CTSSs on minus strand.
#' @param output charater: Path to output files holding CTSSs on both strands.
#' @param outputPlus character: Path to output files holding CTSSs on plus strand.
#' @param outputMinus character: Path to output files holding CTSSs on minus strand.
#' @param genome Seqinfo or character: Genome info passed to rtracklayer::import  (see note).
#'
#' @return TRUE returned invisibly if conversion(s) was succesful, otherwise an error is raised.
#' @export
#'
#' @examples
#' # Find paths to BigWig files
#' data('exampleDesign')
#' bw_plus <- system.file('extdata', exampleDesign$BigWigPlus,
#'                        package = 'CAGEfightR')
#' bw_minus <- system.file('extdata', exampleDesign$BigWigMinus,
#'                         package = 'CAGEfightR')
#'
#' # Designate paths to new files
#' n_samples <- length(bw_plus)
#' beds <- replicate(n=n_samples, tempfile(fileext=".bed"))
#' bg_plus <- replicate(n=n_samples, tempfile(fileext="_plus.bedGraph"))
#' bg_minus <- replicate(n=n_samples, tempfile(fileext="_minus.bedGraph"))
#' conv_plus <- replicate(n=n_samples, tempfile(fileext="_plus.bw"))
#' conv_minus <- replicate(n=n_samples, tempfile(fileext="_minus.bw"))
#'
#' # Convert BigWig to BED
#' convertBigWig2BED(inputPlus=bw_plus,
#'                   inputMinus=bw_minus,
#'                   output=beds)
#'
#' # Convert BED to bedGraph
#' convertBED2BedGraph(input=beds,
#'                     outputPlus=bg_plus,
#'                     outputMinus=bg_minus)
#'
#' # Convert BED to bedGraph
#' mm9 <- SeqinfoForUCSCGenome("mm9")
#' convertBED2BigWig(input=beds,
#'                   outputPlus=conv_plus,
#'                   outputMinus=conv_minus,
#'                   genome=mm9)
#'
#' # Check it's still the same data
#' x <- import(bw_plus[1])
#' y <- import(bg_plus[1])
#' z <- import(conv_plus[1])
#' all(x == y)
#' all(x == z)
#' sum(score(x)) ==  sum(score(y))
#' sum(score(x)) ==  sum(score(z))
convertBED2BigWig <- function(input, outputPlus, outputMinus, genome){
    assert_that(all(checkFileSeries(input)),
                is.character(outputPlus),
                is.character(outputMinus),
                length(input) == length(outputPlus),
                length(input) == length(outputMinus))

    # Only give warn about improper  extensions
    input_exts <- all(checkExtensions(input, ext="bed"))
    outputPlus_exts <- all(checkExtensions(outputPlus, ext="bw"))
    outputMinus_exts <- all(checkExtensions(outputMinus, ext="bw"))

    if(!input_exts){
        warning("Input BED files do not all have the proper .bed extension! ",
                "Attempting import anyway...")
    }

    if(!outputPlus_exts | !outputMinus_exts){
        warning("Output BigWig files do not all have the proper .bw extension! ",
                "Attempting import anyway...")
    }

    # Check genome as well
    if(!methods::is(genome, "Seqinfo")){
        warning("Supplied genome is not a Seqinfo object! ",
                "Attempting import anyway...")
    }

    # Simply loop
    tmp <- BiocParallel::bpmapply(bed2bw,
                                  input,
                                  outputPlus,
                                  outputMinus,
                                  MoreArgs = list(genome = genome))

    # Return invisibly
    invisible(tmp)
}

#' @rdname convertBED2BigWig
#' @export
convertBED2BedGraph <- function(input, outputPlus, outputMinus){
    assert_that(all(checkFileSeries(input)),
                is.character(outputPlus),
                is.character(outputMinus),
                length(input) == length(outputPlus),
                length(input) == length(outputMinus))

    # Only give warn about improper  extensions
    input_exts <- all(checkExtensions(input, ext="bed"))
    outputPlus_exts <- all(checkExtensions(outputPlus, ext="bedGraph"))
    outputMinus_exts <- all(checkExtensions(outputMinus, ext="bedGraph"))

    if(!input_exts){
        warning("Input BED files do not all have the proper .bed extension! ",
                "Attempting import anyway...")
    }

    if(!outputPlus_exts | !outputMinus_exts){
        warning("Output bedGraph files do not all have the proper .bedGraph extension! ",
                "Attempting import anyway...")
    }

    # Simply loop
    tmp <- BiocParallel::bpmapply(bed2bg, input, outputPlus, outputMinus)

    # Return invisibly
    invisible(tmp)
}

#' @rdname convertBED2BigWig
#' @export
convertBedGraph2BigWig <- function(input, output, genome){
    assert_that(all(checkFileSeries(input)),
                is.character(output),
                length(input) == length(output))

    # Only give warn about improper  extensions
    input_exts <- all(checkExtensions(input, ext="bedgraph"))
    output_exts <- all(checkExtensions(output, ext="bw"))

    if(!input_exts){
        warning("Input bedGraph files do not all have the proper .bedGraph extension! ",
                "Attempting import anyway...")
    }

    if(!output_exts){
        warning("Output BigWig files do not all have the proper .bw extension! ",
                "Attempting import anyway...")
    }

    # Check genome as well
    if(!methods::is(genome, "Seqinfo")){
        warning("Supplied genome is not a Seqinfo object! ",
                "Attempting import anyway...")
    }


    # Simply loop
    tmp <-  BiocParallel::bpmapply(bg2bw,
                                   input,
                                   output,
                                   MoreArgs = list(genome = genome))

    # Return invisibly
    invisible(tmp)
}

#' @rdname convertBED2BigWig
#' @export
convertBigWig2BedGraph <- function(input, output){
    assert_that(all(checkFileSeries(input)),
                is.character(output),
                length(input) == length(output))

    # Only give warn about improper  extensions
    input_exts <- all(checkExtensions(input, ext="bw"))
    output_exts <- all(checkExtensions(output, ext="bedGraph"))

    if(!input_exts){
        warning("Input BigWig files do not all have the proper .bw extension! ",
                "Attempting import anyway...")
    }

    if(!output_exts){
        warning("Output bedGraph files do not all have the proper .bedGraph extension! ",
                "Attempting import anyway...")
    }

    # Simply loop
    tmp <-  BiocParallel::bpmapply(bw2bg,
                                   input,
                                   output)

    # Return invisibly
    invisible(tmp)
}

#' @rdname convertBED2BigWig
#' @export
convertBigWig2BED <- function(inputPlus, inputMinus, output){
    assert_that(all(checkFileSeries(inputPlus)),
                all(checkFileSeries(inputMinus)),
                is.character(output),
                length(inputPlus) == length(output),
                length(inputMinus) == length(output))

    # Only give warn about improper  extensions
    inputPlus_exts <- all(checkExtensions(inputPlus, ext="bw"))
    inputMinus_exts <- all(checkExtensions(inputMinus, ext="bw"))
    output_exts <- all(checkExtensions(output, ext="bed"))

    if(!inputPlus_exts | !inputMinus_exts){
        warning("Input BigWig files do not all have the proper .bw extension! ",
                "Attempting import anyway...")
    }

    if(!output_exts){
        warning("Output BED files do not all have the proper .bed extension! ",
                "Attempting import anyway...")
    }

    # Simply loop
    tmp <-  BiocParallel::bpmapply(bwg2bed,
                                   inputPlus,
                                   inputMinus,
                                   output)

    # Return invisibly
    invisible(tmp)
}

#' @rdname convertBED2BigWig
#' @export
convertBedGraph2BED <- function(inputPlus, inputMinus, output){
    assert_that(all(checkFileSeries(inputPlus)),
                all(checkFileSeries(inputMinus)),
                is.character(output),
                length(inputPlus) == length(output),
                length(inputMinus) == length(output))

    # Only give warn about improper  extensions
    inputPlus_exts <- all(checkExtensions(inputPlus, ext="bedGraph"))
    inputMinus_exts <- all(checkExtensions(inputMinus, ext="bedGraph"))
    output_exts <- all(checkExtensions(output, ext="bed"))

    if(!inputPlus_exts | !inputMinus_exts){
        warning("Input BigWig files do not all have the proper .bedGraph extension! ",
                "Attempting import anyway...")
    }

    if(!output_exts){
        warning("Output BED files do not all have the proper .bed extension! ",
                "Attempting import anyway...")
    }

    # Simply loop
    tmp <-  BiocParallel::bpmapply(bwg2bed,
                                   inputPlus,
                                   inputMinus,
                                   output)

    # Return invisibly
    invisible(tmp)
}
