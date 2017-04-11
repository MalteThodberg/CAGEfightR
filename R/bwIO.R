#' Count the number of different genomes in a series of BigWig-files.
#'
#' Find all unique genomes in in a series of BigWig-files.
#'
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data
#'

#' @return SimpleList holding Seqinfo-objects.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family IO functions
#' @import S4Vectors GenomeInfoDb
#' @export
uniqueBigWigGenomes <- function(bwPlus, bwMinus){
	# Get seqinfo
	seqInfoPlus <- lapply(bwPlus, seqinfo)
	seqInfoMinus <- lapply(bwMinus, seqinfo)

	# Check single seqinfo
	seqInfo <- unique(c(seqInfoPlus, seqInfoMinus))

	# Return the genome info
	SimpleList(seqInfo)
}

#' Find a common genome for a series of BigWig files.
#'
#' Finds a common genome for a series of BigWig-files, either using only levels present in all files (intersect) or in any file (union).
#'
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param method character: Either "intersect" or "union".
#'
#' @return Seqinfo-object.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family IO functions
#' @import S4Vectors GenomeInfoDb
#' @export
commonBigWigGenome <- function(bwPlus, bwMinus, method="intersect"){
	# Get seqinfo
	seqInfoPlus <- lapply(bwPlus, seqinfo)
	seqInfoMinus <- lapply(bwMinus, seqinfo)

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

#' Check if BigWig-files are compatible with a genome.
#'
#' Given a genome, checks whether a series of BigWig-files are compatible by checking if all common seqlevels have equal seqlengths.
#'
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param seqInfo Seqinfo: Genome information.
#'
#' @return None: Throws an error if the supplied genome is incompabtible.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @import S4Vectors GenomeInfoDb
#' @family IO functions
#' @export
checkBigWigGenomes <- function(bwPlus, bwMinus, seqInfo){
	# Get seqinfo
	seqInfoPlus <- lapply(bwPlus, seqinfo)
	seqInfoMinus <- lapply(bwMinus, seqinfo)

	# Unique seqinfos
	seqInfos <- unique(c(seqInfoPlus, seqInfoMinus))

	# Check if error is produces
	lapply(seqInfos, merge, y=seqInfo)
}

#' Import a pair of BigWig-files.
#'
#' Reads in CTSS-data stored in a BigWig-file for each strand.
#'
#' @param fnamePlus BigWigFile: BigWig file with plus-strand CTSS data.
#' @param fnameMinus BigWigFile: BigWig file with minus-strand CTSS data.
#' @param seqInfo Seqinfo or NULL: If not NULL, use the supplied Seqinfo as genome.
#'
#' @return GRanges with CTSS info as the score column.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family IO functions
#' @import S4Vectors IRanges GenomicRanges GenomeInfoDb rtracklayer
#' @export
fastReadBigWig <- function(fnamePlus, fnameMinus, seqInfo=NULL){
	# Read files
	countsPlus <- import(fnamePlus)
	strand(countsPlus) <- "+"
	countsMinus <- import(fnameMinus)
	strand(countsMinus) <- "-"

	# To GR
	gr <- c(countsPlus, countsMinus)
	rm(countsPlus, countsMinus)

	# Set new sequence info if necessary
	if(!is.null(seqInfo)){
		# Call error if genomes can't me merged with filename in question!
		chr_map <- match(seqlevels(seqInfo), seqlevels(gr))
		seqinfo(gr, new2old=chr_map, force=TRUE) <- seqInfo
	}

	# Seqinfo
	gr
}

#' Import CTSS data from BigWig-files.
#'
#' This function reads a series of BigWig-files into R as a GRanges. Multiple checks are performed to ensure the genomes of the different BigWig-files are compatible.
#'
#' @param bwPlus BigWigFileList: BigWig files with plus-strand CTSS data.
#' @param bwMinus BigWigFileList: BigWig files with minus-strand CTSS data.
#' @param genome Seqinfo or NULL: If not NULL, use the supplied Seqinfo as genome.
#' @param outFormat character: Output as either a "GRangesList" or "SimpleList".
#' @param biocParallel BiocParallelParam: Settings for parallel backend.
#' @param ... additional arguments passed to commonBigWigGenome used when attempting to merge genomes..
#'
#' @return Either a GRangesList or a SimpleList holding GRanges with CTSS data as the score column.
#' @examples
#' # ADD_EXAMPLES_HERE
#' @family IO functions
#' @import S4Vectors IRanges GenomicRanges GenomeInfoDb rtracklayer
#' @export
readCTSS <- function(bwPlus, bwMinus, genome=NULL, outFormat="GRangesList", biocParallel=bpparam(), ...){
	# Check genome if not provided
	if(is.null(genome)){
		# Check if all files have the same genome
		uniqueSeqInfos <- uniqueBigWigGenomes(bwPlus=bwPlus, bwMinus=bwMinus)

		if(length(uniqueSeqInfos) != 1){
			message("BigWig files do not have the same genome! Attempting to find common reference...")
			genome <- commonBigWigGenome(bwPlus=bwPlus, bwMinus=bwMinus, ...)
		}else{
			message("BigWig files have the same genome...")
			genome <- uniqueSeqInfos[[1]]
		}

	}else if(class(genome) == "character"){
		message("Retrieving genome from BSGenome...")
		genome <- SeqinfoForBSGenome(genome)
		checkBigWigGenomes(bwPlus=bwPlus, bwMinus=bwMinus, seqInfo=genome)

	}else if(class(genome) == "Seqinfo"){
		message("Using supplied genome...")
		checkBigWigGenomes(bwPlus=bwPlus, bwMinus=bwMinus, seqInfo=genome)

	}else{
		stop("Genome must be a Seqinfo-object or character!")

	}

	# Print used genome
	print(genome)

	#Read into R
	message("Reading in CTSS data from BigWig-files...")
	o <- bpmapply(FUN=fastReadBigWig, bwPlus, bwMinus,
								MoreArgs=list(seqInfo=genome),
								BPPARAM=biocParallel)

	if(outFormat == "GRangesList"){
		# Coerce to GRangesList
		message("Preparing final GRangesList...")
		o <- GRangesList(o)
	}else if(outFormat == "SimpleList"){
		message("Return list of GRanges...")
		o <- SimpleList(o)
	}

	# Return
	o
}
