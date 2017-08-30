library(CAGEfightR)
library(rtracklayer)

# Devtools info
# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next:
# 	* Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package

# Load data
base_dir <- "~/R/CAGEfightR/"
setwd("~/Desktop/mm9_nanotubes/")
design <- read.table("design_matrix.txt", header=TRUE, stringsAsFactors=FALSE)

# Better formated design
design$Samples <- gsub(x=design$Samples, pattern=".fastq", replacement="")
design$BigWigPlus <- paste0("mm9.", design$Samples, ".plus.bw")
design$BigWigMinus <- paste0("mm9.", design$Samples, ".minus.bw")
design$Class <- factor(design$Class)

# Format and sort
design <- subset(design, Name != "C548_2", select=-c(Samples))
design <- design[order(design$Class, design$Name),]
rownames(design) <- design$Name
design <- DataFrame(design)

### Load BigWigs

# Set up BigWig files
bw_plus <- BigWigFileList(design$BigWigPlus)
bw_minus <- BigWigFileList(design$BigWigMinus)
names(bw_plus) <- rownames(design)
names(bw_minus) <- rownames(design)

# Read in CTSS and subset
CTSS <- readCTSS(bwPlus=bw_plus, bwMinus=bw_minus, genome="mm9")
#CTSS <- endoapply(CTSS, subset, seqnames %in% c("chr18", "chr19"))
CTSS <- endoapply(CTSS, subset, (seqnames == "chr18" & start > 72500000) | (seqnames == "chr19" & start > 50000000))

# Write small BigWigs to ext data
setwd(base_dir)
CTSS_plus <- endoapply(CTSS, function(x) splitByStrand(x)$"+")
CTSS_minus <- endoapply(CTSS, function(x) splitByStrand(x)$"-")
fname_plus <- file.path("inst/extdata", design$BigWigPlus)
fname_minus <- file.path("inst/extdata", design$BigWigMinus)

mapply(export.bw, CTSS_plus, fname_plus)
mapply(export.bw, CTSS_minus, fname_minus)

# Save design directly
exampleDesign <- design
devtools::use_data(exampleDesign, overwrite = TRUE)
