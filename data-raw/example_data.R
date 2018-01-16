# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next:
# 	* Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package

# Add full data here

library(CAGEtestR)
library(BSgenome)
library(GenomicFeatures)
library(BiocParallel)
library(magrittr)
register(MulticoreParam(3))

library(BSgenome.Mmusculus.UCSC.mm9)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

### Mouse data

bsg <- BSgenome.Mmusculus.UCSC.mm9
si <- seqinfo(bsg)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

# Some example data
design <- read.table("~/Desktop/mm9_nanotubes/design_matrix.txt", header=TRUE, stringsAsFactors=FALSE)

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

# Set up BigWig files
bw_plus <- BigWigFileList(file.path("~/Desktop/mm9_nanotubes", design$BigWigPlus))
bw_minus <- BigWigFileList(file.path("~/Desktop/mm9_nanotubes", design$BigWigMinus))
names(bw_plus) <- rownames(design)
names(bw_minus) <- rownames(design)

# Rewrite design
design <- transform(design,
										BigWigPlus=paste0("mm9.", Name, ".plus.bw"),
										BigWigMinus=paste0("mm9.", Name, ".minus.bw"),
										batch=factor(c("B", "B", "B", "A", "A", "B", "A", "B", "A", "A", "A")))

#### Example R-files

# Test on full
CTSS <- quantifyCTSSs(bw_plus, bw_minus, genome = si); invisible(gc())
colData(CTSS) <- design

# Subset
ctss <- CTSS #subset(CTSS, (seqnames == "chr18" & pos > 72500000) | (seqnames == "chr19" & pos > 50000000))

# For vignette
TSSs <- quickTSSs(ctss)
enhancers <- quickEnhancers(ctss)

# Add total tags to TSSs
TSSs <- calcTotalTags(TSSs)

# Annotate
#TSSs <- assignTxType(swapRanges(TSSs), txdb)
#enhancers <- assignTxType(enhancers, txdb)

# Subset
#enhancers <- subset(enhancers, txType %in% c("intergenic", "intronic"))

# Genes
GE <- TSSs %>%
	assignGeneID(txdb, swap="thick", downstream=100) %>%
	quantifyGenes(genes="geneID")

#### Subset down

ctss <- calcTotalTags(ctss)
ctss <- subset(ctss, (seqnames == "chr18" & start > 72500000) | (seqnames == "chr19" & start > 50000000))

#### Example BigWig

# Reduce down experiment
bws <- colnames(ctss)
names(bws) <- bws
bws <- lapply(bws, swapScores, object=ctss, outputColumn = "score", inputAssay="counts") %>%
	lapply(rowRanges) #%>%
	#lapply(subset, score > 0)

# Files
ori_dir <- getwd()
CTSS_plus <- endoapply(bws, function(x) splitByStrand(x)$"+")
CTSS_minus <- endoapply(bws, function(x) splitByStrand(x)$"-")
fname_plus <- file.path("inst/extdata", design$BigWigPlus)
fname_minus <- file.path("inst/extdata", design$BigWigMinus)

#### Write to files ####

# Files
mapply(export.bw, CTSS_plus, fname_plus)
mapply(export.bw, CTSS_minus, fname_minus)

# R object
exampleDesign <- design
exampleCTSSs <- ctss
exampleUnidirectional <- subsetByOverlaps(TSSs, ctss)
exampleUnidirectional <- subsetBySupport(exampleUnidirectional, minSamples=2)
exampleUnidirectional$support <- NULL
exampleBidirectional <- subsetByOverlaps(enhancers, ctss)

devtools::use_data(exampleDesign, exampleCTSSs, exampleUnidirectional, exampleBidirectional, exampleGenes, overwrite = TRUE)

#
#
#
# # Merge
# rowRanges(TSSs)$clusterType <- "unidirectional"
# rowRanges(enhancers)$clusterType <- "bidirectional"
# RSE <- combineFeatures(TSSs, enhancers, removeIfOverlapping="object1")
# RSE <- subsetBySupport(RSE, minSamples=4)
#
# # Stats
# pooled <- calcTPM(CTSS)
# pooled <- calcPooled(pooled)
#
# # TCs
# tuned <- tuneTagClustering(pooled, searchMethod="exponential")
# TCs <- clusterUnidirectionally(pooled, pooledCutoff=tuned[which.max(tuned$nTCs), 1])
# TCs <- quantifyClusters(CTSS, TCs); gc()
#
# SE <- calcTPM(SE, outputColumn="TCTags"); gc()
# SE <- subsetBySupport(SE, inputAssay="TPM", unexpressed=1, minSamples=2); gc()
#
# # Annotate
# txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
# SE <- assignTxType(SE, txdb, swap="thick")
# SE <- assignTxID(SE, txdb, swap="thick")
# SE <- assignGeneID(SE, txdb, swap="thick", downstream=100)
#
# # Ggene
# GE <- quantifyGenes(SE, genes="geneID")
# GE <- calcTPM(GE, inputAssay="counts", outputAssay="TPM", outputColumn="geneTags")
#
# exampleCTSS <- rowRanges(CTSS)
# exampleTCs <- SE
# exampleGenes <- GE
#
# devtools::use_data(exampleCTSS, exampleTCs, exampleGenes, overwrite = TRUE)
