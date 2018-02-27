## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----bioconductor, eval=FALSE----------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("CAGEfightR")

## ----library, results='hide', message=FALSE--------------------------------
library(CAGEfightR)

## ----github, eval=FALSE----------------------------------------------------
#  devtools::install_github("MalteThodberg/CAGEfightR")

## ----citation, eval=FALSE--------------------------------------------------
#  citation("CAGEfightR")

## ----BigWig_files, results="hide", tidy=TRUE-------------------------------
# Load the example data
data("exampleDesign")
head(exampleDesign)

# Locate files on your harddrive
bw_plus <- system.file("extdata", exampleDesign$BigWigPlus, 
											 package = "CAGEfightR")
bw_minus <- system.file("extdata", exampleDesign$BigWigMinus, 
												package = "CAGEfightR")

# Create two named BigWigFileList-objects:
bw_plus <- BigWigFileList(bw_plus)
bw_minus <- BigWigFileList(bw_minus)
names(bw_plus) <- exampleDesign$Name
names(bw_minus) <- exampleDesign$Name

## ----quickCTSSs, tidy=TRUE-------------------------------------------------
# Get genome information
genomeInfo <- SeqinfoForBSGenome("mm9")

# Quantify CTSSs
CTSSs <- quantifyCTSSs(plusStrand=bw_plus, 
											 minusStrand=bw_minus, 
											 design=exampleDesign, 
											 genome=genomeInfo)

## ----quickTSSs, tidy=TRUE--------------------------------------------------
TSSs <- quickTSSs(CTSSs)

## ----quickEnhancers, tidy=TRUE---------------------------------------------
enhancers <- quickEnhancers(CTSSs)

## ----quickAnnotate, tidy=TRUE----------------------------------------------
# Use the built in annotation for mm9
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

# Annotate both TSSs and enhancers
TSSs <- assignTxType(TSSs, txModels=txdb)
enhancers <- assignTxType(enhancers, txModels=txdb)

## ----quickFilter, tidy=TRUE------------------------------------------------
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))

## ----quickCombine, tidy=TRUE-----------------------------------------------
# Add an identifier column
rowRanges(TSSs)$clusterType <- "TSS"
rowRanges(enhancers)$clusterType <- "enhancer"

# Combine TSSs and enhancers, discarding TSSs if they overlap enhancers
RSE <- combineClusters(TSSs, enhancers, removeIfOverlapping="object1")

## ----quickSupport, tidy=TRUE-----------------------------------------------
# Only keep features with more than 0 counts in more than 1 sample:
RSE <- subsetBySupport(RSE, 
											 inputAssay = "counts", 
											 outputColumn = "support", 
											 unexpressed = 0, 
											 minSamples = 1)

## ----exampleCTSSs, tidy=TRUE-----------------------------------------------
data(exampleCTSSs)
exampleCTSSs

## ----dgCMatrix, tidy=TRUE--------------------------------------------------
head(assay(exampleCTSSs))

## ----GPos, tidy=TRUE-------------------------------------------------------
rowRanges(exampleCTSSs)

## ----calc, tidy=TRUE-------------------------------------------------------
exampleCTSSs <- calcTPM(exampleCTSSs, inputAssay="counts", outputAssay="TPM", outputColumn="subsetTags")

## ----libSizes, tidy=TRUE---------------------------------------------------
# Library sizes
colData(exampleCTSSs)

# TPM values
head(assay(exampleCTSSs, "TPM"))

## ----preCalcTPM, tidy=TRUE-------------------------------------------------
exampleCTSSs <- calcTPM(exampleCTSSs, inputAssay="counts", totalTags="totalTags", outputAssay="TPM")
head(assay(exampleCTSSs, "TPM"))

## ----pooled, tidy=TRUE-----------------------------------------------------
# Library sizes
exampleCTSSs <- calcPooled(exampleCTSSs, inputAssay="TPM")
rowRanges(exampleCTSSs)

## ----support, tidy=TRUE----------------------------------------------------
# Count number of samples with MORE ( > ) than 0 counts:
exampleCTSSs <- calcSupport(exampleCTSSs, inputAssay="counts", outputColumn="support", unexpressed=0)
table(rowRanges(exampleCTSSs)$support)

## ----subset, tidy=TRUE-----------------------------------------------------
supportedCTSSs <- subset(exampleCTSSs, support > 1)
supportedCTSSs <- calcTPM(supportedCTSSs, totalTags="totalTags")
supportedCTSSs <- calcPooled(supportedCTSSs)

## ----naiveTC, tidy=TRUE----------------------------------------------------
naive_TCs <- clusterUnidirectionally(exampleCTSSs, pooledCutoff=0, mergeDist=20)

## ----tuning, tidy=TRUE-----------------------------------------------------
tuned <- tuneTagClustering(exampleCTSSs, mergeDist=20)
tuned

## ----tagClustering, tidy=TRUE----------------------------------------------
optimalCutoff <- tuned[which.max(tuned$nTCs),1]
TCs <- clusterUnidirectionally(exampleCTSSs, pooledCutoff=optimalCutoff, mergeDist=20)

## ----TCanatomy, tidy=TRUE--------------------------------------------------
TCs

## ----bidirectionalClustering, tidy=TRUE------------------------------------
enhancers <- clusterBidirectionally(exampleCTSSs, balanceThreshold=0.95)

## ----enhancerAnatomy, tidy=TRUE--------------------------------------------
enhancers

## ----bidirectionality, tidy=TRUE-------------------------------------------
# Calculate number of bidirectional samples
enhancers <- calcBidirectionality(enhancers, samples=exampleCTSSs)

# Summarize
table(enhancers$bidirectionality)

## ----subsetBidirectionality, tidy=TRUE-------------------------------------
enhancers <- subset(enhancers, bidirectionality > 0)

## ----exampleClusters, tidy=TRUE--------------------------------------------
# Load the example datasets
data(exampleCTSSs)
data(exampleUnidirectional)
data(exampleBidirectional)

## ----quantifyClusters, tidy=TRUE-------------------------------------------
requantified_TSSs <- quantifyClusters(exampleCTSSs, 
																			clusters=rowRanges(exampleUnidirectional), 
																			inputAssay="counts")
requantified_enhancers <- quantifyClusters(exampleCTSSs, 
																					 clusters=rowRanges(exampleBidirectional), 
																					 inputAssay="counts")

## ----supportOnCounts, tidy=TRUE--------------------------------------------
# Only keep enhancers expressed in more than one sample
supported_enhancers <- subsetBySupport(exampleBidirectional, 
																			 inputAssay="counts", 
																			 unexpressed=0, 
																			 minSamples=1)

## ----supportOnTPM, tidy=TRUE-----------------------------------------------
# Calculate TPM using pre-calculated total tags:
exampleUnidirectional <- calcTPM(exampleUnidirectional, totalTags = "totalTags")

# Only TSSs expressed at more than 1 TPM in more than 2 samples
exampleUnidirectional <- subsetBySupport(exampleUnidirectional, 
																				 inputAssay="TPM", 
																				 unexpressed=1, 
																				 minSamples=2)

## ----txdb, tidy=TRUE-------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
txdb

## ----assignTxID, tidy=TRUE-------------------------------------------------
exampleUnidirectional <- assignTxID(exampleUnidirectional, 
																		txModels=txdb, 
																		outputColumn="txID")

## ----multipleTxs, tidy=TRUE------------------------------------------------
rowRanges(exampleUnidirectional)[5:6]

## ----assignTxType, tidy=TRUE-----------------------------------------------
exampleUnidirectional <- assignTxType(exampleUnidirectional, 
																			txModels=txdb, 
																			outputColumn="txTyp")

## ----swappedTxType, tidy=TRUE----------------------------------------------
exampleUnidirectional <- assignTxType(exampleUnidirectional, 
																			txModels=txdb, 
																			outputColumn="peakTxType", 
																			swap="thick")

## ----enhancerTxType, tidy=TRUE---------------------------------------------
# Annotate with TxTypes
exampleBidirectional <- assignTxType(exampleBidirectional, 
																		 txModels=txdb, 
																		 outputColumn="txType")

# Only keep intronic and intergenic enhancers
exampleBidirectional <- subset(exampleBidirectional, 
															 txType %in% c("intron", "intergenic"))

## ----calcIQR, tidy=TRUE----------------------------------------------------
# Recalculate pooled signal
exampleCTSSs <- calcTPM(exampleCTSSs, totalTags = "totalTags")
exampleCTSSs <- calcPooled(exampleCTSSs)

# Calculate shape
exampleUnidirectional <- calcShape(exampleUnidirectional, 
																	 pooled=exampleCTSSs, 
																	 outputColumn = "IQR", 
																	 shapeFunction = shapeIQR, 
																	 lower=0.25, upper=0.75)

## ----histIQR, tidy=TRUE----------------------------------------------------
hist(rowRanges(exampleUnidirectional)$IQR, 
		 breaks=max(rowRanges(exampleUnidirectional)$IQR), 
		 xlim=c(0,100), xlab = "IQR", 
		 col="red")

## ----customShape, tidy=TRUE------------------------------------------------
# Write a function that quantifies the lean of a TSS
shapeLean <- function(x, direction){
	# Coerce to normal vector
	x <- as.vector(x)
	
	# Find highest position:
	i <- which.max(x)
	
	# Calculate sum
	i_total <- sum(x)
	
	# Calculate lean fraction
	if(direction == "right"){
		i_lean <- sum(x[i:length(x)])
	}else if(direction == "left"){
		i_lean <- sum(x[1:i])
	}else{
		stop("direction must be left or right!")
	}
	
	# Calculate lean
	o <- i_lean / i_total
	
	# Return
	o
}

# Calculate lean statistics, 
# additional arguments can be passed to calcShape via "..."
exampleUnidirectional <- calcShape(exampleUnidirectional, exampleCTSSs, 
					outputColumn = "leanRight", shapeFunction=shapeLean, 
					direction="right")

exampleUnidirectional <- calcShape(exampleUnidirectional, exampleCTSSs, 
					outputColumn = "leanLeft", shapeFunction=shapeLean, 
					direction="left")

## ----geneSetup, tidy=TRUE--------------------------------------------------
# Load example TSS 
data(exampleUnidirectional)

# Keep only TCs expressed at more than 1 TPM in more than 2 samples:
exampleUnidirectional <- calcTPM(exampleUnidirectional, totalTags = "totalTags")
exampleUnidirectional <- subsetBySupport(exampleUnidirectional, 
																				 inputAssay="TPM", 
																				 unexpressed=1, 
																				 minSamples=2)

# Use the Bioconductor mm9 UCSC TxXb 
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

## ----geneModels, tidy=TRUE-------------------------------------------------
exampleUnidirectional <- assignGeneID(exampleUnidirectional, 
																			geneModels=txdb, 
																			outputColumn="geneID")

## ----symbols, tidy=TRUE----------------------------------------------------
# Use Bioconductor OrgDb package
library(org.Mm.eg.db)
odb <- org.Mm.eg.db

# Match IDs to symbols
symbols <- mapIds(odb, 
									keys=rowRanges(exampleUnidirectional)$geneID, 
									keytype="ENTREZID",
									column="SYMBOL")

# Add to object
rowRanges(exampleUnidirectional)$symbol <- as.character(symbols)

## ----assignMissing, tidy=TRUE----------------------------------------------
exampleUnidirectional <- assignMissingID(exampleUnidirectional, 
																				 outputColumn="symbol")

## ----quantifyGenes, tidy=TRUE----------------------------------------------
genelevel <- quantifyGenes(exampleUnidirectional, genes="geneID", inputAssay="counts")

## ----geneGRangesList, tidy=TRUE--------------------------------------------
rowRanges(genelevel)

## ----rowData, tidy=TRUE----------------------------------------------------
rowData(genelevel)

## ----calcComposition, tidy=TRUE--------------------------------------------
# Remove TSSs not belonging to any genes
intragenicTSSs <- subset(exampleUnidirectional, !is.na(geneID))

# Calculate composition: The number of samples expressing TSSs above 10% of the total gene expression.
intragenicTSSs <- calcComposition(intragenicTSSs, 
																	outputColumn="composition", 
																	unexpressed=0.1, 
																	genes="geneID")

# Overview of composition values:
table(rowRanges(intragenicTSSs)$composition)

## ----subsetComposition, tidy=TRUE------------------------------------------
# Remove TSSs with a composition score less than 3
intragenicTSSs <- subset(intragenicTSSs, composition > 2)

## ----GViz, tidy=TRUE-------------------------------------------------------
library(Gviz)

data("exampleCTSSs")
data("exampleUnidirectional")
data("exampleBidirectional")

## ----builtCTSStrack, tidy=TRUE---------------------------------------------
# Calculate pooled CTSSs
exampleCTSSs <- calcTPM(exampleCTSSs, totalTags="totalTags")
exampleCTSSs <- calcPooled(exampleCTSSs)

# Built the track
pooled_track <- trackCTSS(exampleCTSSs, name="CTSSs")

## ----plotCTSStrack, tidy=TRUE----------------------------------------------
# Plot the main TSS of the myob gene
plotTracks(pooled_track, from=74601950, to=74602100, chromosome="chr18")

## ----clusterTrack, tidy=TRUE-----------------------------------------------
# Remove columns
exampleUnidirectional$totalTags <- NULL
exampleBidirectional$totalTags <- NULL

# Combine TSSs and enhancers, discarding TSSs if they overlap enhancers
CAGEclusters <- combineClusters(object1=exampleUnidirectional, object2=exampleBidirectional, removeIfOverlapping="object1")

# Only keep features with more than 0 counts in more than 2 samples
CAGEclusters <- subsetBySupport(CAGEclusters, inputAssay = "counts", unexpressed=0, minSamples=2)

# Built track
cluster_track <- trackClusters(CAGEclusters, name="Clusters", col=NA)

## ----plotClustertrack, tidy=TRUE-------------------------------------------
# Plot the main TSS of the myob gene
plotTracks(cluster_track, from=74601950, to=74602100, chromosome="chr18")

## ----prettybrowser, tidy=TRUE----------------------------------------------
# Genome axis tracks
axis_track <- GenomeAxisTrack()

# Transcript model track
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
tx_track <- GeneRegionTrack(txdb, name="Gene Models", col=NA, fill="bisque4", shape="arrow")

# Merge all tracks
all_tracks <- list(axis_track, tx_track, cluster_track, pooled_track)

# Plot the Hdhd2 gene
plotTracks(all_tracks, from=77182700, to=77184000, chromosome="chr18")

# Plot an intergenic enhancer
plotTracks(all_tracks, from=53499000, to=53499600, chromosome="chr19")

## ----parallel, tidy=TRUE, eval=FALSE---------------------------------------
#  library(BiocParallel)
#  
#  # Setup for parallel execustion with 3 threads:
#  register(MulticoreParam(workers=3))
#  
#  # Disable parallel execution
#  register(SerialParam())

## ----sessionInfo, tidy=TRUE------------------------------------------------
sessionInfo()

