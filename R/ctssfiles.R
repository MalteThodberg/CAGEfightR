# library(coRe)
# setup_R()
#
#
# ### Settings
#
# window_width <- 200
#
# ### Paths and files
#
# # Directories
# project_dir <- file.path(run_loc(),
# 												 "seqdata/sandelin/projects/OLAF_YEAST")
#
# rstudio_dir <- file.path(run_loc(),
# 												 "home/nzl922/rstudio_projects/pombe")
#
# personal_dir <- file.path(run_loc(),
# 													"seqdata/sandelin/people/nzl922/pombe")
#
# # Ctss file
# ctss_fnames <- list.files(file.path(project_dir, "CAGE_r26_final/MAPPING/ctss_files"),
# 												 full.names=TRUE)
#
# ctss_samples <- map_chr(ctss_fnames, basename)
# names(ctss_fnames) <- ctss_samples
#
# # Load files
# ctss_grs <- GRangesList(lapply(ctss_fnames, import))
#
# # Calculate TPM
# calcTPM <- function(gr){
# 	# Calculate total sum in millions
# 	s <- sum(score(gr)) / 1e6
#
# 	# Calculate tpm
# 	tpm <- score(gr) / s
#
# 	# Append to gr and return
# 	gr$tpm <- tpm
#
# 	gr
# }
#
# ctss_tpm <- endoapply(ctss_grs, calcTPM)
#
# # Split into strand to calculate coverate
# ctss_plus <- endoapply(ctss_tpm, function(x) subset(x, strand="+"))
# ctss_minus <- endoapply(ctss_tpm, function(x) subset(x, strand="-"))
#
# # Calculate Strand-wise coverage
# coverage_plus <- coverage(ctss_plus, weight="tpm")
# coverage_minus <- coverage(ctss_minus, weight="tpm")
#
# # Peak calling: Find peaks
# TPM_min <- 1
# peaks_plus <- slice(coverage_plus, lower=TPM_min, upper=Inf, rangesOnly=TRUE)
# peaks_minus <- slice(coverage_minus, lower=TPM_min, upper=Inf, rangesOnly=TRUE)
#
# # Peak calling: Merge nearby peaks
# reduced_plus <- reduce(peaks_plus, min.gapwidth=25)
# reduced_minus <- reduce(peaks_minus, min.gapwidth=25)
#
# # Peaks stats: summed TPM
# views_plus <- Views(coverage_plus, reduced_plus)
# views_minus <- Views(coverage_minus, reduced_minus)
#
# sum_plus <- unlist(viewSums(views_plus))
# sum_minus <- unlist(viewSums(views_minus))
#
# # Peaks stats: Dominant TSS
# peaks_plus <- unlist(viewWhichMaxs(views_plus))
# peaks_minus <- unlist(viewWhichMaxs(views_minus))
#
# ranges_plus <- IRanges(start=unlist(start(reduced_plus)) + peaks_plus - 1, width=1)
# ranges_minus <- IRanges(start=unlist(start(reduced_minus)) + peaks_minus - 1, width=1)
#
# ranges_plus <- IRanges(start=peaks_plus, width=1)
# ranges_minus <- IRanges(start=peaks_minus, width=1)
#
# # Peaks shapes: IQR
# vectors_plus <- viewApply(views_plus, as.vector) %>% unlist
# dist_plus <- lapply(vectors_plus, function(x) rep(seq_len(length(x)), x))
# iqr_plus <- sapply(dist_plus, IQR)
# skewness_plus <- sapply(dist_plus, skewness)
# kurtosis_plus <- sapply(dist_plus, kurtosis)
#
# vectors_minus <- viewApply(views_minus, as.vector) %>% unlist
# dist_minus <- lapply(vectors_minus, function(x) rep(seq_len(length(x)), x))
# iqr_minus <- sapply(dist_minus, IQR)
# skewness_minus <- sapply(dist_minus, skewness)
# kurtosis_minus <- sapply(dist_minus, kurtosis)
#
# TCs_plus <- GRanges(reduced_plus, strand="+", peak=ranges_plus, score=sum_plus,
# 										iqr=iqr_plus, skewness=skewness_plus, kurtosis=kurtosis_plus)
# TCs_minus <- GRanges(reduced_minus, strand="+", peak=ranges_minus, score=sum_minus,
# 										iqr=iqr_minus, skewness=skewness_minus, kurtosis=kurtosis_minus)
#
# TCs_both <- c(TCs_plus, TCs_minus)
#
# # Function for counting
# countTCs <- function(features, reads, ignore.strand=TRUE, inter.feature=FALSE){
# 	findOverlaps(query=features, subject=reads) %>%
# 		as.data.frame %>%
# 		as_tibble() %>%
# 		mutate(score=score(reads)[subjectHits]) %>%
# 		group_by(queryHits) %>%
# 		summarise(count=sum(score)) %>%
# 		use_series(count)
# }
#
# # Function for counting
# countTCs <- function(features, reads, ignore.strand, inter.feature){
# 	fo <- findOverlaps(query=features, subject=reads) %>%
# 		as.data.frame %>%
# 		as_tibble()
# 		mutate(score=score(reads)[subjectHits]) %>%
# 		group_by(queryHits) %>%
# 		summarise(count=sum(score)) %>%
# 		use_series(count)
# }
#
# #summarizeOverlaps(features=TCs_both, reads=ctss_grs, mode=countTCs)
#
# gr <- ctss_tpm[[1]]
#
#
# fo <- findOverlaps(query=TCs_both, subject=gr) %>%
# 	as.data.frame %>%
# 	as_tibble() %>%
# 	mutate(score=score(gr)[subjectHits]) %>%
# 	group_by(queryHits) %>%
# 	summarise(count=sum(score))
#
#
#
# # Find overlaps
# fo <- findOverlaps(query=TCs_both, subject=gr)
#
# # Dataframe to be merged
# df <- data_frame(tcIndex=queryHits(fo),
# 					 ctssIndex=subjectHits(fo)) %>%
# 	mutate(ctssScores=score(gr)[ctssIndex]) %>%
# 	group_by(tcIndex) %>%
# 	summarise(tcCount=sum(ctssScores))
#
# # Allocate zero count vector
# tagCounts1 <- rep(0, length(TCs_both))
# tagCounts1[df$tcIndex] <- df$tcCount
#
# # Dataframe to be merged
# df <- data.frame(tcIndex=queryHits(fo),
# 								 ctssIndex=subjectHits(fo))
# df$ctssScores=score(gr)[df$ctssIndex]
# df <- aggregate(ctssScores ~ tcIndex, data=df, FUN=sum)
#
# # Allocate zero count vector
# tagCounts2 <- rep(0, length(TCs_both))
# tagCounts2[df$tcIndex] <- df$ctssScores
#
# countTCs <- function(tcs, ctss){
# 	# Find overlaps
# 	fo <- findOverlaps(query=TCs_both, subject=gr)
#
# 	# Dataframe to be merged
# 	df <- data.frame(tcIndex=queryHits(fo),
# 									 ctssIndex=subjectHits(fo))
# 	df$ctssScores=score(ctss)[df$ctssIndex]
# 	df <- aggregate(ctssScores ~ tcIndex, data=df, FUN=sum)
#
# 	# Allocate zero count vector
# 	tagCounts <- rep(0, length(TCs_both))
# 	tagCounts[df$tcIndex] <- df$ctssScores
#
# 	# Return
# 	tagCounts
# }
#
# countTCs(tcs=TCs_both, ctss=gr)
#
# EM <- lapply(ctss_grs, countTCs, tcs=TCs_both) %>% do.call(cbind,.)
#
# # Stats for TCs
# vectors_plus <- coverage_plus[TCs_plus]
# names(vectors_plus) <- NULL
# vectors_minus <- coverage_minus[TCs_minus]
# names(vectors_minus) <- NULL
#
# peaks_plus <- unlist(lapply(vectors_plus, which.max))
# peaks_minus <- unlist(lapply(vectors_minus, which.max))
#
# ranges_plus <- IRanges(start=start(TCs_plus) + peaks_plus - 1, width=1)
# ranges_minus <- IRanges(start=start(TCs_minus) + peaks_minus - 1, width=1)
#
# TCs_plus$peak <- ranges_plus
# TCs_minus$peak <- ranges_minus
#
# # Version with views, much faster
# views_plus <- Views(coverage_plus, reduced_plus)
# peaks_plus <- unlist(viewWhichMaxs(views_plus))
# sum_plus <- unlist(viewSums(views_plus))
#
# views_minus <- Views(coverage_minus, reduced_minus)
# peaks_minus <- unlist(viewWhichMaxs(views_minus))
# sum_minus <- unlist(viewSums(views_minus))
#
# ### Visualize with gviz
#
# # Subset to more managebale size
# tcs <- subset(TCs_both, score >= 100)
#
# # Get annotation
# library(GenomicFeatures)
# library(biomaRt)
#
# txdb<-makeTxDbFromBiomart(biomart ="fungal_mart",
# 													dataset = "spombe_eg_gene",
# 													host="fungi.ensembl.org")
#
# # Tracks
# options(ucscChromosomeNames=FALSE)
# trackinfo <- GenomeAxisTrack()
# txs <- GeneRegionTrack(txdb,geneSymbols=TRUE)
# plustrack <- DataTrack(subset(gr, strand == "+"), name="+", data="tpm")
# minustrack <- DataTrack(subset(gr, strand == "-"), name="-", data="tpm")
# tctrack <- AnnotationTrack(tcs, name='TCs', shape='box')
#
# plotTracks(list(trackinfo, txs, plustrack, minustrack, tctrack),
# 					 chromosome="I", from=1813740, to=1815796, extend.right=100, extend.left=100,
# 					 transcriptAnnotation="symbol",
# 					 type="histogram")
#
#
#
#
