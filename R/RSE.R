# #' Collect CAGE data into SummarizedExperiment
# #'
# #' Collects CTSS-GRanges, Tag Cluster (TCs), Expression Matrix (EM) and study design into a SummarizedExperiment. Additional options for documenting the object are available.
# #'
# #' @param ctss GRangesList: CTSS-GRanges.
# #' @param tcs GRanges: TCs
# #' @param em matrix, integer: EM
# #' @param design data.frame: Study note
# #'
# #' @note This function does not sort or match the input, so arguments must be provided with matching columns, rows, etc.
# #'
# #' @return RangedSummarizedExperiment
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @export
# summariseCAGE <- function(ctss, design, ctssCutoff, clusteringFun, ...){
# 	### TO DO
# 	# Option for keepnig ctss in final object
# 	# Get column names from design rownames
# 	# Option for time stamping object
# 	# Option for author stamp
#
# 	tmp <- dim(em)
# 	message("Assembling CAGE experiment of ", tmp[1], " TCs in ", tmp[2], " samples...")
#
# 	# Preserve CTSS files
# 	sample_data <- S4Vectors::DataFrame(ctss=ctss, design)
# 	rownames(sample_data) <- rownames(design)
# 	colnames(em) <- rownames(design)
#
# 	# Assemble
# 	RSE <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=em),
# 															rowRanges=tcs,
# 															colData=sample_data)
#
# 	# Return
# 	RSE
# }
