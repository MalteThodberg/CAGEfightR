#' Example CAGE Data
#'
#' Subset of the CAGE dataset from the paper 'Identification of Gene
#' Transcription Start Sites and Enhancers Responding to Pulmonary Carbon
#' Nanotube Exposure in Vivo'. CTSS data from subsets of chr18 and chr19 across
#' 3 mouse (mm9 ) samples are included. Datasets can be loaded with the data
#' function.
#'
#' @format Example data from various stages of CAGEfightR: \describe{
#'   \item{exampleDesign}{DataFrame: Description of samples, including .bw
#'   filenames} \item{exampleCTSS}{RangedSummarizedExperiment: CTSSs}
#'   \item{exampleUnidirectional}{RangedSummarizedExperiment: Unidirectional or
#'   Tag Clusters}
#'   \item{exampleBidirectionalCluster}{RangedSummarizedExperiment:
#'   Bidirectional clusters} \item{exampleGenes}{RangedSummarizedExperiment:
#'   Genes} }
#' @source \url{http://pubs.acs.org/doi/abs/10.1021/acsnano.6b07533}
#' @examples
#' data(exampleDesign)
#' data(exampleCTSSs)
#' data(exampleUnidirectional)
#' data(exampleBidirectional)
#' data(exampleGenes)
"exampleDesign"

#' @rdname exampleDesign
"exampleCTSSs"

#' @rdname exampleDesign
"exampleUnidirectional"

#' @rdname exampleDesign
"exampleBidirectional"

#' @rdname exampleDesign
"exampleGenes"

