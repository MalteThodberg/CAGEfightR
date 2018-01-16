#' Example CAGE Data
#'
#' Subset of the CAGE dataset from the paper "Identification of Gene Transcription Start Sites and Enhancers Responding to Pulmonary Carbon Nanotube Exposure in Vivo".
#' CTSS data from subsets of chr18 and chr19 across 11 mouse (mm9 )samples are included. Datasets can be loaded with the data function.
#'
#' @format Example data from various stages of CAGEfightR:
#' \describe{
#'   \item{exampleDesign}{DataFrame: Description of samples, including .bw filenames}
#'   \item{exampleCTSS}{RangedSummarizedExperiment: CTSSs}
#'   \item{exampleUnidirectional}{RangedSummarizedExperiment: Tag clusters (TSSs)}
#'   \item{exampleBidirectionalCluster}{RangedSummarizedExperiment: Bidirectional clusters (enhancers)}
#'   \item{exampleGenes}{RangedSummarizedExperiment: Genes}
#' }
#' @source \url{http://pubs.acs.org/doi/abs/10.1021/acsnano.6b07533}
"exampleDesign"

#' @rdname exampleDesign
"exampleCTSSs"

#' @rdname exampleDesign
"exampleUnidirectional"

#' @rdname exampleDesign
"exampleBidirectional"

#' @rdname exampleDesign
"exampleGenes"

