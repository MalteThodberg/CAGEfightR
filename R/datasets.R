#' Example CAGE Data
#'
#' Subset of the CAGE dataset from the paper "Identification of Gene Transcription Start Sites and Enhancers Responding to Pulmonary Carbon Nanotube Exposure in Vivo".
#' CTSS data from chr18 and chr19 across 11 mouse (mm9 )samples are included.
#'
#' @format Example data from various stages of CAGEfightR:
#' \describe{
#'   \item{exampleDesign}{DataFrame: Description of samples, including .bw filenames}
#'   \item{exampleCTSS}{GRangesList: CTSS data for each sample}
#'   \item{exampleCoverage}{GRanges: Coverage in TPM across all samples}
#'   \item{exampleTagCluster}{GRanges: Subset of Tag Clusters}
#'   \item{exampleBidirectionalCluster}{GRanges: Subset of bidirectional clusters}
#' }
#' @source \url{http://pubs.acs.org/doi/abs/10.1021/acsnano.6b07533}
"exampleDesign"

#' @rdname exampleDesign
"exampleCTSS"

#' @rdname exampleDesign
"exampleCoverage"

#' @rdname exampleDesign
"exampleTagClusters"

#' @rdname exampleDesign
"exampleBidirectionalClusters"

