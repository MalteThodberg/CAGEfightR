% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find.R
\name{findLinks}
\alias{findLinks}
\alias{findLinks,GRanges-method}
\alias{findLinks,RangedSummarizedExperiment-method}
\title{Find nearby pairs of clusters and calculate pairwise correlations.}
\usage{
findLinks(object, ...)

\S4method{findLinks}{GRanges}(object, maxDist = 10000L, directional = NULL)

\S4method{findLinks}{RangedSummarizedExperiment}(
  object,
  inputAssay,
  maxDist = 10000L,
  directional = NULL,
  corFun = stats::cor.test,
  vals = c("estimate", "p.value"),
  ...
)
}
\arguments{
\item{object}{GRanges or RangedSummarizedExperiment: Clusters, possibly with
expression for calculating correlations.}

\item{...}{additional arguments passed to methods or ultimately corFun.}

\item{maxDist}{integer: Maximum distance between links.}

\item{directional}{character: Name of a column in object holding a grouping
of the clusters. This must be a factor with two levels. The first level is
used as the basis for calculating orientation (see below).}

\item{inputAssay}{character: Name of assay holding expression values (if
object is a RangedSummarizedExperiment)}

\item{corFun}{function: Function for calculating pairwise correlations. See
notes for supplying custom functions.}

\item{vals}{character: Statistics extracted from the results produced by
corFun. See notes for supplying custom functions.}
}
\value{
A GInteractions holding the links, along with the distance between
  them and correlation estimate and p-value calculated from their expression.
  If a directional analysis was performed, the two anchors are always
  connecting members of the two classes and the orientation of the second
  anchor relative to the first is additionaly calculated (e.g. whether an
  enhancers is upstream or downstream of the TSS).
}
\description{
Finds all links or pairs of clusters within a certain distance of each other
and then calculates the correlation between them. The links found can be
restricted to only be between two classes, for example TSSs to enhancers.
}
\details{
A custom function for calculation correlations can be supplied by
  the user. The output of this function must be a named list or vector of
  numeric values. The names of the vals to be extracted should be supplied to
  vals.
}
\examples{
library(InteractionSet)

# Subset to highly expressed unidirectional clusters
TCs <- subset(exampleUnidirectional, score > 10)

# Find links within a certain distance
findLinks(TCs, inputAssay="counts", maxDist=10000L)

# To find TSS-to-enhancer type links, first merge the clusters:
colData(exampleBidirectional) <- colData(TCs)
rowRanges(TCs)$clusterType <- "TSS"
rowRanges(exampleBidirectional)$clusterType <- "Enhancer"
SE <- combineClusters(TCs, exampleBidirectional, removeIfOverlapping="object1")
rowRanges(SE)$clusterType <- factor(rowRanges(SE)$clusterType, levels=c("TSS", "Enhancer"))

# Calculate kendall correlations of TPM values:
SE <- calcTPM(SE, totalTags="totalTags")
findLinks(SE, inputAssay="TPM", maxDist=10000L, directional="clusterType", method="kendall")
}
\seealso{
Other Spatial functions: 
\code{\link{findStretches}()},
\code{\link{trackLinks}()}
}
\concept{Spatial functions}
