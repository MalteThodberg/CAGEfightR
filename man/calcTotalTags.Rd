% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.R
\name{calcTotalTags}
\alias{calcTotalTags}
\title{Calculate the total number of CAGE tags across samples.}
\usage{
calcTotalTags(object, inputAssay = "counts", outputColumn = "totalTags")
}
\arguments{
\item{object}{RangedSummarizedExperiment: CAGE data quantified at CTSS,
cluster or gene-level.}

\item{inputAssay}{character: Name of assay holding input expression values.}

\item{outputColumn}{character: Name of column in colData to hold number of
total tags.}
}
\value{
object with total tags per library added as a column in
  colData.
}
\description{
For each CAGE library, calculate the total number of tags.
}
\examples{
data(exampleUnidirectional)
calcTotalTags(exampleUnidirectional)
}
\seealso{
Other Calculation functions: 
\code{\link{calcBidirectionality}()},
\code{\link{calcComposition}()},
\code{\link{calcPooled}()},
\code{\link{calcShape}()},
\code{\link{calcSupport}()},
\code{\link{calcTPM}()},
\code{\link{subsetByBidirectionality}()},
\code{\link{subsetByComposition}()},
\code{\link{subsetBySupport}()}
}
\concept{Calculation functions}
