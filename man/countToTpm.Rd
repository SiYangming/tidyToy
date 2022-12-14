% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAexpression.R
\name{countToTpm}
\alias{countToTpm}
\title{Convert count matrix to tpm expression}
\usage{
countToTpm(counts, effLen)
}
\arguments{
\item{counts}{A count matrix from NGS}

\item{effLen}{effective length from NGS}
}
\value{
A TPM expression matrix
}
\description{
Convert count matrix to tpm expression
}
\examples{
cnts <- c(4250, 3300, 200, 1750, 50, 0)
lens <- c(900, 1020, 2000, 770, 3000, 1777)
countDf <- data.frame(count = cnts, length = lens)
# assume a mean(FLD) = 203.7
countDf$effLength <- countDf$length - 203.7 + 1
countDf$tpm <- with(countDf, countToTpm(count, effLength))
countDf$tpm
}
\author{
Yangming si
}
