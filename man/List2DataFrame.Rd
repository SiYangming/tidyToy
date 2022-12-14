% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataSummary.R
\name{List2DataFrame}
\alias{List2DataFrame}
\title{Convert a list to data frame by row or column}
\usage{
List2DataFrame(x, by = "row")
}
\arguments{
\item{by}{row pr col}

\item{pValue}{a vector of p value}
}
\value{
formatted p value
}
\description{
Convert a list to data frame by row or column
}
\examples{
l <- list(rnorm(10), runif(10), rnorm(10, 1))
df <- List2DataFrame(l)
}
\author{
Yangming si
}
