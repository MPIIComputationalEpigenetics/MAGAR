% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{qtlDistanceScatterplot}
\alias{qtlDistanceScatterplot}
\title{qtlDistanceScatterplot}
\usage{
qtlDistanceScatterplot(meth.qtl.result, out.dir = NULL, out.name = NULL)
}
\arguments{
\item{meth.qtl.result}{An object of type \code{\link{MethQTLResult-class}} containing called methQTL}

\item{out.dir}{If specified, the plot is stored as a pdf in this directory}

\item{out.name}{Optional name for the resulting plot}
}
\value{
An object of type \code{ggplot} comparing the distance between CpG and SNP. Negative values indicate that the
       SNP is downstream of the CpG.
}
\description{
Computes a scatterplot between CpG-SNP distance with both effect size and p-value
}
\examples{
meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
qtlDistanceScatterplot(meth.qtl.res)
}
\author{
Michael Scherer
}
