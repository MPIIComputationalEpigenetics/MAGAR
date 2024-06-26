% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpretation.R
\name{qtlAnnotationEnrichment}
\alias{qtlAnnotationEnrichment}
\title{qtlAnnotationEnrichment}
\usage{
qtlAnnotationEnrichment(
  meth.qtl.res,
  type = "SNP",
  annotation = "cpgislands",
  assembly = "hg19"
)
}
\arguments{
\item{meth.qtl.res}{An object of type \code{\link{MethQTLResult-class}} or a list of such objects.}

\item{type}{The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
or \code{'cor.block'}}

\item{annotation}{The genomic annotation to be used. Can be the ones available in \code{\link{rnb.region.types}} or
\code{"ctcf", "distal", "dnase", "proximal", "tfbs", "tss"}}

\item{assembly}{The assembly used. Only \code{"hg19"} (default) and \code{"hg38"} supported.}
}
\value{
A list of two p-values named \code{'enrichment'} for overrepresentation and \code{'depletion'} for underrepresentation
}
\description{
This functions performs enrichment analysis using the Fisher's test for the methQTLs detected
with respect to different genomic annotations.
}
\details{
We use all data points that have been used to calculate methQTLs as the background
and compare the overlaps with the annotation of interest in comparison to the methQTLs that
have been computed in case a \code{\link{MethQTLResult-class}} is provided. If a list of \code{\link{MethQTLResult-class}} objects
is provided, the intersection between the methQTLs from all objects in the list is compared with the union of all interactions
that have been tested.
}
\examples{
meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
res <- qtlAnnotationEnrichment(meth.qtl.res)
}
\author{
Michael Scherer
}
