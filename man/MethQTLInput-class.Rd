% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethQTLInput-class.R
\docType{class}
\name{MethQTLInput-class}
\alias{MethQTLInput-class}
\title{MethQTLInput-class}
\description{
Class storing methQTL input data, such as DNA methylation and genotyping data, as well as sample metadata
}
\details{
This class is the basis for computing methQTLs in the methQTL-package. It stores all the relevant information
including methylation data and genotype data as a matrix or HDF5Matrix, the phenotypic data as a data frame
and the genomic annotation of both the methylation sites and the SNP data.
}
\section{Slots}{

\describe{
\item{\code{meth.data}}{The methylation data as a numeric matrix of beta values or as an object of type \code{\link{HDF5Matrix}}}
\item{\code{geno.data}}{The genotyping data as a numeric matrix of SNP genotypes (0=homozygote reference,
   1=heterozygote, 2=homozygote alternative allele) or as an object of type \code{\link{HDF5Matrix}}}
\item{\code{pheno.data}}{Phenotypic data describing the samples used in the study. Matches the dimensions of
   both \code{meth.data} and \code{geno.data}}
\item{\code{anno.meth}}{Genomic annotation of the methylation sites as a \code{data.frame}. Has the same number
   of rows as \code{meth.data}.}
\item{\code{anno.geno}}{Genomic annotation of the SNPs as a \code{data.frame}. Has the same number of rows as
   \code{geno.data}.}
\item{\code{samples}}{The sample identifiers used both for \code{meth.data} and \code{geno.data}, and as the rownames of
   \code{pheno.data}.}
\item{\code{assembly}}{The genome assembly used.}
\item{\code{platform}}{The platform used to compute the methylation data.}
\item{\code{disk.dump}}{Flag indicating if the matrices are stored on disk rather than in memory.}
\item{\code{imputed}}{Flag indicating if genotype dataset has been imputed.}
}
}

\section{Methods}{

\describe{
\item{\code{\link[=getMethData,methQTL-method]{getMeth}}}{Returns the methylation matrix.}
\item{\code{\link[=getGeno,methQTL-method]{getGeno}}}{Returns the genotyping matrix.}
\item{\code{\link[=getPheno,methQTL-method]{getPheno}}}{Returns the phenotypic information.}
\item{\code{\link[=getAnno,methQTL-method]{getAnno}}}{Returns the genomic annotation.}
\item{\code{\link[=saveMethQTLInput,MethQTLInput-method]{saveMethQTLInput}}}{Stores the object on disk.}
\item{\code{\link[=imputeMeth,methQTL-method]{imputeMeth}}}{Imputes the DNA methylation data matrix}
}
}

\author{
Michael Scherer
}
