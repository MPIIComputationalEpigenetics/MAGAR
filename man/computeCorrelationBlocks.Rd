% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_blocks.R
\name{computeCorrelationBlocks}
\alias{computeCorrelationBlocks}
\title{computeCorrelationBlocks}
\usage{
computeCorrelationBlocks(
  meth.data,
  annotation,
  cor.threshold = qtlGetOption("cluster.cor.threshold"),
  sd.gauss = qtlGetOption("standard.deviation.gauss"),
  absolute.cutoff = qtlGetOption("absolute.distance.cutoff"),
  max.cpgs = qtlGetOption("max.cpgs"),
  assembly = "hg19",
  chromosome = "chr1"
)
}
\arguments{
\item{meth.data}{A \code{data.frame} containing the methylation data with CpGs in the rows and samples in the columns.}

\item{annotation}{The genomic annotations of the CpG positions.}

\item{cor.threshold}{The correlation threshold used to discard edges from the correlation-based network.}

\item{sd.gauss}{Standard deviation of the Gauss distribution used to weight the distance}

\item{absolute.cutoff}{Absolute distance cutoff after which no methQTL interaction is to be considered.}

\item{max.cpgs}{Maximum number of CpGs used in the computation (used to save memory). 40,000 is a reasonable
default for machines with ~128GB of main memory. Should be smaller for smaller machines and larger
for larger ones.}

\item{assembly}{The assembly used}

\item{chromosome}{The chromosome for which correlation block calling is to be performed}
}
\value{
A list representing the clustering of CpGs into correlation blocks. Each element is a cluster, which contains
   row indices of the DNA methylation matrix that correspond to this cluster.
}
\description{
This function computes CpG correlation blocks from correlations of CpGs across samples by Louvian
clustering.
}
\details{
This method performs clustering of the correlation matrix obtaind from the DNA methylation matrix. Correlations
   are computed for each pair of CpGs across all the samples. We then compute a similarity matrix from this correlation
   matrix and set correlations lower than the given threshold to 0. In the next step, we weight the correlations
   by the distance between the CpGs: smaller distances get higher weights according to Gaussian distribution with
   mean 0 and standard deviation as specified above. Furthermore, similarities of CpGs that are further than
   \code{absolute.distance.cutoff} away from one another are discarded.

   We then compute the associated weighted, undirected graph from the similarity matrix and execute Louvain clustering
   on the graph. The resulting clusters of CpGs are returned.
}
\examples{
meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
meth.data <- getMethData(meth.qtl)
anno.meth <- getAnno(meth.qtl,"meth")
cor.blocks <- computeCorrelationBlocks(meth.data[seq(1,10),],annotation=anno.meth[seq(1,10),])
}
\author{
Michael Scherer
}
