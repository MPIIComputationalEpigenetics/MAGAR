##########################################################################################
# compute_methQTL.R
# created: 2020-03-26
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods for calling correlation blocks from DNA methylation data.
##########################################################################################

#'computeCorrelationBlocks
#'
#'This function computes CpG correlation blocks from correlations of CpGs across samples by Louvian
#'clustering.
#'
#'@param    meth.data A \code{data.frame} containing the methylation data with CpGs in the rows and samples in the columns.
#'@param    annotation The genomic annotations of the CpG positions.
#'@param    cor.threshold The correlation threshold used to discard edges from the correlation-based network.
#'@param    sd.gauss Standard deviation of the Gauss distribution used to weight the distance
#'@param    absolute.cutoff Absolute distance cutoff after which no methQTL interaction is to be considered.
#'@param    max.cpgs Maximum number of CpGs used in the computation (used to save memory). 40,000 is a reasonable
#'        default for machines with ~128GB of main memory. Should be smaller for smaller machines and larger
#'        for larger ones.
#'@param    assembly The assembly used
#'@param    chromosome The chromosome for which correlation block calling is to be performed
#'@return    A list representing the clustering of CpGs into correlation blocks. Each element is a cluster, which contains
#'    row indices of the DNA methylation matrix that correspond to this cluster.
#'@details    This method performs clustering of the correlation matrix obtaind from the DNA methylation matrix. Correlations
#'    are computed for each pair of CpGs across all the samples. We then compute a similarity matrix from this correlation
#'    matrix and set correlations lower than the given threshold to 0. In the next step, we weight the correlations
#'    by the distance between the CpGs: smaller distances get higher weights according to Gaussian distribution with
#'    mean 0 and standard deviation as specified above. Furthermore, similarities of CpGs that are further than
#'    \code{absolute.distance.cutoff} away from one another are discarded.
#'
#'    We then compute the associated weighted, undirected graph from the similarity matrix and execute Louvain clustering
#'    on the graph. The resulting clusters of CpGs are returned.
#'
#'@author    Michael Scherer
#'@export
#'@importFrom igraph graph.adjacency cluster_louvain groups
#'@importFrom stats cor dnorm
#'@import    bigstatsr
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'meth.data <- getMethData(meth.qtl)
#'anno.meth <- getAnno(meth.qtl,"meth")
#'cor.blocks <- computeCorrelationBlocks(meth.data[seq(1,10),],annotation=anno.meth[seq(1,10),])
computeCorrelationBlocks <- function(meth.data,
                                        annotation,
                                        cor.threshold=qtlGetOption("cluster.cor.threshold"),
                                        sd.gauss=qtlGetOption("standard.deviation.gauss"),
                                        absolute.cutoff=qtlGetOption("absolute.distance.cutoff"),
                                        max.cpgs=qtlGetOption("max.cpgs"),
                                        assembly="hg19",
                                        chromosome="chr1"){
    logger.start("Compute correlation blocks")
    if(nrow(annotation)>max.cpgs){
    logger.info(paste("Split workload, since facing",nrow(annotation),"CpGs (Maximum is",max.cpgs,")"))
    bin.split <- round(nrow(annotation)/2)
    return(c(computeCorrelationBlocks(meth.data=meth.data[seq(1,bin.split),],
                                        annotation=annotation[seq(1,bin.split),],
                                        cor.threshold = cor.threshold,
                                        sd.gauss = sd.gauss,
                                        absolute.cutoff = absolute.cutoff,
                                        max.cpgs = max.cpgs,
                                        assembly=assembly,
                                        chromosome=chromosome),
            lapply(computeCorrelationBlocks(meth.data=meth.data[(bin.split+1):nrow(annotation),],
                                        annotation=annotation[(bin.split+1):nrow(annotation),],
                                        cor.threshold = cor.threshold,
                                        sd.gauss = sd.gauss,
                                        absolute.cutoff = absolute.cutoff,
                                        max.cpgs = max.cpgs,
                                        assembly=assembly,
                                        chromosome=chromosome),
                    function(x) x+bin.split
                    )
            ))
    }
    logger.start("Compute correlation matrix")
    if(qtlGetOption("correlation.type")=="pearson"){
        cor.all <- big_cor(as_FBM(t(as.matrix(meth.data)),type="double"))
    }else{
        cor.all <- cor(t(as.matrix(meth.data)),qtlGetOption("correlation.type"))
    }
    rm(meth.data)
    logger.completed()
    cor.all <- cor.all[,,drop= FALSE]
    if(qtlGetOption("hdf5dump")){
        cor.all <- writeHDF5Array(cor.all)
    }
    rep.vals <- cor.all<cor.threshold
    if(qtlGetOption("hdf5dump")){
        rep.vals <- writeHDF5Array(rep.vals)
    }
    cor.all[rep.vals] <- 0
    genomic.positions <- annotation$Start
    logger.start("Compute pairwise distances")
    gc()
    pairwise.distance <- abs(as.data.frame(lapply(genomic.positions,function(x)x-genomic.positions)))
    logger.completed()
    rep.vals <- pairwise.distance>absolute.cutoff
    if(qtlGetOption("hdf5dump")){
        rep.vals <- writeHDF5Array(rep.vals)
    }
    cor.all[rep.vals] <- 0
    gc()
    logger.start("Weight distances")
    if(qtlGetOption("hdf5dump")){
    weighted.distances <- matrix(nrow=nrow(cor.all),ncol=ncol(cor.all))
    weighted.distances <- writeHDF5Array(weighted.distances)
    chunk.size <- 10000
    i <- 1
    while(i < nrow(cor.all)){
        if((i + chunk.size)>nrow(cor.all)){
        do.work <- i:nrow(cor.all)
        weighted.distances[do.work,] <- as.matrix(cor.all[do.work,])*
            dnorm(as.matrix(pairwise.distance[do.work,]),0,sd.gauss)
        break
        }
        do.work <- i:(i+chunk.size)
        weighted.distances[do.work,] <- as.matrix(cor.all[do.work,])*
        dnorm(as.matrix(pairwise.distance[do.work,]),0,sd.gauss)
        i <- i+chunk.size+1
    }
    }else{
        weighted.distances <- cor.all*dnorm(as.matrix(pairwise.distance),0,sd.gauss)
    }
    logger.completed()
    colnames(weighted.distances) <- as.character(seq(1,ncol(weighted.distances)))
    rownames(weighted.distances) <- as.character(seq(1,nrow(weighted.distances)))
    rm(rep.vals)
    rm(cor.all)
    gc()
    logger.start("Compute graph")
    graph.ad <- graph.adjacency(as.matrix(weighted.distances),weighted = TRUE,mode = "undirected",diag= FALSE)
    logger.completed()
    logger.start("Compute clustering")
    clust <- cluster_louvain(graph.ad)
    rm(weighted.distances)
    gc()
    logger.completed()
    logger.completed()
    return(lapply(groups(clust),function(x)as.numeric(x)))
}
