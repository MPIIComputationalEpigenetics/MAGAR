##########################################################################################
# MethQTLResult-class.R
# created: 2019-08-27
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# MethQTLResult class definition
##########################################################################################

#'MethQTLResult-class
#'
#'Class storing methQTL analysis results and the associated genomic annotations
#'
#'@details
#'This class stores the results of the methQTL analysis. It stores a \code{data.frame} with the methQTL results,
#'and associated genomic annotations for both the methylation sites and SNPs.
#'
#'@section Slots:
#'\describe{
#'\item{\code{result.frame}}{The methQTL results as a \code{data.frame}}
#'\item{\code{anno.meth}}{Genomic annotation of the methylation sites as a \code{data.frame}.}
#'\item{\code{anno.geno}}{Genomic annotation of the SNPs as a \code{data.frame}.}
#'\item{\code{correlation.blocks}}{Correlation blocks determined from the methylation matrix.}
#'\item{\code{method}}{The method used to call methQTL.}
#'\item{\code{rep.type}}{Method used to determine representative CpGs from correlation blocks.}
#'\item{\code{chr}}{Optional argument specifying if methQTL were called on a single chromosome.}
#'}
#'@section Methods:
#'\describe{
#'\item{\code{\link[=getResult,MethQTLResult-method]{getResult}}}{Returns the methQTL results.}
#'\item{\code{\link[=getAnno,MethQTLResult-method]{getAnno}}}{Returns the genomic annotation.}
#'}
#'
#'@name MethQTLResult-class
#'@rdname MethQTLResult-class
#'@author    Michael Scherer
#'@exportClass MethQTLResult

setClass("MethQTLResult",
        representation(
            result.frame="data.frame",
            anno.meth="data.frame",
            anno.geno="data.frame",
            correlation.blocks="list",
            method="character",
            rep.type="character",
            chr="characterOrNULL"
        ),
        prototype(
            result.frame=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            correlation.blocks=list(),
            method="classical.linear",
            rep.type="row.medians",
            chr=NULL
        ),
        package="MAGAR")

# CONSTRUCTOR
setMethod("initialize","MethQTLResult",
            function(.Object,
            result.frame=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            correlation.blocks=list(),
            method="classical.linear",
            rep.type="row.medians",
            chr=NULL
            ){
            .Object@result.frame <- result.frame
            .Object@anno.meth <- anno.meth
            .Object@anno.geno <- anno.geno
            .Object@correlation.blocks=correlation.blocks
            .Object@method <- method
            .Object@rep.type <- rep.type
            .Object@chr <- chr

            .Object
            })

##########################################################################################
# GETTERS
##########################################################################################
if(!isGeneric("getResult")) setGeneric("getResult",function(object,...) standardGeneric("getResult"))

#'getResult
#'
#'Returns the methQTL results stores in the object.
#'
#'@param    object An of type \code{\link{MethQTLResult-class}}.
#'@param    cor.blocks Correlation blocks as obtained using \code{getCorrelationBlocks}. Please note that the
#'        correlation blocks need to contain the CpG identifiers, so the \code{\link{MethQTLInput-class}} object
#'        needs to be provided to \code{getCorrelationBlocks}.
#'@param    na.rm Flag indicating if rows containing NA values are to be removed from the result.
#'@return    The methQTL results as a \code{data.frame} with each row being a methQTL.
#'@rdname getResult
#'@docType methods
#'@aliases getResult,MethQTLResult-method
#'@aliases getResult
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'head(getResult(meth.qtl.res))
setMethod("getResult",signature(object="MethQTLResult"),
            function(object,cor.blocks=NULL,na.rm=FALSE){
            ret <- object@result.frame
        if(na.rm){
            keep.lines <- apply(ret,1,function(line){
                any(!is.na(line))
            })
            ret <- ret[keep.lines,]
        }
            if(!is.null(cor.blocks)){
                if(is.list(cor.blocks[[1]])){
                cor.blocks.assigned <- list()
                for(i in seq(1,length(cor.blocks))){
                    chr <- cor.blocks[[i]]
                    for(j in seq(1,length(chr))){
                    block <- chr[[j]]
                    cpg <- unique(as.character(ret$CpG)[as.character(ret$CpG) %in% block])
                    if(length(cpg)>0){
                        cor.blocks.assigned[[unique(cpg)]] <- block
                    }
                    }
                }
                }else{
                cor.blocks.assigned <- list()
                for(i in seq(1,length(cor.blocks))){
                    block <- cor.blocks[[i]]
                    cpg <- unique(as.character(ret$CpG)[as.character(ret$CpG) %in% block])
                    if(length(cpg)>0){
                    cor.blocks.assigned[[cpg]] <- block
                    }
                }
                }
                ret$CorrelationBlock <- cor.blocks.assigned[as.character(ret$CpG)]
            }
            return(ret)
            }
)

if(!isGeneric("getResultGWASMap")) setGeneric("getResultGWASMap",
                                                function(object,...) standardGeneric("getResultGWASMap"))

#'getResultGWASMap
#'
#'Returns the methQTL results in the format used as input to GWAS-map and stores in the object.
#'
#'@param    object An of type \code{\link{MethQTLResult-class}}.
#'@param    meth.qtl An object of type \code{\link{MethQTLInput-class}} containing further information about the QTLs
#'@return    The methQTL results as a \code{data.frame} with each row being a methQTL.
#'@rdname getResultGWASMap
#'@docType methods
#'@aliases getResultGWASMap,MethQTLResult-method
#'@aliases getResultGWASMap
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'head(getResultGWASMap(meth.qtl.res,meth.qtl))
setMethod("getResultGWASMap",signature(object="MethQTLResult"),
            function(object,meth.qtl){
            ret <- object@result.frame
#            keep.lines <- apply(ret,1,function(line){
#                any(!is.na(line))
#            })
#            ret <- ret[keep.lines,]
        anno.geno <- getAnno(meth.qtl,"geno")[as.character(ret$SNP),]
        ret$ReferenceAllele <- anno.geno$Allele.1
        ret$EffectiveAllele <- anno.geno$Allele.2
            ret$EffectiveAlleleFrequency <- anno.geno$Allele.2.Freq
            return(ret)
            }
)

if(!isGeneric("getAnno")) setGeneric("getAnno",
                                    function(object,...) standardGeneric("getAnno"))

#'getAnno
#'@rdname getAnno
#'@docType methods
#'@aliases getAnno,methQTL-method
#'@aliases getAnno
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'head(getAnno(meth.qtl.res,"meth"))
#'head(getAnno(meth.qtl.res,"geno"))
setMethod("getAnno",signature(object="MethQTLResult"),
            function(object,type="meth"){
            if(type=="meth"){
                return(object@anno.meth)
            }else if(type=="geno"){
                return(object@anno.geno)
            }else{
                stop("Invalid value for type: needs to be 'meth' or 'geno'")
            }
            }
)

if(!isGeneric("getCorrelationBlocks")) setGeneric("getCorrelationBlocks",
                                                    function(object) standardGeneric("getCorrelationBlocks"))

#'getCorrelationBlocks
#'
#'Returns the correlation blocks defined for the given dataset
#'
#'@param    object An object of class \code{\link{MethQTLResult-class}}.
#'@return    A \code{list} object containing the correlation blocks.
#'@rdname getCorrelationBlocks
#'@docType methods
#'@aliases getCorrelationBlocks,MethQTLResult-method
#'@aliases getCorrelationBlocks
#'@return
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'head(getCorrelationBlocks(meth.qtl.res))
setMethod("getCorrelationBlocks",signature(object="MethQTLResult"),
            function(object){
            cor.blocks <- object@correlation.blocks
            if(is.list(cor.blocks[[1]])){
                anno <- getAnno(object)
                ret <- list()
                if(grepl("chr",names(cor.blocks))){
                chr.ids <- names(cor.blocks)
                }else{
                chr.ids <- paste0("chr:",seq(1,length(cor.blocks)))
                }
                for(chr.id in chr.ids){
                chr <- cor.blocks[[chr.id]]
                chr.id <- ifelse(chr.id=="chr23","X",ifelse(chr.id=="chr24","Y",chr.id))
                anno.chr <- anno[anno$Chromosome%in%chr.id,]
                res <- lapply(chr,function(cg){
                    row.names(anno.chr)[cg]
                })
                ret[[chr.id]] <- res
                }
                ret <- ret[order(as.numeric(gsub("chr","",names(ret))))]
            }else{
                anno.chr <- getAnno(object)
                anno.chr <- anno.chr[anno.chr$Chromosome%in%object@chr,]
                ret <- lapply(cor.blocks,function(cg){
                row.names(anno.chr)[cg]
                })
            }
            return(ret)
            }
)

setMethod("show","MethQTLResult",
            function(object){
            ret.str <- list()
            ret.str[1] <- "Object of class MethQTLResult\n"
            ret.str[2] <- paste("\t Contains",nrow(object@result.frame),"methQTL\n")
            if(length(object@correlation.blocks)>0){
                if(is.list(object@correlation.blocks[[1]])){
                ret.str[3] <- paste("\t Contains",sum(lengths(object@correlation.blocks)),"correlation blocks\n")
                }else{
                ret.str[3] <- paste("\t Contains",length(object@correlation.blocks),"correlation blocks\n")
                }
            }else{
                ret.str[3] <- "\t Contains 0 correlation blocks\n"
            }
            ret.str[4] <- paste("\t methQTL called using",object@method,"\n")
            ret.str[5] <- paste("\t representative CpGs computed with",object@rep.type,"\n")
            if(!is.null(object@chr)){
                ret.str[6] <- paste("\t methQTL called for chromosome",object@chr,"\n")
            }
            cat(do.call(paste0,ret.str))
            }
)

if(!isGeneric("filterPval")) setGeneric("filterPval",
                                        function(object,...)standardGeneric("filterPval"))

#'filterPval
#'
#'This functions filters the methQTL results according to a given p-value cutoff
#'
#'@param    object The \code{\link{MethQTLResult-class}} object to be filtered
#'@param    p.val.cutoff The p-value cutoff to be employed
#'@return    The filtered \code{\link{MethQTLResult-class}} object
#'@rdname filterPval
#'@docType methods
#'@aliases filterPval,MethQTLResult-method
#'@aliases filterPval
#'@author    Michael Scherer
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'meth.qtl.res <- filterPval(meth.qtl.res)
#'meth.qtl.res
setMethod("filterPval","MethQTLResult",
            function(object,p.val.cutoff=0.01){
            res <- object@result.frame
            res <- res[res$p.val.adj.fdr <= p.val.cutoff,]
            object@result.frame <- res
            return(object)
            }
)

if(!isGeneric("saveMethQTLResult")) setGeneric("saveMethQTLResult",
    function(object,...)standardGeneric("saveMethQTLResult"))

#'saveMethQTLResult
#'
#'This functions stores a MethQTLInput object in disk.
#'
#'@param    object The \code{\link{MethQTLResult-class}} object to be stored on disk.
#'@param    path A path to a non-existing directory for files to be stored.
#'@return    None
#'
#'@rdname saveMethQTLResult
#'@docType methods
#'@aliases saveMethQTLResult,methQTL-method
#'@aliases saveMethQTLResult
#'@author    Michael Scherer
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'saveMethQTLResult(meth.qtl.res,"MethQTLResult")
setMethod("saveMethQTLResult","MethQTLResult",
            function(object,path){
            if(file.exists(path)){
                if(dir.exists(path)){
                path <- file.path(path,"MethQTLResult")
                if(file.exists(path)){
                    stop("Will not overwrite existing data")
                }
                dir.create(path)
                }else{
                stop("Will not overwrite existing data")
                }
            }else{
                dir.create(path)
            }
            saveRDS(object@result.frame,file=file.path(path,"result_frame.RDS"))
            saveRDS(object@anno.meth,file=file.path(path,"anno_meth.RDS"))
            saveRDS(object@anno.geno,file=file.path(path,"anno_geno.RDS"))
            saveRDS(object@correlation.blocks,file=file.path(path,"correlation_blocks.RDS"))
            object@result.frame <- data.frame()
            object@anno.meth <- data.frame()
            object@anno.geno <- data.frame()
            object@correlation.blocks <- list()
            save(object,file=file.path(path,"MethQTLResult.RData"))
            }
)

#'loadMethQTLResult
#'
#'This functions load a \code{\link{MethQTLResult-class}} object from disk.
#'
#'@param    path Path to the directory that has been created by \code{saveMethQTLResult,MethQTLResult-method}.
#'@return    The object of type \code{\link{MethQTLResult-class}} that has been stored on disk.
#'@author    Michael Scherer
#'@export
#'@examples
#'meth.qtl.res <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'meth.qtl.res
loadMethQTLResult <- function(path){
    if(any(!(file.exists(file.path(path,"result_frame.RDS"))),
        !file.exists(file.path(path,"anno_meth.RDS")),
        !file.exists(file.path(path,"anno_geno.RDS")),
        !file.exists(file.path(path,"correlation_blocks.RDS")))){
    stop("Invalid value for path. Potentially not a directory saved with saveMethQTLResult")
    }
    load_env<-new.env(parent=emptyenv())
    load(file.path(path, "MethQTLResult.RData"),envir=load_env)
    object <- get("object",load_env)
    result.frame <- readRDS(file.path(path,"result_frame.RDS"))
    anno.meth <- readRDS(file.path(path,"anno_meth.RDS"))
    anno.geno <- readRDS(file.path(path,"anno_geno.RDS"))
    correlation.blocks <- readRDS(file.path(path,"correlation_blocks.RDS"))
    object@result.frame <- result.frame
    object@anno.meth <- anno.meth
    object@anno.geno <- anno.geno
    object@correlation.blocks <- correlation.blocks
    return(object)
}

#'joinMethQTLResult
#'
#'This function combines a list of \code{\link{MethQTLResult-class}} objects.
#'
#'@param    obj.list A list of \code{\link{MethQTLResult-class}} objects to be joined
#'@return    An object of type \code{\link{MethQTLResult-class}} containing the combined information
#'@author    Michael Scherer
#'@export
#'@examples
#'meth.qtl.res.1 <- loadMethQTLResult(system.file("extdata","MethQTLResult_chr18",package="MAGAR"))
#'meth.qtl.res.2 <- meth.qtl.res.1
#'meth.qtl.res <- joinMethQTLResult(list(meth.qtl.res.1,meth.qtl.res.2))
joinMethQTLResult <- function(obj.list){
    if(any(!unlist(lapply(obj.list,function(x)inherits(x,"MethQTLResult"))))){
    logger.error("Objects needs to be of type MethQTLResult")
    }
    result.frame <- c()
    anno.meth <- c()
    anno.geno <- c()
    correlation.blocks <- list()
    methods <- c()
    rep.types <- c()
    for(obj in obj.list){
    result.frame <- rbind(result.frame,getResult(obj))
    anno.meth <- rbind(anno.meth,getAnno(obj))
    anno.geno <- rbind(anno.geno,getAnno(obj,"geno"))
    if(!is.null(obj@chr)){
        correlation.blocks[[obj@chr]] <- obj@correlation.blocks
    }else{
        correlation.blocks <- c(correlation.blocks,obj@correlation.blocks)
    }
    methods <- c(methods,obj@method)
    if(any(methods != obj@method)){
        logger.error("Incompatible methQTL calling methods")
    }
    rep.types <- c(rep.types,obj@rep.type)
    if(any(rep.types != obj@rep.type)){
        logger.error("Incompatible representative CpG computation methods")
    }
    }
    if(!is.null(result.frame)&&nrow(result.frame)>0){
    result.frame <- data.frame(result.frame[order(result.frame[,1]),])
    anno.meth <- data.frame(anno.meth[order(anno.meth[,1]),])
    anno.geno <- data.frame(anno.geno[order(anno.geno[,1]),])
    }
    ret.obj <- new("MethQTLResult",
                result.frame=result.frame,
                anno.meth=anno.meth,
                anno.geno=anno.geno,
                correlation.blocks=correlation.blocks,
                method=methods[1],
                rep.type=rep.types[1],
                chr=NULL)
    return(ret.obj)
}
