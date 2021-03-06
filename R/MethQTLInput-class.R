##########################################################################################
# MethQTLInput-class.R
# created: 2019-08-27
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# MethQTLInput class definition
##########################################################################################

##########################################################################################
# CLASS DEFINITIONS
##########################################################################################
setClassUnion("matrixOrHDF",c("matrix","HDF5Matrix"))
setClassUnion("characterOrNULL",c("character","NULL"))
setClassUnion("GRangesOrNULL",c("GRanges","NULL"))

#'MethQTLInput-class
#'
#'Class storing methQTL input data, such as DNA methylation and genotyping data, as well as sample metadata
#'
#'@details
#'This class is the basis for computing methQTLs in the methQTL-package. It stores all the relevant information
#'including methylation data and genotype data as a matrix or HDF5Matrix, the phenotypic data as a data frame
#'and the genomic annotation of both the methylation sites and the SNP data.
#'
#'@section Slots:
#'\describe{
#'\item{\code{meth.data}}{The methylation data as a numeric matrix of beta values or as an object of type \code{\link{HDF5Matrix}}}
#'\item{\code{geno.data}}{The genotyping data as a numeric matrix of SNP genotypes (0=homozygote reference,
#'    1=heterozygote, 2=homozygote alternative allele) or as an object of type \code{\link{HDF5Matrix}}}
#'\item{\code{pheno.data}}{Phenotypic data describing the samples used in the study. Matches the dimensions of
#'    both \code{meth.data} and \code{geno.data}}
#'\item{\code{anno.meth}}{Genomic annotation of the methylation sites as a \code{data.frame}. Has the same number
#'    of rows as \code{meth.data}.}
#'\item{\code{anno.geno}}{Genomic annotation of the SNPs as a \code{data.frame}. Has the same number of rows as
#'    \code{geno.data}.}
#'\item{\code{samples}}{The sample identifiers used both for \code{meth.data} and \code{geno.data}, and as the rownames of
#'    \code{pheno.data}.}
#'\item{\code{assembly}}{The genome assembly used.}
#'\item{\code{platform}}{The platform used to compute the methylation data.}
#'\item{\code{disk.dump}}{Flag indicating if the matrices are stored on disk rather than in memory.}
#'\item{\code{imputed}}{Flag indicating if genotype dataset has been imputed.}
#'}
#'@section Methods:
#'\describe{
#'\item{\code{\link[=getMethData,methQTL-method]{getMeth}}}{Returns the methylation matrix.}
#'\item{\code{\link[=getGeno,methQTL-method]{getGeno}}}{Returns the genotyping matrix.}
#'\item{\code{\link[=getPheno,methQTL-method]{getPheno}}}{Returns the phenotypic information.}
#'\item{\code{\link[=getAnno,methQTL-method]{getAnno}}}{Returns the genomic annotation.}
#'\item{\code{\link[=saveMethQTLInput,MethQTLInput-method]{saveMethQTLInput}}}{Stores the object on disk.}
#'\item{\code{\link[=imputeMeth,methQTL-method]{imputeMeth}}}{Imputes the DNA methylation data matrix}
#'}
#'
#'@name MethQTLInput-class
#'@rdname MethQTLInput-class
#'@author    Michael Scherer
#'@exportClass MethQTLInput
setClass("MethQTLInput",
        representation(
            meth.data="matrixOrHDF",
            geno.data="matrixOrHDF",
            pheno.data="data.frame",
            anno.meth="data.frame",
            anno.geno="data.frame",
            samples="characterOrNULL",
            assembly="character",
            disk.dump="logical",
            imputed="logical",
            platform="character"),
        prototype(
            meth.data=matrix(nrow=0,ncol=0),
            geno.data=matrix(nrow=0,ncol=0),
            pheno.data=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            samples=c(),
            assembly="hg19",
            disk.dump=FALSE,
            imputed=FALSE,
            platform="probesEPIC"),
        package="MAGAR")

# CONSTRUCTOR
setMethod("initialize","MethQTLInput",
            function(.Object,
            meth.data=matrix(nrow=0,ncol=0),
            geno.data=matrix(nrow=0,ncol=0),
            pheno.data=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            samples=c(),
            assembly="hg19",
            disk.dump=FALSE,
            imputed=FALSE,
            platform="probesEPIC"){
            if(length(samples) != ncol(meth.data) |
                length(samples) != ncol(geno.data) |
                length(samples) != nrow(pheno.data)){
                stop("Samples do not match dimension of the matrices.")
            }
            .Object@meth.data <- meth.data
            .Object@geno.data <- geno.data
            .Object@pheno.data <- pheno.data
            .Object@anno.meth <- anno.meth
            .Object@anno.geno <- anno.geno
            .Object@samples <- samples
            .Object@assembly <- assembly
            .Object@disk.dump <- disk.dump
            .Object@imputed <- imputed
            .Object@platform <- platform

            .Object
            })

#'get.value
#'
#'This functions returns the values for a particular sample and position from the matrix given in the first argument
#'
#'@param    mat The matrix from which data is to be extracted. Can be either methylation of genotyping data.
#'@param    site The sites to be selected either as a numeric or logical vector. If \code{NULL} all sites are returned.
#'@param    sample The samples to be selected either as a numeric or logical vector. If \code{NULL} all samples are returned.
#'@return    The selected values from the matrix either as a matrix or \code{HDF5Matrix}.
#'@author    Michael Scherer
#'@noRd
get.value <- function(mat,site=NULL,sample=NULL){
    if(!is.element(class(site),c("NULL","integer","numeric","logical"))){
    stop("Invalid value for site, needs to be numeric, logical or NULL")
    }
    if(!is.element(class(sample),c("NULL","integer","numeric","logical","character"))){
    stop("Invalid value for sample, needs to be numeric, character, logical or NULL")
    }
    if(is.null(site)){
    if(is.null(sample)){
        res <- mat
    }else{
        res <- mat[,sample]
    }
    }else{
    if(is.null(sample)){
        res <- mat[site,]
    }else{
        res <- mat[site,sample]
    }
    }
    return(res)
}
##########################################################################################
# GETTERS
##########################################################################################

if(!isGeneric("getMethData")) setGeneric("getMethData",
                                        function(object,...) standardGeneric("getMethData"))

#'getMethData
#'
#'Returns methylation information for the given dataset.
#'
#'@param    object An object of class \code{\link{MethQTLInput-class}}.
#'@param    site The sites to be selected either as a numeric or logical vector. If \code{NULL} all sites are returned.
#'@param    sample The samples to be selected either as a numeric or logical vector. If \code{NULL} all samples are returned.
#'@return    The methylation matrix either as a matrix of \code{\link{HDF5Matrix}}.
#'
#'@rdname getMethData
#'@docType methods
#'@aliases getMethData,methQTL-method
#'@aliases getMethData
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'head(getMethData(meth.qtl))
setMethod("getMethData",signature(object="MethQTLInput"),
            function(object,site=NULL,sample=NULL){
                ret.mat <- get.value(object@meth.data,site=site,sample=sample)
                if(object@disk.dump){
                colnames(ret.mat) <- getSamples(object)
                }
                return(ret.mat)
    }
)

if(!isGeneric("getGeno")) setGeneric("getGeno",
                                    function(object,...) standardGeneric("getGeno"))

#'getGeno
#'
#'Returns genotyping information for the given dataset.
#'
#'@param    object An object of class \code{\link{MethQTLInput-class}}.
#'@param    site The sites to be selected either as a numeric or logical vector. If \code{NULL} all sites are returned.
#'@param    sample The samples to be selected either as a numeric or logical vector. If \code{NULL} all samples are returned.
#'@return    The genotyping matrix either as a matrix of \code{\link{HDF5Matrix}}.
#'
#'@rdname getGeno
#'@docType methods
#'@aliases getGeno,methQTL-method
#'@aliases getGeno
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'head(getGeno(meth.qtl))
setMethod("getGeno",signature(object="MethQTLInput"),
        function(object,site=NULL,sample=NULL){
            ret.mat <- get.value(object@geno.data,site=site,sample=sample)
            if(object@disk.dump){
            colnames(ret.mat) <- getSamples(object)
            }
            return(ret.mat)
        }
)

if(!isGeneric("getPheno")) setGeneric("getPheno",
                                        function(object) standardGeneric("getPheno"))

#'getPheno
#'
#'Returns phenotypic information for the given dataset.
#'
#'@param    object An object of class \code{\link{MethQTLInput-class}}.
#'@return    The phenotypic data either as a \code{data.frame}.
#'
#'@rdname getPheno
#'@docType methods
#'@aliases getPheno,methQTL-method
#'@aliases getPheno
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'head(getPheno(meth.qtl))
setMethod("getPheno",signature(object="MethQTLInput"),
            function(object){
            return(object@pheno.data)
            }
)

if(!isGeneric("getAnno")) setGeneric("getAnno",
                                    function(object,...) standardGeneric("getAnno"))

#'getAnno
#'
#'Returns genomic annotation information for the given dataset.
#'
#'@param    object An object of class \code{\link{MethQTLInput-class}} or \code{\link{MethQTLResult-class}}.
#'@param    type The type of annotation to be returned. Can either be \code{'meth'} or \code{'geno'} for methylation,
#'and genotyping information, respectively.
#'@return    The genomic annotation as a \code{data.frame}.
#'
#'@rdname getAnno
#'@docType methods
#'@aliases getAnno,methQTL-method
#'@aliases getAnno
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'head(getAnno(meth.qtl,"meth"))
#'head(getAnno(meth.qtl,"geno"))
setMethod("getAnno",signature(object="MethQTLInput"),
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

if(!isGeneric("getSamples")) setGeneric("getSamples",
                                        function(object) standardGeneric("getSamples"))

#'getSamples
#'
#'Returns the samples of the given dataset.
#'
#'@param    object An object of class \code{\link{MethQTLInput-class}}.
#'@return    The samples of the dataset as a character vector.
#'
#'@rdname getSamples
#'@docType methods
#'@aliases getSamples,methQTL-method
#'@aliases getSamples
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'getSamples(meth.qtl)
setMethod("getSamples",signature(object="MethQTLInput"),
            function(object){
            return(object@samples)
            }
)

if(!isGeneric("imputeMeth")) setGeneric("imputeMeth",
                                        function(object) standardGeneric("imputeMeth"))

#'imputeMeth
#'
#'Replaces missing values in the DNA methylation data matrix by imputed values
#'
#'@param    object An object of class \code{\link{MethQTLInput-class}}.
#'@return    The object with imputed values.
#'
#'@rdname imputeMeth
#'@import    impute
#'@docType methods
#'@aliases imputeMeth,methQTL-method
#'@aliases imputeMeth
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'meth.qtl.imp <- imputeMeth(meth.qtl)
setMethod("imputeMeth",signature(object="MethQTLInput"),
            function(object){
            rnb.xml2options(qtlGetOption("rnbeads.options"))
            if(!object@disk.dump){
                object@meth.data <- rnb.execute.imputation(object@meth.data)
            }else{
                meth.mat <- as.matrix(object@meth.data)
                meth.mat <- rnb.execute.imputation(meth.mat)
                object@meth.data <- writeHDF5Array(meth.mat)
            }
            return(object)
            }
)

setMethod("show","MethQTLInput",
    function(object){
    ret.str <- list()
    ret.str[1] <- "Object of class MethQTLInput\n"
    ret.str[2] <- paste("\t Contains",length(object@samples),"samples\n")
    ret.str[3] <- paste("\t Methylation data for",nrow(object@meth.data),"CpGs\n")
    ret.str[4] <- paste("\t Genotyping data for",nrow(object@geno.data),"SNPs\n")
    ret.str[5] <- paste("\t Genome assembly:",object@assembly,"\n")
    cat(do.call(paste0,ret.str))
    }
)

if(!isGeneric("saveMethQTLInput")) setGeneric("saveMethQTLInput",
                                        function(object,...)standardGeneric("saveMethQTLInput"))

#'saveMethQTLInput
#'
#'This functions stores a MethQTLInput object in disk.
#'
#'@param    object The \code{\link{MethQTLInput-class}} object to be stored on disk.
#'@param    path A path to a non-existing directory for files to be stored.
#'@return    None
#'
#'@rdname saveMethQTLInput
#'@docType methods
#'@aliases saveMethQTLInput,methQTL-method
#'@aliases saveMethQTLInput
#'@author    Michael Scherer
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'saveMethQTLInput(meth.qtl,"MethQTLInput")
setMethod("saveMethQTLInput","MethQTLInput",
            function(object,path){
            if(file.exists(path)){
                if(dir.exists(path)){
                path <- file.path(path,"methQTL")
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
            if(!object@disk.dump){
                saveRDS(object@meth.data,file=file.path(path,"meth_data.RDS"))
                saveRDS(object@geno.data,file=file.path(path,"geno_data.RDS"))
            }else{
                writeHDF5Array(object@meth.data,filepath = file.path(path,"meth_data.h5"),name="meth.data")
                writeHDF5Array(object@geno.data,filepath = file.path(path,"geno_data.h5"),name="geno.data")
            }
            saveRDS(object@anno.meth,file=file.path(path,"anno_meth.RDS"))
            saveRDS(object@anno.geno,file=file.path(path,"anno_geno.RDS"))
            saveRDS(object@pheno.data,file=file.path(path,"pheno_data.RDS"))
            object@meth.data <- matrix(nrow = 0,ncol = 0)
            object@geno.data <- matrix(nrow = 0,ncol = 0)
            object@anno.meth <- data.frame()
            object@anno.geno <- data.frame()
            object@pheno.data <- data.frame()
            save(object,file=file.path(path,"MethQTLInput.RData"))
            }
)

#'loadMethQTLInput
#'
#'This functions load a \code{\link{MethQTLInput-class}} object from disk.
#'
#'@param    path Path to the directory that has been created by \code{\link{saveMethQTLInput,MethQTLInput-method}}.
#'@return    The object of type \code{\link{MethQTLInput-class}} that has been stored on disk.
#'@author    Michael Scherer
#'@export
#'@examples
#'meth.qtl <- loadMethQTLInput(system.file("extdata","reduced_methQTL",package="MAGAR"))
#'meth.qtl
loadMethQTLInput <- function(path){
    if(any(!(file.exists(file.path(path,"meth_data.RDS"))||file.exists(file.path(path,"meth_data.h5"))),
        !(file.exists(file.path(path,"geno_data.RDS"))||file.exists(file.path(path,"geno_data.h5"))),
        !file.exists(file.path(path,"anno_meth.RDS")),
        !file.exists(file.path(path,"anno_geno.RDS")),
        !file.exists(file.path(path,"pheno_data.RDS")))){
    stop("Invalid value for path. Potentially not a directory saved with saveMethQTLInput")
    }
    load_env<-new.env(parent=emptyenv())
    load(file.path(path, "MethQTLInput.RData"),envir=load_env)
    object <- get("object",load_env)
    is.dumped <- object@disk.dump
    if(!is.dumped){
    meth.data <- readRDS(file.path(path,"meth_data.RDS"))
    geno.data <- readRDS(file.path(path,"geno_data.RDS"))
    }else{
    meth.data <- HDF5Array(filepath=file.path(path,"meth_data.h5"),name="meth.data")
    geno.data <- HDF5Array(filepath=file.path(path,"geno_data.h5"),name="geno.data")
    }
    anno.meth <- readRDS(file.path(path,"anno_meth.RDS"))
    anno.geno <- readRDS(file.path(path,"anno_geno.RDS"))
    pheno.data <- readRDS(file.path(path,"pheno_data.RDS"))
    object@meth.data <- meth.data
    object@geno.data <- geno.data
    object@anno.meth <- anno.meth
    object@anno.geno <- anno.geno
    object@pheno.data <- pheno.data
    return(object)
}
