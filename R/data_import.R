##########################################################################################
# data_import.R
# created: 2019-08-28
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods to perform data import for methQTL analysis.
##########################################################################################

#' doImport
#'
#' Performs input for the given DNA methylation and genotyping data.
#'
#' @param data.location Named character vector specifying the data location. The names correspond to:
#'          \describe{
#'            \item{idat.dir}{Path to the idat-folder for raw DNA methylation data.}
#'            \item{geno.dir}{Path to the genotyping data file directory, either containing PLINK files
#'             (i.e. with \code{bed}, \code{bim} and \code{fam} files), or imputed, dosage files (i.e. with
#'             \code{dos} and \code{txt} files)}
#'          }
#' @param s.anno Path to the sample annotation sheet. If \code{NULL}, the program searches for potential sample
#'          annotation sheets in the data location directories.
#' @param assembly.meth Assembly used for the DNA methylation data. Typically is \code{"hg19"} for Illumina BeadArray
#'          data.
#' @param assembly.geno Assembly used for the genotyping data. If it is not the same as \code{assembly.geno}, the
#'          positions will be matched using liftOver.
#' @param tab.sep The table separator used for the sample annotation sheet.
#' @param s.id.col The column name of the sample annotation sheet that specifies the sample identifier.
#' @param out.folder The output directory to store diagnostic plots
#' @param ... Futher parameters passed to e.g. \code{\link{doGenoImport}}
#' @return An object of type \code{\link{methQTLInput-class}} with the methylation and genotyping information added.
#' @details Import of DNA methylation and genotyping data is done separately:
#'          \describe{
#'            \item{DNA methylation data}{DNA methylation data is imported using the \code{\link{RnBeads}} package. We
#'            use a default option setting commonly used for DNA methylation data obtained from the Illumina BeadArray
#'            series. If you want to specify further options, we refer to the \code{\link{rnb.options}}.}
#'            \item{Genotyping data}{Genotyping data is processed using PLINK. We focus on genotyping data generated
#'            with the Illumina BeadArray series and use default options. For further option settings, you should consult
#'            the \code{\link{qtlSetOption}} documentation.}
#'          }
#'
#'          If \code{assembly.meth} and \code{assembly.geno} do not match, we use the liftOver function to match the
#'          older assembly to the newer one and discard those sites, that are only covered by the newer version.
#'
#' @author Michael Scherer
#' @export
#' @import HDF5Array

doImport <- function(data.location,
                      s.anno=NULL,
                      assembly.meth="hg19",
                      assembly.geno="hg19",
                      tab.sep=",",
                      s.id.col="sample_id",
                      out.folder=tempdir(),
                      ...){

  if(!dir.exists(out.folder)){
    stop("Output directory does not exist.")
  }
  if(dir.exists(file.path(out.folder,"rnbeads_preprocessing"))){
    stop("RnBeads directory (rnbeads_preprocessing) does already exist. Please remove it and restart.")
  }
  logger.start("Import methQTL data")
  if(is.null(s.anno)){
    for(i in 1:length(data.location)){
      all.files <- list.files(data.location,full.names=T)
      s.anno <- all.files[grepl("*annotation.csv",all.files)]
      if(is.null(s.anno)){
        s.anno <- all.files[grepl("*annotation.tsv",all.files)]
        tab.sep <- "\t"
        if(is.null(s.anno) && i==length(data.location)){
          stop("Could not find sample annotation sheet")
        }else if(!is.null(anno)){
          break
        }
      }else{
        break
      }
    }
  }
  pheno.data <- read.table(s.anno,sep=tab.sep,header = T)
  if(ncol(pheno.data)==1){
    logger.warning("Pheno data does only contain one column. Did you specify the wrong 'sep.tab'?")
  }
  geno.import <- doGenoImport(data.location,pheno.data,s.id.col,out.folder,...)
  pheno.data <- geno.import$pheno.data
  s.anno <- file.path(out.folder,ifelse(tab.sep==",","sample_annotation.csv","sample_annotation.tsv"))
  write.table(pheno.data,s.anno,sep=tab.sep)
  meth.import <- doMethImport(data.location,
                                assembly.meth,
                                s.anno,
                                s.id.col,
                                tab.sep,
                                snp.location=geno.import$annotation,
                                out.folder=out.folder)
  s.names <- intersect(meth.import$samples,geno.import$samples)
  sel.meth <- match(s.names,meth.import$samples)
  sel.geno <- match(s.names,geno.import$samples)
  if(any(is.na(sel.geno))||any(is.na(sel.meth))){
	  stop("Sample IDs for methylation and genotyping data do not match")
  }
  meth.data <- meth.import$data[,sel.meth]
  geno.data <- geno.import$data[,sel.geno]
  if(qtlGetOption("hdf5dump")){
  	meth.data <- as(meth.data,"HDF5Matrix")
	  geno.data <- as(geno.data,"HDF5Matrix")
  }
  pheno.dat.names <- as.character(pheno.data[,s.id.col])
  if(length(unique(pheno.dat.names))!=length(pheno.dat.names)){
    stop("Sample identifiers in the sample annotation sheet are not unique!")
  }
  pheno.data <- pheno.data[pheno.dat.names%in%s.names,]
  if(is.null(s.names) || (length(unique(s.names)) < length(s.names))){
    stop("Invalid value for s.id.col, needs to specify unique identfiers in the sample annotation sheet")
  }
#  if(any(!(s.names%in%meth.import$samples) || any!(s.names%in%geno.import$samples)){
#    stop("Samples are not present in all of the datasets")
#  }
  row.names(pheno.data) <- s.names
#  pheno.data <- pheno.data[,!(colnames(pheno.data) %in% s.id.col)]
  dataset.import <- new("methQTLInput",
    meth.data=meth.data,
    geno.data=geno.data,
    anno.meth=meth.import$annotation,
    anno.geno=geno.import$annotation,
    pheno.data=pheno.data,
    samples=s.names,
    assembly=assembly.meth,
    disk.dump=qtlGetOption("hdf5dump"),
    imputed=geno.import$imputed,
    platform=meth.import$platform,
    segmentation=meth.import$segmentation
  )
  if(assembly.meth != assembly.geno){
    dataset.import <- match.assemblies(dataset.import)
  }
  rm(meth.data)
  rm(geno.data)
  gc()
  logger.completed()
  return(dataset.import)
}

#' doMethImport
#'
#' This function performs import of DNA methylation data and basic preprocessing steps implemented in the \code{\link{RnBeads}} package.
#'
#' @param data.location The location of the data files. See \code{\link{do.import}} for further details.
#' @param assembly The assembly used.
#' @param s.anno Path to the sample annotation sheet
#' @param s.id.col Column name in the sample annotation sheet specifying the sample identifiers.
#' @param tab.sep Table separator used.
#' @param snp.location Locations of SNPs called in the genotyped processing to be removed from the list of CpGs.
#' @param out.folder If not \code{NULL} a directory in which intermediate results are to be written
#' @param ... Further parameters passed to, e.g., \code{\link{qtlRunSegmentation}}, \code{\link{doGenoImportIDAT}}
#' @return A list with five elements:
#'          \describe{
#'            \item{sample}{The samples used in the dataset as a character vector}
#'            \item{data}{The DNA methylation data matrix as the results of import and preprocessing}
#'            \item{annotation}{The genomic annotation of the sites present in \code{data}}
#'            \item{platform}{The platform used, e.g. \code{'probesEPIC'}}
#'            \item{segmentation}{If segmentation was performed, the segmentation as a \code{GRanges} object}
#'          }
#' @details This function execute the import and preprocessing modules of RnBeads. In the default setting, a common
#'          option setting for Illumina BeadArray data is used (described in \code{inst/extdata/rnbeads_options.xml}).
#'          For further option setting, you should use the \code{rnbeads.options} option of the methQTL package to
#'          specify another XML file with custom RnBeads options.
#' @author Michael Scherer
#' @noRd
#' @import RnBeads
doMethImport <- function(data.location,assembly="hg19",s.anno,s.id.col,tab.sep=",",snp.location=NULL,out.folder=NULL,...){
  logger.start("Processing DNA methylation data")
  rnb.xml2options(qtlGetOption("rnbeads.options"))
  rnb.options(identifiers.column=s.id.col,
              assembly=assembly,
              import.table.separator=tab.sep)
  if(qtlGetOption("meth.data.type")%in%c("GEO","rnb.set")){
    data.s <- data.location[["idat.dir"]]
  }else if(qtlGetOption("meth.data.type")=="data.files"){
    data.s <- c(s.anno,data.location[["idat.dir"]])
  }else{
    data.s <- c(data.location[["idat.dir"]],s.anno)
  }
  rnb.imp <- rnb.execute.import(data.source = data.s, data.type = qtlGetOption("meth.data.type"))
  if(qtlGetOption("rnbeads.qc")){
    if(dir.exists(qtlGetOption("rnbeads.report"))){
      rnb.report <- file.path(qtlGetOption("rnbeads.report"),"rnbeads_QC")
    }else if(!is.null(out.folder)){
      rnb.report <- file.path(out.folder,"rnbeads_QC")
    }else{
      rnb.report <- file.path(tempdir(),"rnbeads_QC")
    }
    rnb.run.qc(rnb.imp,rnb.report)
  }
  if(dir.exists(qtlGetOption("rnbeads.report"))){
    rnb.report <- file.path(qtlGetOption("rnbeads.report"),"rnbeads_preprocessing")
  }else if(!is.null(out.folder)){
    rnb.report <- file.path(out.folder,"rnbeads_preprocessing")
  }else{
    rnb.report <- file.path(tempdir(),"rnbeads_preprocessing")
  }
  rnb.imp <- rnb.run.preprocessing(rnb.imp,rnb.report)$rnb.set
  if(!is.null(out.folder)){
    save.rnb.set(rnb.imp,file.path(out.folder,"rnbSet_preprocessed"))
  }
  if(!is.null(snp.location)){
    snp.location <- GRanges(Rle(snp.location$Chromosome),IRanges(start=as.numeric(snp.location$Start),
                                                                 end=as.numeric(snp.location$Start)))
    cpg.annotation <- makeGRangesFromDataFrame(annotation(rnb.imp))
    op <- findOverlaps(cpg.annotation,snp.location,ignore.strand=F,maxgap = 1)
    rem.sites <- rep(FALSE,length(cpg.annotation))
    rem.sites[queryHits(op)] <- TRUE
    logger.start(paste("Removing",sum(rem.sites),"CpGs overlapping with SNPs"))
    rnb.imp <- remove.sites(rnb.imp,rem.sites)
    logger.completed()
  }
  segmentation <- qtlRunSegmentation(rnb.imp,out.folder,...)
  s.names <- samples(rnb.imp)
  meth.data <- meth(rnb.imp)
  if(qtlGetOption("hdf5dump")){
    meth.data <- writeHDF5Array(meth.data)
  }
  anno.meth <- annotation(rnb.imp)
  gc()
  logger.completed()
  return(list(samples=s.names,data=meth.data,annotation=anno.meth,platform=rnb.imp@target,segmentation=segmentation))
}

#' doGenoImport
#'
#' This routine performs data import and processing of genotyping data using the \code{PLINK} program.
#'
#' @param data.location Path to the directory, where the PLINK files are stored. For more information, see \code{\link{doImport}}
#' @param s.anno The sample annotation sheet.
#' @param s.id.col The column name of the sample annotation sheet specifying the sample identifiers.
#' @param out.folder The output folder for storing diagnostic plots.
#' @param ... Futher parameters passed to \code{\link{doGenoImportImputed}} or \code{\link{doGenoImportIDAT}} depending on the
#'          option \code{'geno.data.type'}.
#' @return A list of three elements:
#'        \describe{
#'          \item{data}{The processed genotyping data as a data.frame}
#'          \item{annotation}{The genomic annotation of the SNPs in \code{data}}
#'          \item{pheno.data}{The, potentially modified, sample annotation sheet}
#'        }
#' @author Michael Scherer
#' @export
#' @import data.table
#' @import snpStats
doGenoImport <- function(data.location,s.anno,s.id.col,out.folder,...){
  logger.start("Processing genotyping data")
  snp.loc <- data.location[["geno.dir"]]
  all.files <- list.files(snp.loc,full.names=T)
  data.type <- qtlGetOption("geno.data.type")
  if(data.type=="plink"){
    bed.file <- all.files[grepl(".bed",all.files)]
    bim.file <- all.files[grepl(".bim",all.files)]
    fam.file <- all.files[grepl(".fam",all.files)]
    if(any(c(length(bed.file)==0,length(bim.file)==0,length(fam.file)==0))){
      dos.file <- all.files[grepl(".dos",all.files)]
      id.map <- all.files[grepl(".txt",all.files)]
      if(any(c(length(dos.file)==0,length(id.map)==0))){
        stop("Incompatible input to genotyping processing, needs to be in PLINK format (i.e. .bed, .bim and .fam files)")
      }
      return(doGenoImportImputed(dos.file,id.map,s.anno,s.id.col,out.folder,...))
    }
  }else if(data.type=="idat"){
    res <- doGenoImportIDAT(snp.loc,s.anno,s.id.col,out.folder,...)
    bed.file <- res["bed.file"]
    bim.file <- res["bim.file"]
    fam.file <- res["fam.file"]
  }
  if(qtlGetOption("impute.geno.data")){
    if(any(grepl("_",s.anno[,s.id.col]))){
	stop("Underscores are not allowed in the sampleID if imputation is to be performed")
    }
    res <- doImputation(bed.file,
                         bim.file,
                         fam.file,
                         out.folder)
    bed.file <- res["bed.file"]
    bim.file <- res["bim.file"]
    fam.file <- res["fam.file"]
  }
  snp.dat <- read.plink(bed = bed.file,bim = bim.file,fam = fam.file)
  fam <- snp.dat$fam
  s.anno <- s.anno[as.character(s.anno[,s.id.col]) %in% row.names(fam),]
  if(is.null(s.anno)){
    stop("Sample ids and IDs of PLINK files do not match.")
  }
  keep.frame <- fam[as.character(s.anno[,s.id.col]),c("pedigree","member")]
  keep.file <- file.path(out.folder,"kept_samples.txt")
  proc.data <- file.path(out.folder,"processed_snp_data")
  write.table(keep.frame,keep.file,sep="\t",row.names=F,col.names=F,quote=F)
  plink.file <- gsub(".bed","",bed.file)
  cmd <- paste(qtlGetOption("plink.path"),"--bfile",plink.file,"--keep",keep.file,"--hwe",qtlGetOption("hardy.weinberg.p"),
               "--maf",qtlGetOption("minor.allele.frequency"),"--mind",qtlGetOption("missing.values.samples"),
               "--geno",qtlGetOption("plink.geno"),"--make-bed --out",proc.data)
  system(cmd)
  snp.dat <- read.plink(bed=paste0(proc.data,".bed"),bim=paste0(proc.data,".bim"),fam=paste0(proc.data,".fam"))
  snp.mat <- t(as(snp.dat$genotypes,"numeric"))
  anno.geno <- snp.dat$map
  colnames(anno.geno) <- c("Chromosome","Name","cM","Start","Allele.1","Allele.2")
  ids <- as.character(anno.geno$Name)
  if(!any(grepl("rs*",ids))){
    logger.start("Matching ids in dbSNP")
    if(is.null(qtlGetOption("db.snp.ref"))){
      logger.warning("Please provide a valid path to dbSNP. Cannot map ids to SNPs.")
    }else{
      db.snp <- fread(qtlGetOption("db.snp.ref"))
      db.snp.gr <- makeGRangesFromDataFrame(db.snp,seqnames="#CHROM",start.field="POS",end.field="POS")
      seqlevelsStyle(db.snp.gr) <- "UCSC"
      anno.geno.gr <- makeGRangesFromDataFrame(anno.geno,seqnames="Chromosome",start.field = "Start",end.field = "Start")
      seqlevelsStyle(anno.geno.gr) <- "UCSC"
      op <- findOverlaps(anno.geno.gr,db.snp.gr,select="first")
      logger.info(paste("Replacing",length(op),"IDs with dbSNP IDs."))
      ids[!is.na(op)] <- as.character(db.snp$ID)[op[!is.na(op)]]
      ids <- make.unique(ids)
    }
    logger.completed()
  }
  row.names(anno.geno) <- ids
  anno.geno <- anno.geno[,c("Chromosome","Start","cM","Allele.1","Allele.2")]
  if(!any(grepl("chr*",anno.geno$Chromosome))){
    anno.geno$Chromosome <- paste0("chr",anno.geno$Chromosome)
  }
  logger.start("Computing allele frequencies")
  maj.allele.frequencies <- apply(snp.mat,1,function(x){
    x <- x[!is.na(x)]
    (2*sum(x==0)+sum(x==1))/(2*length(x))
  })
  anno.geno$Allele.1.Freq <- maj.allele.frequencies
  anno.geno$Allele.2.Freq <- 1-maj.allele.frequencies
  logger.completed()
  if(qtlGetOption("recode.allele.frequencies")){
    logger.start("Recoding allele frequencies")
    allele.frequencies <- maj.allele.frequencies<0.5
    allele.frequencies[is.na(allele.frequencies)] <- FALSE
    snp.mat[allele.frequencies,] <- 2-(snp.mat[allele.frequencies,])
    temp <- anno.geno$Allele.2[allele.frequencies]
    anno.geno$Allele.2[allele.frequencies] <- anno.geno$Allele.1[allele.frequencies]
    anno.geno$Allele.1[allele.frequencies] <- temp
    anno.geno$Allele.1.Freq[allele.frequencies] <- 1-maj.allele.frequencies[allele.frequencies]
    anno.geno$Allele.2.Freq[allele.frequencies] <- maj.allele.frequencies[allele.frequencies]
    logger.completed()
  }else{
    snp.mat[snp.mat==2] <- 3
    snp.mat[snp.mat==0] <- 2
    snp.mat[snp.mat==3] <- 0
  }
  snp.mat <- snp.mat[,as.character(s.anno[,s.id.col])]
  pca.obj <- prcomp(t(na.omit(snp.mat)))
  sum.pca <- summary(pca.obj)$importance
  logger.start("Compute genotype PCA")
  to.plot <- data.frame(PC1=pca.obj$x[,1],PC2=pca.obj$x[,2])
  plot <- ggplot(to.plot,aes(x=PC1,y=PC2))+geom_point()+xlab(paste0("PC1 (",round(sum.pca[2,1],3)*100,"% variance explained)"))+
    ylab(paste0("PC2 (",round(sum.pca[2,2],3)*100,"% variance explained)"))+theme_bw()+
    theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
          axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
          axis.text = element_text(size=15,color="black"))
  ggsave(file.path(out.folder,"genetics_PCA.pdf"),plot)
  logger.completed()
  n.comps <- qtlGetOption('n.prin.comp')
  if(!is.null(n.comps)){
    loads <- pca.obj$x
    if(n.comps>ncol(loads)){
      n.comps <- ncol(loads)
      logger.info(paste("Too many PCs specified, setting to maximum number",n.comps))
    }
    ncol.anno <- ncol(s.anno)
    s.anno <- data.frame(s.anno,loads[,1:n.comps])
    colnames(s.anno)[(ncol.anno+1):ncol(s.anno)] <- paste0("PC",1:n.comps)
  }
  if(qtlGetOption("hdf5dump")){
    snp.mat <- writeHDF5Array(snp.mat)
  }
  logger.completed()
  return(list(data=snp.mat,annotation=anno.geno,pheno.data=s.anno,samples=s.anno[,s.id.col],imputed=FALSE))
}

#' doGenoImportImputed
#'
#' This function executes import of imputed genotyping data
#'
#' @param dos.file A dosage file containing the dosage information for the imputed data. Can be gzipped.
#' @param id.map A path to a file mapping the row names of the dosage file to SNP (rs) identifiers.
#' @param s.anno The sample annotation sheet as a data frame.
#' @param s.id.col A column name in the sample annotation sheet specifying the sample identifers.
#' @param out.folder A path to a folder to which the results are to be stored.
#' @param tab.sep The table separator used for \code{dos.file} and \code{id.map}.
#' @return A list of three elements:
#'        \describe{
#'          \item{data}{The processed genotyping data as a data.frame}
#'          \item{annotation}{The genomic annotation of the SNPs in \code{data}}
#'          \item{pheno.data}{The, potentially modified, sample annotation sheet}
#'        }
#' @details Genotyping data can also be imported from imputed genotype data. However, data needs to be already
#'  preprocessed, such that for instance SNPs outside of the Hardy-Weinberg equilibrium, sites with too many
#'  missing values or SNPs with a minimum minor allele frequency are removed/kept, since the data cannot be further
#'  processed by \code{plink}. We expect as input a tabular file that contains SNPs in the rows and samples in the
#'  columns, where each entry is a continous numeric value containing the dosage.
#' @author Michael Scherer
#' @export
doGenoImportImputed <- function(dos.file,
                            id.map,
                            s.anno,
                            s.id.col,
                            out.folder,
                            tab.sep=" "){
  logger.start("Processing imputed data")
  snp.dat <- read.table(dos.file,sep = tab.sep,header=T)
  if(row.names(snp.dat)[1] == "1"){
    row.names(snp.dat) <- as.character(snp.dat[,1])
    snp.dat <- snp.dat[,-1]
  }
  ids.rows <- read.table(id.map,sep=tab.sep,header=T)
  match.ids <- match(row.names(snp.dat),ids.rows$SNP)
  r.names <- ids.rows$rs.id[match.ids]
  match.unique <- match(unique(r.names),r.names)
  snp.dat <- snp.dat[match.unique,]
  anno.geno <- lapply(row.names(snp.dat),function(x)strsplit(x,"_"))
  anno.geno <- t(as.data.frame(anno.geno))
  colnames(anno.geno) <- c("Chromosome","Start","Allele.1","Allele.2")
  anno.geno <- as.data.frame(anno.geno)
  anno.geno$Name <- unique(r.names)
  row.names(snp.dat) <- unique(r.names)
  s.anno <- s.anno[as.character(s.anno[,s.id.col]) %in% colnames(snp.dat),]
  if(is.null(s.anno)){
    stop("Sample ids and IDs of genotypes do not match.")
  }
  # We do not use plink filtering for imputed data
  snp.dat <- snp.dat[,as.character(s.anno[,s.id.col])]
  pca.obj <- prcomp(t(na.omit(snp.dat)))
  sum.pca <- summary(pca.obj)$importance
  logger.start("Compute genotype PCA")
  to.plot <- data.frame(PC1=pca.obj$x[,1],PC2=pca.obj$x[,2])
  plot <- ggplot(to.plot,aes(x=PC1,y=PC2))+geom_point()+xlab(paste0("PC1 (",round(sum.pca[2,1],3)*100,"% variance explained)"))+
    ylab(paste0("PC2 (",round(sum.pca[2,2],3)*100,"% variance explained)"))+theme_bw()+
    theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
         axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
         axis.text = element_text(size=15,color="black"))
  ggsave(file.path(out.folder,"genetics_PCA.pdf"),plot)
  n.comps <- qtlGetOption('n.prin.comp')
  if(!is.null(n.comps)){
    loads <- pca.obj$x
    if(n.comps>ncol(loads)){
      n.comps <- ncol(loads)
      logger.info(paste("Too many PCs specified, setting to maximum number",n.comps))
    }
    ncol.anno <- ncol(s.anno)
    s.anno <- data.frame(s.anno,loads[,1:n.comps])
    colnames(s.anno)[(ncol.anno+1):ncol(s.anno)] <- paste0("PC",1:n.comps)
  }
  logger.completed()
  if(!any(grepl("chr*",anno.geno$Chromosome))){
    anno.geno$Chromosome <- paste0("chr",anno.geno$Chromosome)
  }
  row.names(anno.geno) <- anno.geno$Name
  anno.geno$Start <- as.numeric(as.character(anno.geno$Start))
  # Recode major and minor alleles
  snp.dat <- as.matrix(snp.dat)
  maj.allele.frequencies <- apply(snp.dat,1,function(x){
    x <- x[!is.na(x)]
    (2*sum(x==0)+sum(x==1))/(2*length(x))
  })
  anno.geno$Allele.1.Freq <- maj.allele.frequencies
  anno.geno$Allele.2.Freq <- 1-maj.allele.frequencies
  if(qtlGetOption("recode.allele.frequencies")){
    allele.frequencies <- maj.allele.frequencies<0.5
    allele.frequencies[is.na(allele.frequencies)] <- FALSE
    snp.dat[allele.frequencies,] <- 2-(snp.dat[allele.frequencies,])
    temp <- anno.geno$Allele.2[allele.frequencies]
    anno.geno$Allele.2[allele.frequencies] <- anno.geno$Allele.1[allele.frequencies]
    anno.geno$Allele.1[allele.frequencies] <- temp
    anno.geno$Allele.1.Freq[allele.frequencies] <- 1-maj.allele.frequencies[allele.frequencies]
    anno.geno$Allele.2.Freq[allele.frequencies] <- maj.allele.frequencies[allele.frequencies]
  }
  logger.completed()
   if(qtlGetOption("hdf5dump")){
    snp.dat <- writeHDF5Array(snp.dat)
  }
return(list(data=snp.dat,annotation=anno.geno,pheno.data=s.anno,samples=s.anno[,s.id.col],imputed=TRUE))
}

match.assemblies <- function(meth.qtl){
  print("Not yet implemented")
  return(meth.qtl)
}

#' doGenoImportIDAT
#'
#' This function imports genotyping data from IDAT files using the \code{'crlmm'} R-package
#'
#' @param idat.files A path to the directory where the IDAT files are stored
#' @param s.anno The sample annotation sheet as a \code{'data.frame'}
#' @param s.id.col The column in the sample annotation sheet specifying the sample identifiers
#' @param out.dir The output directory
#' @param store.crlmm Flag indicating, whether the intermediate \code{CNSet} object is to be stored
#'    on disk. On this object, quality control on using the \code{CRLMM} package can be performed.
#' @param idat.platform The array platform to be used. Check those available using the crlmm function \code{\link{validCdfNames}}.
#'          Additionally, for the Illumina OmniExpress 12 v1.0, and Infinium Omni2.5-8 v1.4 we
#'          provide two custom annotations, which can be used using the 'OmniExpress' and
#'          'Omni2-5' options, respectively.
#' @param call.method The genotype calling method passed to \code{\link{genotype.Illumina}}
#' @param gender.col Optional parameter specifying the column name in the sample annotation sheet specifying
#'     the individual's sexes. If not specified, the package will search for such a column.
#' @return A vector with three elements: \describe{
#'   \item{\code{"bed.file"}}{Path to the BED file for PLINK}
#'   \item{\code{"bim.file"}}{Path to the BIM file for PLINK}
#'   \item{\code{"fam.file"}}{Path to the FAM file for PLINK}
#' }
#' @author Michael Scherer
#' @export
#' @import crlmm
#' @import data.table
#' @import methQTL.data
#' @import ff
doGenoImportIDAT <- function(idat.files,
                                   s.anno,
                                   s.id.col,
                                   out.dir,
                                   idat.platform="humanomni258v1a",
                                   call.method="krlmm",
                                   gender.col=NULL,
				   store.crlmm=F){
  logger.start("Importing genotyping IDAT files")
  if(!("GenoSentrixPosition"%in%colnames(s.anno))){
    stop("Missing required column 'GenoSentrixPosition' in the sample annotation sheet")
  }
  array.names <- file.path(idat.files,as.character(s.anno[,"GenoSentrixPosition"]))
  test.exist <- file.exists(paste0(array.names,"_Grn.idat")) & file.exists(paste0(array.names,"_Red.idat"))
  if(any(!test.exist)){
    stop(paste0("Missing idat files, e.g. ",paste(s.anno[!test.exist,"GenoSentrixPosition"],collapse=" ")))
  }
  array.info <- list(barcode=NULL,position="GenoSentrixPosition")
  batch.info <- rep("1",nrow(s.anno))
  loadNamespace("ff")
  if(is.null(gender.col)){
    if(any(c("gender","sex")%in%tolower(colnames(s.anno)))){
      logger.info("Parsing sex information")
      gender.col <- colnames(s.anno)[which(tolower(colnames(s.anno))%in%c("gender","sex"))]
    }else{
      logger.info("No sex information found")
    }
  }
  if(!is.null(gender.col)){
    sex <- s.anno[,gender.col]
    if(any(grepl("f|F",sex))){
      sex <- ifelse(grepl("f|F",sex),2,1)
    }
  }else{
    sex <- rep(NA,nrow(s.anno))
  }
  if(idat.platform=="OmniExpress"){
	if(!requireNamespace("methQTL.data")){
		stop("Missing required package methQTL.data for idat.platform OmniExpress")
	}
	my.anno <- readRDS(system.file("extdata/omni_express_annotation.rds",package="methQTL.data"))
	genome <- "hg19"
        crlmm.obj <- genotype.Illumina(sampleSheet=s.anno,
                         arrayNames=array.names,
                         arrayInfoColNames=array.info,
                         cdfName="nopackage",
			 anno=my.anno,
			 genome=genome,
                         batch=batch.info,
                         call.method=call.method,
                         copynumber=FALSE,
                         fitMixture=FALSE,
                         gender=sex)
	annot <- my.anno@data
  }else if(idat.platform=="Omni2-5"){
	if(!requireNamespace("methQTL.data")){
		stop("Missing required package methQTL.data for idat.platform Omni2-5")
	}
	my.anno <- readRDS(system.file("extdata/omni_25_annotation.rds",package="methQTL.data"))
	genome <- "hg19"
        crlmm.obj <- genotype.Illumina(sampleSheet=s.anno,
                         arrayNames=array.names,
                         arrayInfoColNames=array.info,
                         cdfName="nopackage",
			 anno=my.anno,
			 genome=genome,
                         batch=batch.info,
                         call.method=call.method,
                         copynumber=FALSE,
                         fitMixture=FALSE,
                         gender=sex)
	annot <- my.anno@data
  }else{
	  crlmm.obj <- genotype.Illumina(sampleSheet=s.anno,
		                         arrayNames=array.names,
		                         arrayInfoColNames=array.info,
		                         cdfName=idat.platform,
		                         batch=batch.info,
		                         call.method=call.method,
		                         copynumber=FALSE,
		                         fitMixture=FALSE,
		                         gender=sex)
	  anno.data <- system.file("extdata/annotation.rda",package = paste0(idat.platform,"Crlmm"))
	  if(!file.exists(anno.data)){
	    stop(paste("Missing required annotation data for BeadArray",idat.platform))
	  }
	  load(anno.data)
  }
  if(store.crlmm){
    saveRDS(crlmm.obj,file.path(out.dir,"CNSet.rds"))
  }
  f.dat <- featureData(crlmm.obj)
  assembly <- crlmm.obj@genome
  if(any(!(featureNames(f.dat)%in%annot$Name))){
    stop("Provided array annotation not sufficient, Aborting")
  }
  alleles <- annot$SourceSeq[match(featureNames(f.dat),as.character(annot$Name))]
  pos.mark <- regexpr("\\/",alleles)
  allele.A <- substr(alleles,pos.mark-1,pos.mark-1)
  allele.B <- substr(alleles,pos.mark+1,pos.mark+1)
  chroms <- chromosome(f.dat)
  is.chr.0 <- chroms%in%c(0,25)
  chroms[chroms==23] <- "X"
  chroms[chroms==24] <- "Y"
  snp.mat <- t(calls(crlmm.obj)[,,drop=F])
  snp.mat <- snp.mat[,!is.chr.0]
  chroms <- chroms[!is.chr.0]
  position <- position(f.dat)[!is.chr.0]
  allele.A <- allele.A[!is.chr.0]
  allele.B <- allele.B[!is.chr.0]
  snp.names <- featureNames(f.dat)[!is.chr.0]
  if(!is.null(qtlGetOption("db.snp.ref"))){
    logger.start("Matching reference allele according to dbSNP")
    db.snp <- fread(qtlGetOption("db.snp.ref"))
    db.snp.gr <- makeGRangesFromDataFrame(db.snp,seqnames="#CHROM",start.field="POS",end.field="POS")
    anno.gr <- GRanges(Rle(chroms),IRanges(start=position,end=position))
    op <- findOverlaps(anno.gr,db.snp.gr)
    is.mismatch <- allele.A[queryHits(op)]!=db.snp$REF[subjectHits(op)]
    logger.info(paste("Recoding",sum(is.mismatch),"variants"))
    temp <- allele.A[queryHits(op)][is.mismatch]
    allele.A[queryHits(op)][is.mismatch] <- allele.B[queryHits(op)][is.mismatch]
    allele.B[queryHits(op)][is.mismatch] <- temp
    snp.mat[,queryHits(op)][,is.mismatch] <- 4-(snp.mat[,queryHits(op)][,is.mismatch])
    rem.sites <- rep(FALSE,length(allele.A))
    rem.sites[queryHits(op)] <- allele.A[queryHits(op)]!=db.snp$REF[subjectHits(op)]
    logger.info(paste("Removing",sum(rem.sites),"variants, since neither ALT nor REF matches dbSNP"))
    snp.mat <- snp.mat[,!rem.sites]
    chroms <- chroms[!rem.sites]
    position <- position[!rem.sites]
    allele.A <- allele.A[!rem.sites]
    allele.B <- allele.B[!rem.sites]
    snp.names <- snp.names[!rem.sites]
    logger.completed()
  }
  snp.mat <- new("SnpMatrix",snp.mat)
  ids <- as.character(s.anno[,s.id.col])
  row.names(snp.mat) <- ids
#  snp.dat <- data.frame(chromosome=chroms,
#	  position=position,
#	  genetic.distance=rep(NA,length(chroms)),
#	  allele.1=allele.A,
#	  allele.2=allele.B)
#  row.names(snp.dat) <- snp.names
  write.plink(file.path(out.dir,"plink"),
              snps=snp.mat,
              pedigree=ids,
              id=ids,
              sex=sex,
#              snp.data=snp.dat,
	      chromosome=chroms,
	      position=position,
	      allele.1=allele.A,
	      allele.2=allele.B
              )
  # Sort data by chrosome
  cmd <- paste(qtlGetOption('plink.path'),"--bfile",file.path(out.dir,"plink"),"--make-bed --out",file.path(out.dir,"plink_sorted"))
  system(cmd)
  cmd <- paste("rm -rf",file.path(out.dir,"plink.*"))
  system(cmd)
  logger.completed()
  return(c(
    bed.file=file.path(out.dir,"plink_sorted.bed"),
    bim.file=file.path(out.dir,"plink_sorted.bim"),
    fam.file=file.path(out.dir,"plink_sorted.fam")
  ))
}

#' qtlRunSegmentation
#'
#' This function performs DNA methylation based segmentation using the 'epicPMDdetect' package
#'
#' @param rnb.set An object of type \code{\link{RnBSet-class}} with required DNA methylation information
#' @param out.folder The output folder to store intermediate results
#' @param train.chr The chromosome on which the HMM should be trained
#' @return A \code{GRanges} object with the segmentation performed
#' @details The 'epicPMDdetect' package has been created by Malte Gross
#' @author Michael Scherer
#' @export
qtlRunSegmentation <- function(rnb.set,
				out.folder,
				train.chr="chr2"){
  if(qtlGetOption("use.segmentation")){
    if(requireNamespace("epicPMDdetect")){
      logger.start("Start segmentation")
      gr <- epicPMDdetect::buildMethGrangesFromRnbSet(rnb.set)
      segmentation <- epicPMDdetect::segmentPMDsKNN(gr,training.chr.sel=train.chr)
      segmentation <- segmentation[values(segmentation)$type%in%c("PMD","notPMD")]
      logger.completed()
    }else{
      stop("Please install the 'epicPMDdetect' package, which is required for computing segmentations")
    }
  }else{
    segmentation <- NULL
  }
  return(segmentation)
}
