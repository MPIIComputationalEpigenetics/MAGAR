##########################################################################################
# interpretation.R
# created: 2020-01-31
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# functions to interpret methQTL results including LOLA and GO enrichments
##########################################################################################

#' qtlLOLAEnrichment
#'
#' This functions performn LOLA enrichment analysis for the methQTL sites or the sites that are
#' shared across all methQTLs in the input list.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param lola.db The location of a LOLA DB already downloaded
#' @param assembly The assembly used (defaul \code{'hg19'})
#' @return The LOLA enrichment result
#' @details We use all data points that have been used to calculate methQTLs as the background
#'  and compare the overlaps with the annotation of interest in comparison to the methQTLs that
#'  have been computed in case a \code{\link{methQTLResult-class}} is provided. If a list of \code{\link{methQTLResult-class}} objects
#'  is provided, the intersection between the methQTLs from all objects in the list is compared with the union of all interactions
#'  that have been tested.
#' @author Michael Scherer
#' @export
qtlLOLAEnrichment <- function(meth.qtl.res,type="SNP",lola.db=NULL,assembly="hg19"){
  if(requireNamespace("LOLA")){
    stats <- getOverlapUniverse(meth.qtl.res,type)
    all.input <- stats$all.input
    all.qtl <- stats$all.qtl
    if(is.null(lola.db)){
      lola.db <- downloadLolaDbs(tempdir())[[assembly]]
    }
    lola.db <- loadRegionDB(lola.db)
    lola.res <- runLOLA(userSets=all.qtl,userUniverse = all.input,lola.db)
    return(list(lola.res=lola.res,lola.db=lola.db))
  }
}

#' qtlAnnotationEnrichment
#'
#' This functions performs enrichment analysis using the Fisher's test for the methQTLs detected
#' with respect to different genomic annotations.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects.
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param annotation The genomic annotation to be used. Can be the ones available in \code{\link{rnb.region.types}} or
#'     \code{"ctcf", "distal", "dnase", "proximal", "tfbs", "tss"}
#' @return A list of two p-values named \code{'enrichment'} for overrepresentation and \code{'depletion'} for underrepresentation
#' @details We use all data points that have been used to calculate methQTLs as the background
#'    and compare the overlaps with the annotation of interest in comparison to the methQTLs that
#'    have been computed in case a \code{\link{methQTLResult-class}} is provided. If a list of \code{\link{methQTLResult-class}} objects
#'    is provided, the intersection between the methQTLs from all objects in the list is compared with the union of all interactions
#'    that have been tested.
#' @author Michael Scherer
#' @export
qtlAnnotationEnrichment <- function(meth.qtl.res,
                                      type="SNP",
                                      annotation="cpgislands"){
  ensembl.anno <- c("ctcf","distal","dnase","proximal","tfbs","tss")
  if(!(annotation%in%c("cpgislands","genes","promoters",ensembl.anno))){
    stop(paste("Invalid value for annotation, needs to be",c("cpgislands","genes","promoters",ensembl.anno)))
  }
  if(annotation%in%c("cpgislands","genes","promoters")){
    anno.type <- annotation
  }else{
    anno.type <- paste0("ensembleRegBuildBP",annotation)
    if(!anno.type%in%rnb.region.types()){
      rnb.load.annotation.from.db(anno.type)
    }
  }
  if(!(type%in%c("SNP","CpG","cor.block"))){
    stop("Invalif value for type, needs to be 'SNP','CpG', or 'cor.block'")
  }
  stats <- getOverlapUniverse(meth.qtl.res,type)
  all.input <- stats$all.input
  all.qtl <- stats$all.qtl
  all.input <- all.input[!(names(all.input)%in%names(all.qtl))]
  annotation <- unlist(rnb.get.annotation(anno.type))
  tps <- length(findOverlaps(all.qtl,annotation))
  fps <- length(all.qtl)-tps
  fns <- length(findOverlaps(all.input,annotation))
  tns <- length(all.input)-fns
  gr <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
  le <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="less")$p.value
  or <- (tps/fps)/(fns/tns)
  return(list(enrichment=gr,depletion=le,OddsRatio=or))
}

#' qtlBaseSubstitutionEnrichment
#'
#' This function tests for enrichment of a specific base substitution in the methQTL interactions.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects.
#' @param merge Flag indicating if 5' and 3' substitutions are to be merged or to be analyzed separately.
#' @return A list with one element for each potential base substitution containing the enrichment p-value.
#' @details The names of the list are e.g. \code{A_G}, which refers to a replacement of the reference base \code{A}
#'   with an \code{A}. Enrichment is computed using Fisher's exact test, using all SNP that have been used
#'   as input as the background.
#' @author Michael Scherer
#' @import plyr
#' @export
qtlBaseSubstitutionEnrichment <- function(meth.qtl.res,
					merge=F){
  stats <- getOverlapUniverse(meth.qtl.res,type="SNP")
  if(is.list(meth.qtl.res)){
    anno.cpgs <- overlapInputs(meth.qtl.res,type="SNP")
  }else{
    anno.cpgs <- getAnno(meth.qtl.res,"geno")
  }
  sel.anno <- anno.cpgs[names(stats$all.qtl),]
  subs.qtl <- paste(sel.anno[,"Allele.1"],sel.anno[,"Allele.2"],sep="_")
  sel.anno <- anno.cpgs[names(stats$all.input),]
  subs.input <- paste(sel.anno[,"Allele.1"],sel.anno[,"Allele.2"],sep="_")
  all.interactions <- as.character(sapply(c("A","C","G","T"),function(b){
    paste(c("A","C","G","T")[!c("A","C","G","T")%in%b],b,sep ="_")
  }))
  if(merge){
    subs.map <- c("A_C"="A_C/T_G",
			"T_G"="A_C/T_G",
			"A_G"="A_G/T_C",
			"T_C"="A_G/T_C",
			"A_T"="A_T/T_A",
			"T_A"="A_T/T_A",
			"C_A"="C_A/G_T",
			"G_T"="C_A/G_T",
			"C_G"="C_G/G_C",
			"G_C"="C_G/G_C",
			"C_T"="C_T/G_A",
			"G_A"="C_T/G_A"
			)
    subs.qtl <- subs.map[subs.qtl]
    sub.qtl <- unname(subs.qtl)
    subs.input <- subs.map[subs.input]
    subs.input <- unname(subs.input)
    all.interactions <- subs.map[all.interactions]
    all.interactions <- unique(all.interactions)
    all.interactions <- unname(all.interactions)
  }
  cu.qtl <- count(subs.qtl)
  cu.input <- count(subs.input)
  enr.all <- sapply(all.interactions, function(int){
    tps <- cu.qtl[cu.qtl$x%in%int,"freq"]
    fps <- sum(cu.qtl[!(cu.qtl$x%in%int),"freq"])
    fns <- cu.input[cu.input$x%in%int,"freq"]
    tns <- sum(cu.input[!(cu.input$x%in%int),"freq"])
    p.val <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
    or <- (tps/fps)/(fns/tns)
    return(c(p.value=p.val,OddsRatio=or))
  })
  return(enr.all)
}

#' qtlTFBSMotifEnrichment
#'
#' This function performs TFBS enrichment analysis for the methQTL SNPs/CpGs detected and returns overrepresented
#' binding motifs.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param size Motif enrichment is only supported for genomic regions. Therefore, we resize the invididual methQTL to
#'     genomic regions using a width of this size around the site of interest.
#' @param assembly The assembly used. Only \code{"hg19"} and \code{"hg38"} supported
#' @param subsample Integer specifying how many of the regions are to be subsamples from the universe.
#' @param out.dir The output directory in which resulting plots will be stored.
#' @param ... Further parameters passed to \code{\link[motifRG]{findMotifFgBg}}
#' @return A plot describing the TFB motif enrichment
#' @details This function is in part based on the tutorial for Motif discovery in https://compgenomr.github.io/book/motif-discovery.html.
#' We use all data points that have been used to calculate methQTLs as the background
#'  and compare the overlaps with the annotation of interest in comparison to the methQTLs that
#'  have been computed in case a \code{\link{methQTLResult-class}} is provided. If a list of \code{\link{methQTLResult-class}} objects
#'  is provided, the intersection between the methQTLs from all objects in the list is compared with the union of all interactions
#'  that have been tested.
#' @author Michael Scherer
#' @export
qtlTFBSMotifEnrichment <- function(meth.qtl.res,
	type="SNP",
	size=500,
	assembly="hg19",
	subsample=100000,
	out.dir=getwd(),
	...){
  if(!assembly%in%c("hg19","hg38")){
    stop("Motif enrichment only supported for 'hg19' and 'hg38'")
  }
  if(requireNamespace("motifRG")&requireNamespace("TFBSTools")&requireNamespace("JASPAR2018")){
    stats <- getOverlapUniverse(meth.qtl.res,type)
    all.input <- stats$all.input
    all.input <- resize(all.input,width=size,fix="center")
    all.qtl <- stats$all.qtl
    all.qtl <- resize(all.qtl,width=size,fix="center")
    if(assembly%in%"hg19"){
	all.input <- getSeq(BSgenome.Hsapiens.UCSC.hg19, all.input)
	all.qtl <- getSeq(BSgenome.Hsapiens.UCSC.hg19, all.qtl)
    }else{
	all.input <- getSeq(BSgenome.Hsapiens.UCSC.hg38, all.input)
	all.qtl <- getSeq(BSgenome.Hsapiens.UCSC.hg38, all.qtl)
    }
    all.input <- all.input[sample(1:length(all.input),subsample)]
    motifs <- findMotifFgBg(
      fg.seq = all.qtl,
      bg.seq = all.input,
      ...
    )
    refined.motifs = lapply(motifs$motifs, function(x){
      refinePWMMotifExtend(motifs = x@match$pattern, seqs = all.qtl)
    })
    for(i in 1:length(refined.motifs)){
      mot <- refined.motifs[[i]]
      motif.name <- names(refined.motifs)[i]
      pwm.mot <- PWMatrix(ID='unk',
			profileMatrix=mot$model$prob)
      pwm.mot.lib <- getMatrixSet(
			JASPAR2018,
			opts=list(collection="CORE",species="Homo sapiens",matrixtype="PWM")
		)
      pwm.mot.sim <- PWMSimilarity(pwm.mot.lib,pwm.mot,method="Pearson")
      info.mat <- t(as.data.frame(lapply(pwm.mot.lib,function(x){
         c(ID(x),name(x))
      })))
      info.mat <- data.frame(info.mat,similarity=pwm.mot.sim[info.mat[,1]])
      info.mat <- info.mat[order(info.mat$similarity,decreasing=T),]
      pdf(file.path(out.dir,paste0(motif.name,".pdf")))
      plot.new()
      seqLogo::seqLogo(mot$model$prob)
      text(x=0.5,y=1,label=paste("Most similar to",info.mat$X2[1],"Similarity:",round(info.mat$similarity[1],3)))
      dev.off()
    }
  }
}
