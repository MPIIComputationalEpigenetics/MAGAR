##########################################################################################
# imputation.R
# created: 2020-04-22
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods for imputing genotype data using the Michigan imputation server.
##########################################################################################

#' do.imputation
#' 
#' Function to perform imputation from PLINK data
#' 
#' @param bed.file Path to the PLINK BED file
#' @param bim.file Path to the PLINK BIM file
#' @param fam.file Path to the PLINK FAM file
#' @param out.dir The output directory
#' @return A vector with three elements: \describe{
#'   \item{\code{"bed.file"}}{Path to the BED file for PLINK}
#'   \item{\code{"bim.file"}}{Path to the BIM file for PLINK}
#'   \item{\code{"fam.file"}}{Path to the FAM file for PLINK}
#' }
#' @author Michael Scherer
#' @export
do.imputation <- function(bed.file,
                          bim.file,
                          fam.file,
                          out.dir){
  snp.dat <- read.plink(bed = bed.file,bim = bim.file,fam = fam.file)
  all.chroms <- unique(snp.dat$map$chromosome)
  if(!any(grepl("chr",all.chroms))){
    all.chroms <- paste0("chr",all.chroms)
  }
  all.chroms.files <- file.path(out.dir,all.chroms)
  if(!all(file.exists(all.chroms.files))){
    sapply(all.chroms.files,dir.create)
  }
  proc.data <- file.path(out.dir,"imputed_data")
  plink.loc <- gsub(".bed","",bed.file)
  cmd <- paste(qtl.getOption("plink.path"),"--bfile",plink.loc,
               "--hwe",qtl.getOption("hardy.weinberg.p"),
               "--maf",qtl.getOption("minor.allele.frequency"),"--mind",qtl.getOption("missing.values.samples"),
               "--make-bed --recode vcf --out",proc.data)
  system(cmd)
  # sort and split by chromosome
  vcftools.path <- qtl.getOption("vcftools.path")
  if(is.null(vcftools.path)){
    stop("Path to a function version of VCFtools needs to be specified using the option 'vcftools'")
  }
  cmd <- paste("perl",file.path(vcftools.path,"vcf-sort"),paste0(proc.data,".vcf"),">",paste0(proc.data,"_sorted.vcf; rm -rf"),paste0(proc.data,".vcf"))
  system(cmd)
  cmd <- paste(qtl.getOption("bgzip.path"),paste0(proc.data,"_sorted.vcf"),"\n",qtl.getOption("tabix.path"),"-p vcf",paste0(proc.data,"_sorted.vcf.gz"))
  system(cmd)
  for(chr in all.chroms){
    cmd <- paste(qtl.getOption("tabix.path"),"-h",paste0(proc.data,"_sorted.vcf.gz"),gsub("chr","",chr),">",paste0(out.dir,"/",chr,".vcf"))
    system(cmd)
    cmd <- paste(qtl.getOption("bgzip.path"),paste0(out.dir,"/",chr,".vcf"))
    system(cmd)
    cmd <- paste0("curl -H 'X-Auth-Token: ",
        qtl.getOption("imputation.user.token"),
        "' -F 'input-files=@",paste0(out.dir,"/",chr,".vcf.gz'"),
        " -F 'input-refpanel=",qtl.getOption("imputation.reference.panel"),
        "' -F 'input-phasing=", qtl.getOption("imputation.phasing.method"),
        "' https://imputationserver.sph.umich.edu/api/v2/jobs/submit/minimac4")
    res <- system(cmd,intern=T)
    res <- unlist(strsplit(res,"\""))
    j.id <- res[unlist(lapply(res,function(x)grepl("job-",x)))]
    run <- T
    logger.start("Waiting for imputation jobs to finish. This can take up to several days and you can check the state of the jobs in you user profile on https://imputationserver.sph.umich.edu/index.html#!pages/jobs")
    while(run){
      Sys.sleep(100)
      cmd <- paste0("curl -H 'X-Auth-Token:",
      qtl.getOption("imputation.user.token"), "' https://imputationserver.sph.umich.edu/api/v2/jobs/",j.id,"/status")
      rep <- system(cmd,intern=T)
      rep <- fromJSON(rep)
      run <- !rep$complete
    }
    logger.completed()
    imp.resu <- paste0("curl -H 'X-Auth-Token:",
      qtl.getOption("imputation.user.token"),"' https://imputationserver.sph.umich.edu/results/",j.id,"/local/",chr,".zip --output temp.zip")
    imp.resu <- system(imp.resu)
    cat("Enter password for zip archive send per mail ")
    pwd <- readLines("stdin",n=1)
    cmd <- paste("unzip","-P",pwd,"temp.zip")
    system(cmd)
    cmd <- paste(qtl.getOption("plink.path"),"--vcf",paste0(chr,".dose.vcf.gz"),"--recode --make-bed --out",paste0(out.dir,"/",chr))
    system(cmd)
    write.table(paste0(out.dir,"/",chr,c(".bed",".bim",".fam"),collapse="\t"),file.path(out.dir,"allfiles.txt"),append=T,col.names=F,row.names=F,quote=F)
  }
  all.files <- read.table(file.path(wd,"allfiles.txt"))
  for(fi in 1:nrow(all.files)){
    fi <- all.files[fi,]
    cmd <- paste(qtl.getOption("plink.path"),"--bfile",file.path(out.dir,unlist(fi[1])),
      "--exclude", file.path(out.dir,"all_dups.dups"),"--make-bed --out",file.path(out.dir,paste0(unlist(fi[1]))))
    system(cmd)
    cmd <- paste(qtl.getOption("plink.path"),"--bfile",file.path(out.dir,paste0(unlist(fi[1]))),
                 " --flip", file.path(out.dir,"imputed_data-merge.missnp"),"--make-bed --out",file.path(out.dir,unlist(fi[1])))
    system(cmd)
  }
  
  cmd <- paste(qtl.getOption("plink.path"),"--merge-list", file.path(out.dir,"allfiles.txt"),"--missing-genotype N --make-bed --out",proc.data)
  system(cmd)
  
  proc.data <- file.path(out,"imputed_data")
  return(c(
    bed.file=paste0(proc.data,".bed"),
    bim.file=paste0(proc.data,".bim"),
    fam.file=paste0(proc.data,".fam"),
  ))
}