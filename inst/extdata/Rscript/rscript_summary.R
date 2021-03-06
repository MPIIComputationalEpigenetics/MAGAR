suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(MAGAR))

ap <- ArgumentParser()
ap$add_argument("-o","--output",action="store",help="Output directory")
cmd.args <- ap$parse_args()

logger.start("Combining results")
methQTL.result.files <- list.files(cmd.args$output,pattern = "MethQTLResult_",full.names = TRUE)
methQTL.results <- lapply(methQTL.result.files,loadMethQTLResult)
methQTL.results <- joinMethQTLResult(methQTL.results)
unlink(methQTL.result.files,recursive = TRUE)
logger.completed()

logger.start("Saving results")
path.save <- file.path(cmd.args$output,paste0("MethQTLResult"))
saveMethQTLResult(methQTL.results,path.save)
logger.completed()
