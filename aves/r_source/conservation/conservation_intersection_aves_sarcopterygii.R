base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'

library(rphast)

source(file.path(base_path, "aves/r_source/common/maf_helper.R"))
source(file.path(base_path, "aves/r_source/common/format_files.R"))

# source("~/2018/common/config/paths_config_modneuPhastCons100way.R")


load.data <- function(chr.id)
{
  config <- environment.config()
  conservation.sarcopterygii.base.path <- config$phastCons$local_base_path
  conservation.aves.base.path <- config$phastCons$local_base_path
  
  setwd(conservation.sarcopterygii.base.path)
  file.name <- paste("chr", chr.id, "_sarcopterygii.bed", sep="")
  data.sarcopterygii <- read.table(file.name, sep="\t")
  colnames(data.sarcopterygii) <- c("chr", "start", "end")
  
  setwd(conservation.aves.base.path)
  file.name <- paste("chr", chr.id, "_aves.bed", sep="")
  data.aves <- read.table(file.name, sep="\t")
  colnames(data.aves) <- c("chr", "start", "end")
  
  result <- list(sarcopterygii=data.sarcopterygii,
                 aves = data.aves)
  
  return(result)
    
}


main <- function() {

  config <- environment.config()
  
  #chr.id <- 22
  for(chr.id in c(30:33))
  {
    data <- load.data(chr.id)
  
    seq.name <- paste('chr', chr.id, sep="")
    
    feat.1 <- bed.to.feat(seq.name, data$sarcopterygii)
    feat.2 <- bed.to.feat(seq.name, data$aves)
    
    intersection.regions <- coverage.feat(feat.1, feat.2, 
                                          or=FALSE, get.feats=TRUE)
    
    output <- feat.to.bed(chr.id, intersection.regions)
    
    
    conservation.intersection.output.base.path <- config$phastConsIntersection$local_base_path
    file.name <- paste0(config$phastConsIntersection$local_file_name_head, chr.id, config$phastConsIntersection$local_file_name_tail)

    setwd(conservation.intersection.output.base.path)
    #file.name <- paste("chr", chr.id, "_intersection_sarcopterygii_aves.bed", sep="")
    
    write.table(output, file.name, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  
  }
  
}