base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'
# source("/home/rstudio/2018/common/config/paths_config_modneuPhastCons100way.R")
source(file.path(base_path, 'aves/r_source/common/format_files.R'))
source(file.path(base_path, 'aves/r_source/common/maf_helper.R'))

run.obsPhyloP <- function(align, neutral.mod, split.feats) {
  
  #align.copy <- copy.msa(align)
  
  obsPhyloP <- phyloP(neutral.mod, 
                      msa=align, 
                      mode="ACC",
                      features=split.feats, 
                      branches = "AVES")    
  
  return(obsPhyloP)
  
}


main <- function(chr.id, split.length)
{
  
  config <- environment.config("common")
  
  print(chr.id)
  
  # download.alignment.from.S3(chr.id)
    
  align <- load.alignment(chr.id, NULL, "aves")
  
  # file.name.tail <- "_mammals.maf"
  # remote.base.folder.path <- remote.mammals.align.base.folder.path
  # local.base.folder.path <- mammals.align.base.path
  # 
  # align <- load.msa(chr.id, file.name.tail, 
  #                   remote.base.folder.path, 
  #                   local.base.folder.path)
  
  # local.file.name.gz <- paste("chr", chr.id, "_mammals.maf.gz", sep="")
  # local.file.name <- paste("chr", chr.id, "_mammals.maf", sep="")
  # align <- read.msa.local(base.path, local.file.name.gz, local.file.name)
    
  # setwd(neutral.model.base.path)
  # neutral.model.file.name <- "hg38.phastCons100way_label.mod"
  # neutral.mod <- read.tm(neutral.model.file.name)
  
  neutral.model.base.path <- config$neutral_model$model_base_path
  neutral.model.file.name <- config$neutral_model$model_global_file_name
  setwd(neutral.model.base.path)
  neutral.mod <- read.tm(neutral.model.file.name)
  
  phastConsElements.base.path <- config$phastConsIntersection$local_base_path
  phastCons.file.name <- paste0(config$phastConsIntersection$local_file_name_head, chr.id, config$phastConsIntersection$local_file_name_tail)
  setwd(phastConsElements.base.path)
  feats.bed <- read.table(phastCons.file.name, sep="\t")
  colnames(feats.bed) <- c("chr", "start", "end")
  seq.name <- paste("chr", chr.id, sep="")
  feats <- bed.to.feat(seq.name, feats.bed)
  
  # setwd(conservation.intersection.output.base.path)
  # feat.file.name <- paste("chr", chr.id, "_intersection_sarcopterygii_mammals.bed", sep="")
  # feats.bed <- read.table(feat.file.name, sep=" ")
  # colnames(feats.bed) <- c("chr", "start", "end")
  # seq.name <- paste("chr", chr.id, sep="")
  # feats <- bed.to.feat(seq.name, feats.bed)
  
  split.feats.base.path <- config$acceleration$split_feats$local_base_path
  split.feats.file.name.head <- config$acceleration$split_feats$local_file_name_head
  split.feats.file.name.tail <- config$acceleration$split_feats$local_file_name_tail
  
  file.name <- paste0(split.feats.file.name.head, chr.id, "_", split.length, split.feats.file.name.tail)
  #split de feats de long split.length
  setwd(split.feats.base.path)
  
  if(file.exists(file.name))
  {
    split.feats <- read.table(file.name, sep="\t", header = TRUE)
  }  else  {
    split.feats <- split.feat(feats, f=split.length, drop=TRUE)
    write.table(split.feats, file.name, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
  phyloP.result <- run.obsPhyloP(align, neutral.mod, split.feats) 
  
  #ordeno los resultados de phyloP por score
  phyloP.result.sort <- phyloP.result[order(-phyloP.result$score), ]
  
  #escribo los resultados de obs.phyloP:
  obsphyloP.output.base.path <- config$acceleration$obs_phyloP$local_base_path
  out.file.name <- paste0(config$acceleration$obs_phyloP$local_file_name_head, chr.id, 
                          "_", split.length, config$acceleration$obs_phyloP$local_file_name_tail)
  
  setwd(obsphyloP.output.base.path)
  write.table(phyloP.result.sort, out.file.name, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

}

#for(chr.id in c(1:28, 30:33))
for(chr.id in c(31:31))
{
  download.alignment.from.S3(chr.id)
  
  #split.length <- 25
  #main(chr.id, split.length)
  #gc()

  split.length <- 50
  main(chr.id, split.length)
  gc()
}





