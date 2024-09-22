library(yaml)
base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'
source(file.path(base_path, '/aves/r_source/common/aws_base.R'))

environment.config <- function(species.group = 'aves'){
  
  config.file.name <- paste0('config_', species.group,'.yaml')
  config.file.path <- file.path(base_path, 'aves', 'r_source', 'config', 
                                config.file.name)
  # config <- read_yaml('/home/rstudio/2019/aves/r_source/config/config_aves.yaml')
  #config <- read_yaml('/home/rstudio/2019/aves/r_source/config/config_sarcopterygii.yaml')
  config <- read_yaml(config.file.path)
  env.name <- config$environment
  env.config <- config[[env.name]]
  return(env.config)
}



download.alignment.from.S3 <- function(chr){
  
  config <- environment.config()
  
  account.key <- config$aws_s3$key
  account.secret <- config$aws_s3$secret
  remote.base.folder.path <- config$aws_s3$remote_maf_base_path
  file.name.gz <- paste0(config$aws_s3$remote_maf_file_name_head, chr, config$aws_s3$remote_maf_file_name_tail)
  local.base.folder.path <-  config$maf$local_base_path
  local.file.name.gz <- file.name.gz
  local.file.name <- paste0(config$maf$local_maf_file_name_head, chr, config$maf$local_maf_file_name_tail)
  bucket.name <- config$aws_s3$bucket_name
  
  download.from.s3(account.key, account.secret,
                   remote.base.folder.path, file.name.gz,
                   local.base.folder.path, local.file.name.gz,
                   local.file.name, bucket.name)
  
  
}


load.alignment <- function(chr, species = NULL, config.branch = NULL) {
  
  if (is.null(config.branch)) {
    config <- environment.config()  
  } else {
    config <- environment.config(config.branch)  
  }
  
  align.file.base.folder <- config$maf$local_base_path
  align.file.name <- paste0(config$maf$local_maf_file_name_head, chr, config$maf$local_maf_file_name_tail)
  
  setwd(align.file.base.folder)
  
  if(! file.exists(align.file.name)){
    download.alignment.from.S3(chr)  
  }
  
  if(is.null(species)) {
    align <- read.msa(align.file.name, pointer.only=TRUE)
  } else {
    align <- read.msa(align.file.name, seqnames=species, pointer.only=TRUE)
  }
  
  # 
  # file.name.tail <- "_mammals.maf"
  # remote.base.folder.path <- remote.mammals.align.base.folder.path
  # local.base.folder.path <- mammals.align.base.path
  # # align <- load.alignment.base(chr.id, file.name.tail, 
  # #                              remote.base.folder.path, 
  # #                              local.base.folder.path)
  # 
  # align <- load.msa(chr.id, file.name.tail, 
  #                   remote.base.folder.path, 
  #                   local.base.folder.path)
  
  return(align)  
  
}

