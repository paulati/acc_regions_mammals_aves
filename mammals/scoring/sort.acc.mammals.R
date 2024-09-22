source("/home/rstudio/2018/scoring/scoring.config.mammals.R")

concat.acc.elements <- function(chr){

  mammals.25.acc.file.name <- paste(mammals.25.acc.file.name.head, chr, mammals.25.acc.file.name.tail, sep="")
  setwd(mammals.25.acc.base.path)
  data.25.mammals.acc <- read.table(mammals.25.acc.file.name, sep="\t")
  colnames(data.25.mammals.acc) <- c("chr", "start", "end")
  
  mammals.25.obsPhyloP.file.name <- paste(mammals.25.obsPhyloP.file.name.head, chr, mammals.25.obsPhyloP.file.name.tail, sep="")
  setwd(mammals.25.obsPhyloP.base.path)
  data.25.mammals.obsPhyloP <- read.table(mammals.25.obsPhyloP.file.name, sep="\t", header = TRUE)
  
  acc.25.obs.lnlrt <- merge(data.25.mammals.obsPhyloP, data.25.mammals.acc, by=c("chr", "start", "end"))
  
  #acc.25.obs.lnlrt.sort <- acc.25.obs.lnlrt[order(acc.25.obs.lnlrt$lnlratio), ]
  #acc.25.obs.lnlrt.sort.bed <- acc.25.obs.lnlrt.sort[, c("chr", "start", "end")]
  
  mammals.50.acc.file.name <- paste(mammals.50.acc.file.name.head, chr, mammals.50.acc.file.name.tail, sep="")
  setwd(mammals.50.acc.base.path)
  data.50.mammals.acc <- read.table(mammals.50.acc.file.name, sep="\t")
  colnames(data.50.mammals.acc) <- c("chr", "start", "end")
  
  mammals.50.obsPhyloP.file.name <- paste(mammals.50.obsPhyloP.file.name.head, chr, mammals.50.obsPhyloP.file.name.tail, sep="")
  setwd(mammals.50.obsPhyloP.base.path)
  data.50.mammals.obsPhyloP <- read.table(mammals.50.obsPhyloP.file.name, sep="\t", header = TRUE)
  
  acc.50.obs.lnlrt <- merge(data.50.mammals.obsPhyloP, data.50.mammals.acc, by=c("chr", "start", "end"))
  
  #acc.50.obs.lnlrt.sort <- acc.50.obs.lnlrt[order(acc.50.obs.lnlrt$lnlratio), ]
  #acc.50.obs.lnlrt.sort.bed <- acc.50.obs.lnlrt.sort[, c("chr", "start", "end")]

  acc.25.50.obs.lnlrt <- rbind(acc.25.obs.lnlrt, acc.50.obs.lnlrt)
  
  return (acc.25.50.obs.lnlrt)
  
}


chr <- 1
result <- concat.acc.elements(chr)

for(chr in 2:22)
{
  tmp <- concat.acc.elements(chr)
  result <- rbind(result, tmp)
}

result.sort <- result[order(result$lnlratio), ]
result.sort.bed <- result.sort[, c("chr", "start", "end")]
setwd(mammals.scoring.base.path)
write.table(result.sort.bed, "mammals_acc.bed", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)


