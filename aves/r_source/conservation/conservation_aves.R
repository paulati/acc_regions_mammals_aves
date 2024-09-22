base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'

library(rphast)
source(file.path(base_path, "aves/r_source/common/maf_helper.R"))


required.species.feats <- function(align)
{
  alignment.species <- names.msa(align)
  
  falcon.eagle <- c()
  if ("falChe1" %in% alignment.species) {
    falcon.eagle <- c(falcon.eagle, "falChe1")
  } 
  if("falPer1" %in% alignment.species) {
    falcon.eagle <- c(falcon.eagle, "falPer1")
  }
  if("halLeu1" %in% alignment.species) {
    falcon.eagle <- c(falcon.eagle, "halLeu1")
  }
  
  zebrafinch.medium_ground_finch <- c()
  if("taeGut2" %in% alignment.species) {
    zebrafinch.medium_ground_finch <- c(zebrafinch.medium_ground_finch, "taeGut2")
  }
  if("geoFor1" %in% alignment.species) {
    zebrafinch.medium_ground_finch <- c(zebrafinch.medium_ground_finch, "geoFor1")
  }
  
  ostrich.tinamou <- c()
  if("strCam1" %in% alignment.species) {
    ostrich.tinamou <- c(ostrich.tinamou, "strCam1")
  }
  if("tinGut2" %in% alignment.species) {
    ostrich.tinamou <- c(ostrich.tinamou, "tinGut2")
  }
  
  #required 11 birds
  
  # #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
  req.regions.1 <- informative.regions.msa(align, min.numspec=1, spec=falcon.eagle, 
                                        refseq="galGal6", gaps.inf=FALSE) 

  req.regions.2 <- informative.regions.msa(align, min.numspec=1, spec=zebrafinch.medium_ground_finch, 
                                           refseq="galGal6", gaps.inf=FALSE) 
  
  req.regions.3 <- informative.regions.msa(align, min.numspec=1, spec=ostrich.tinamou, 
                                           refseq="galGal6", gaps.inf=FALSE) 
  
  req.regions.4 <- informative.regions.msa(align, min.numspec=11, spec=NULL, 
                                           refseq="galGal6", gaps.inf=FALSE)
  
  
  # spec: the default value of NULL implies use all species in the alignment.
  #req.regions.2 <- informative.regions.msa(align.2, min.numspec=10, spec=NULL, 
  #                                       refseq="oryzias_latipes", gaps.inf=FALSE) 
  
  #coverage.feat: Any features object passed into this function which is stored as a pointer 
  #to an object stored in C may be reordered (sorted) by this function.
  # result <- coverage.feat(req.regions, get.feats=TRUE)
  
  intersection.regions <- coverage.feat(req.regions.1, 
                                     req.regions.2, 
                                     req.regions.3, 
                                     req.regions.4, 
                                     or=FALSE, get.feats=TRUE)
  
  
  
  
  intersection.regions.order <- sort(intersection.regions, decreasing = FALSE)
  
  return(intersection.regions.order)
  
  
  # Primate subset: al menos hg38 y otoGar3. (No sé si incorporar a otro primate, como Rhesus)
  # Euarchontoglires subset: al menos mm10, oryCun2, ochPri3
  # Laurasiatheria subset: al menos susScr3, turTru2, bosTau8, felCat8, myoLuc2
  # Afrotheria subset: al menos loxAfr3 o echTel2 
  # Mammal subset: al menos dasNov3 o monDom5 o macEug2, ornAna1.
  
  # cons.req.primates <- c("hg38", "otoGar3" )
  # cons.req.euarchontoglires <- c("mm10", "oryCun2", "ochPri3")
  # cons.req.laurasiatheria <- c("susScr3", "turTru2", "bosTau8", "felCat8", "myoLuc2")
  # cons.req.mammal <- c("ornAna1")
  # 
  # cons.req <- c(cons.req.primates, cons.req.euarchontoglires, cons.req.laurasiatheria, cons.req.mammal)
  # 
  # #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
  # req.regions.1 <- informative.regions.msa(align, min.numspec=11, spec=cons.req, 
  #                                          refseq="hg38", gaps.inf=FALSE)    
  # 
  # cons.opc.afrotheria <- c("loxAfr3", "echTel2")
  # req.regions.2 <- informative.regions.msa(align, min.numspec=1, spec=cons.opc.afrotheria, 
  #                                          refseq="hg38", gaps.inf=FALSE)    
  # 
  # cons.opc.mammal <- c("dasNov3", "monDom5", "macEug2")
  # req.regions.3 <- informative.regions.msa(align, min.numspec=1, spec=cons.opc.mammal, 
  #                                          refseq="hg38", gaps.inf=FALSE)    
  # 
  # 
  # #coverage.feat: Any features object passed into this function which is stored as a pointer 
  # #to an object stored in C may be reordered (sorted) by this function.
  # intersection.regions <- coverage.feat(req.regions.1, 
  #                                       req.regions.2, 
  #                                       req.regions.3, 
  #                                       or=FALSE, get.feats=TRUE)
  # 
  # intersection.regions.order <- sort(intersection.regions, decreasing = FALSE)
  # 
  # return(intersection.regions.order)
  
}


main <- function()
{
  
  config <- environment.config()
  
  for(chr.id in c(30:30))
  {
    print(chr.id)
    
    download.alignment.from.S3(chr.id)
    
    align <- load.alignment(chr.id)
    
    # "galGal6",  
    # species <- c("tinGut2", "picPub1", "araMac1", "pseHum1", "taeGut2", "serCan1", "corBra1", 
    #              "halLeu1", "falPer1", "aptFor1", "phaCar1", "phoRub1", "cucCan1", "tytAlb1", 
    #              "anaPla1", "opiHoa1", "nipNip1", "gavSte1", "eurHel1", "strCam1", "geoFor1" )
    
    # exclude.species <- paste(c("galGal6.chr"), c(1:29,30:33), sep = "")
    
    #align.2 <- load.alignment(chr.id, species)
    
    #libero mem, rompe con chr 1 y 2
    gc()
    
    required.feats <- required.species.feats(align)
    
    neutral.model.base.path <- config$neutral_model$model_base_path
    neutral.model.file.name <- config$neutral_model$model_global_file_name
    setwd(neutral.model.base.path)
    neutralMod <- read.tm(neutral.model.file.name)
    
    #libero mem, rompe con chr 1 y 2
    gc()
    
    elements <- phastCons(align, neutralMod, expected.length=45,
                          target.coverage=0.3, rho=0.3, viterbi=TRUE)
    
    #cons.elements es un objeto de tipo feat
    cons.elements <- elements$most.conserved
    
    cons.reqs <- coverage.feat(required.feats, cons.elements, or=FALSE, get.feats=TRUE)
    
    # q <- nrow(intersection.cons.reqs)
    q <- nrow(cons.reqs)
    bed.data <- data.frame(chr=character(q),
                           start=integer(q),
                           end=integer(q))
    
    char.name <- paste("chr", chr.id, sep="")
    
    bed.data$chr <- rep(char.name, q)
    bed.data$start <- cons.reqs$start
    bed.data$end <- cons.reqs$end
    phastConsElements.out.base.path <- config$phastCons$local_base_path
    out.file.name <- paste0(config$phastCons$local_file_name_head, chr.id, config$phastCons$local_file_name_tail)
    setwd(phastConsElements.out.base.path)
    write.table(bed.data, out.file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    
    
    rm()
    gc()
    
  }  
  
}