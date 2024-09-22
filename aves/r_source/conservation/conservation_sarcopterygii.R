base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'

library(rphast)
source(file.path(base_path, "aves/r_source/common/maf_helper.R"))


required.species.feats <- function(align)
{
  
  alignment.species <- names.msa(align)
  
  cons.req.species <- c()
  if ("allMis1" %in% alignment.species) {
    cons.req.species <- c(cons.req.species, "allMis1")
  } 
  if ("anoCar2" %in% alignment.species) {
    cons.req.species <- c(cons.req.species, "anoCar2")
  } 
  
  cons.req.turtles <- c()
  if ("cheMyd1" %in% alignment.species) {
    cons.req.turtles <- c(cons.req.turtles, "cheMyd1")
  } 
  if ("chrPic2" %in% alignment.species) {
    cons.req.turtles <- c(cons.req.turtles, "chrPic2")
  } 
  if ("pelSin1" %in% alignment.species) {
    cons.req.turtles <- c(cons.req.turtles, "pelSin1")
  } 
  if ("apaSpi1" %in% alignment.species) {
    cons.req.turtles <- c(cons.req.turtles, "apaSpi1")
  } 
  
  # cons.req.species <- c("allMis1", "anoCar2")
  #cons.req.turtles <- c("cheMyd1", "chrPic2", "pelSin1", "apaSpi1")
  #esta no esta en el alineamiento 
  
  #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
  req.regions <- informative.regions.msa(align, min.numspec=2, spec=cons.req.species, 
                                         refseq="galGal6", gaps.inf=FALSE)    
  
  #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.  
  turtle.regions <- informative.regions.msa(align, min.numspec=1, spec=cons.req.turtles, 
                                            refseq="galGal6", gaps.inf=FALSE)    
  
  #calculo la interseccion entre req.regions y turtle.regions
  #coverage.feat: Any features object passed into this function which is stored as a pointer 
  #to an object stored in C may be reordered (sorted) by this function.
  intersection.regions <- coverage.feat(req.regions, turtle.regions, or=FALSE, get.feats=TRUE)
  
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
  
  for(chr.id in c(30:33))
  {
    print(chr.id)
    
    download.alignment.from.S3(chr.id)
    
    align <- load.alignment(chr.id)
    
    # "astyanax_mexicanus", 
    # species <- c("oryzias_latipes", "oryzias_latipes_hni", "oryzias_melastigma",
    #               "gambusia_affinis", "xiphophorus_maculatus", "poecilia_reticulata",
    #               "amphiprion_percula", "astatotilapia_calliptera", "maylandia_zebra",
    #               "cynoglossus_semilaevis", "scophthalmus_maximus", "gasterosteus_aculeatus", 
    #               "takifugu_rubripes", "tetraodon_nigroviridis", "seriola_dumerili", 
    #               "mastacembelus_armatus", "xiphophorus_couchianus", "anabas_testudineus",
    #               "oryzias_latipes_hsok", "danio_rerio", 
    #               "lepisosteus_oculatus", "mola_mola" )
    # # exclude.species <- paste(c("astyanax_mexicanus"), c(1:24), sep = ".")
    # align.2 <- load.alignment(chr.id, species)
    # 
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