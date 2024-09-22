library(GenomicRanges)
library(GenomicFeatures)
#install.packages("RMariaDB")
library(RMariaDB)

library(rtracklayer)
ucscGenomes()[ , "db"]

## Display the list of tables supported by makeTxDbFromUCSC():
#supportedUCSCtables()
supportedUCSCtables(genome="galGal6")

## Retrieve a full transcript dataset for Yeast from UCSC:
txdb1 <- makeTxDbFromUCSC(genome="galGal6", tablename="ncbiRefSeqCurated")
txdb1



library(AnnotationDbi)
db_base_path <- "/home/rstudio/2019/aves/data/functional_regions"
db_file_name <- "galGal6_ncbiRefSeqCurated.sqlite"
setwd(db_base_path)

saveDb(txdb1, db_file_name) 

#testeo la base:
txdb <- loadDb(db_file_name)
exons(txdb)
exonicParts(txdb)

