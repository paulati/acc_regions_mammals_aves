
base.path <- "/paula/2018/acc_regions/scoring/data/norm_lnlrt_filter"
setwd(base.path)

file.25 <- "acc_25.bed"
file.50 <- "acc_50.bed"

data.25 <- read.table(file.25, sep = "\t")
data.50 <- read.table(file.50, sep = "\t")

library(GenomicRanges)

ir.25 <- IRanges(start = data.25$V2, end = data.25$V3, names = data.25$V1)
gr.25 <- GRanges(seqnames = names(ir.25), ir.25)
ir.50 <- IRanges(start = data.50$V2, end = data.50$V3, names = data.50$V1)
gr.50 <- GRanges(seqnames = names(ir.50), ir.50)

gr.all <- c(gr.25, gr.50)
gr.all.clean <- reduce(gr.all)

length(ir.25)
length(ir.50)
length(ir.all) == length(ir.25) + length(ir.50)
length(ir.all)
length(ir.all.clean)

data.acc <- as.data.frame(gr.all.clean)[, 1:3]

write.table(data.acc, "/paula/2018/acc_regions/scoring/data/norm_lnlrt_filter/acc_all.bed", 
            sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
