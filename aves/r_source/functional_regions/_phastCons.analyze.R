base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'

library(GenomicRanges)

collapse.neighbors <- function(gr, distance.threshold, length.threshold) {
  
  gr.reduced <- reduce(gr, drop.empty.ranges=FALSE, min.gapwidth=distance.threshold, 
                          with.revmap=FALSE, with.inframe.attrib=FALSE)
  # min.gapwidth Ranges separated by a gap of at least min.gapwidth positions are not merged
  
  gr.reduced.length.filter <- gr.reduced[width(gr.reduced) >= length.threshold]
  
  return(gr.reduced.length.filter)
  
}

load.phastCons <- function(base.path) {
  #setwd(base.path)
  #chrs <- c(1:22)

  bed_files <- list.files(base.path, full.names = TRUE)
  
  result.df <- data.frame(chr = character(),
                       start = numeric(),
                       end = numeric(),
                       stringsAsFactors = FALSE)
  
  for(file.name in bed_files) {
  
    #file.name <- paste0("chr", chr, "_intersection_sarcopterygii_mammals.bed")    
    data <- read.table(file.name, sep="\t", header=FALSE, stringsAsFactors = FALSE)  
    colnames(data) <- c("chr", "start", "end")
    result.df <- rbind(result.df, data)
  
  }
  
  result <- GRanges(result.df$chr, IRanges(start = result.df$star, end = result.df$end))
  return(result)
  
}

build.functional.elements <- function(phastCons.base.paths) {

  #data.base.paths <- "/paula/2019/aves/data/phastConsSarcopterygiiAves"

  gr <- load.phastCons(phastCons.base.paths)

  distance.threshold <- 20

  length.threshold <- 100

  gr.util <- collapse.neighbors(gr, distance.threshold, length.threshold)

  gr.util$id <- c(1:length(gr.util))
  
  gr.util.df <- as.data.frame(gr.util)
  
  out.base.path <- file.path(base_path, 'aves', 'data')
  out.file.path <- file.path(out.base.path, "intersection_sarcopterygii_aves_functional_elements.bed")
  write.table(gr.util.df, out.file.path, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)  
  
  return(gr.util)
  
}


# start

gr.functional.elements <- build.functional.elements(file.path(base_path, "aves/data/phastCons"))
functional.elements.file.path <- file.path(base_path, "aves/data/intersection_sarcopterygii_aves_functional_elements.bed")
gr.functional.elements.df <- read.table(functional.elements.file.path, sep="\t", header = TRUE, stringsAsFactors = FALSE)
gr.functional.elements <- GRanges(gr.functional.elements.df$seqnames, IRanges(gr.functional.elements.df$start, gr.functional.elements.df$end))
gr.functional.elements$id <- gr.functional.elements.df$id

# nonCod.acc.elements.file.path <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join/join_filtered_elements_norm_nonCod.csv"
# nonCod.data <- read.table(nonCod.acc.elements.file.path, sep="\t", header = TRUE)
# nonCod.gr <- GRanges(seqnames = nonCod.data$bed_chr, IRanges(start = nonCod.data$bed_start, end = nonCod.data$bed_end))

# overlaps <- findOverlaps(query = nonCod.gr, subject = gr.functional.elements, type = "any")

# uso los indices de overlaps para indexar el data frame nonCod.data
# nonCod.data.functional <- nonCod.data[queryHits(overlaps), ]

# pego los datos de noCod.gr con los de functional.elements
nonCod.data.functional.count <- nrow(nonCod.data.functional)
# nonCod.data.functional$phastCons_chr <- rep("", nonCod.data.functional.count)
# nonCod.data.functional$phastCons_start <- rep(0, nonCod.data.functional.count)
# nonCod.data.functional$phastCons_end <- rep(0, nonCod.data.functional.count)
nonCod.data.functional$phastCons_id <- rep("", nonCod.data.functional.count)


for(i in c(1:nonCod.data.functional.count)) {
  subject.index <- subjectHits(overlaps)[i]
  phastCons.data <- gr.functional.elements[subject.index]
  # chr <- as.character(seqnames(phastCons.data))
  # start <- as.numeric(start(phastCons.data))
  # end <- as.numeric(end(phastCons.data))
  id <- as.character(phastCons.data$id)

  # nonCod.data.functional[i, "phastCons_chr"] <- chr
  # nonCod.data.functional[i, "phastCons_start"] <- start
  # nonCod.data.functional[i, "phastCons_end"] <- end
  nonCod.data.functional[i, "phastCons_id"] <- id
}


# cuento cuantos elementos acelerados hay por elemento funcional:
# elements.count.functional.element <- as.data.frame(table(nonCod.data.functional$phastCons_id))
# colnames(elements.count.functional.element) <- c("phastCons_id", "acc_elems_count")

# armo una tabla con datos de elemento funcional (phastCons)
functional.element.count <- nrow(elements.count.functional.element)
functional.element.data <- data.frame( id = numeric(functional.element.count),
                                       acc_elements_count = numeric(functional.element.count),
                                       functional_element_shift_count = numeric(functional.element.count),
                                       functional_element_gap_count = numeric(functional.element.count),
                                       acc_elements_shift_counts = character(functional.element.count), # lista de shifts por elementos
                                       acc_elements_gap_counts = character(functional.element.count), # lista de gaps por elementos
                                       stringsAsFactors = FALSE )


functional.element.data$id <- elements.count.functional.element$phastCons_id
functional.element.data$acc_elements_count <- elements.count.functional.element$acc_elems_count

#id <- 10045
get.functional.element.counts <- function(id, nonCod.data.functional) {
  
  nonCod <- nonCod.data.functional[nonCod.data.functional$phastCons_id == id, ]
  nonCod.gr <- GRanges(nonCod$bed_chr, IRanges(nonCod$bed_start, nonCod$bed_end))
  
  acc_elements_shift_counts <- paste(nonCod$shift_count, collapse = " ")
  acc_elements_gap_counts <- paste(nonCod$gap_count, collapse = " ")
  
  nonCod.gr.within <- findOverlaps(nonCod.gr, nonCod.gr, type="within")
  nonCod.gr.equal <- findOverlaps(nonCod.gr, nonCod.gr, type="equal")
  nonCod.gr.within.ne <- setdiff(nonCod.gr.within, nonCod.gr.equal)
  
  if(length(nonCod.gr.within.ne) > 0) {
  
    to.delete.indexes <- queryHits(nonCod.gr.within.ne)
    shifts.functional.element <- nonCod$shift_count[-to.delete.indexes]
    functional_element_shift_count <- sum(shifts.functional.element)
    gaps.functional.element <- nonCod$gap_count[-to.delete.indexes]
    functional_element_gap_count <- sum(gaps.functional.element)
    
  } else {
    functional_element_shift_count <- sum(nonCod$shift_count)
    functional_element_gap_count <- sum(nonCod$gap_count)
  }

  result <- list(
    acc_elements_shift_counts = acc_elements_shift_counts, 
    acc_elements_gap_counts = acc_elements_gap_counts, 
    functional_element_shift_count = functional_element_shift_count,
    functional_element_gap_count = functional_element_gap_count)
  
  return(result)
}

count.data <- sapply(functional.element.data$id, function(x) get.functional.element.counts(x, nonCod.data.functional))
colnames(count.data) <- functional.element.data$id
count.data <- t(count.data)
functional.element.data$acc_elements_shift_counts <- as.character(count.data[, "acc_elements_shift_counts"])
functional.element.data$acc_elements_gap_counts <- as.character(count.data[, "acc_elements_gap_counts"])
functional.element.data$functional_element_shift_count <- as.numeric(count.data[, "functional_element_shift_count"])
functional.element.data$functional_element_gap_count <- as.numeric(count.data[, "functional_element_gap_count"])
# agrego los datos de localizacion del phastcons:
functional.element.data <- merge(functional.element.data, gr.functional.elements.df,
                                 by= c("id"), all.x = TRUE)

############ filtro

functional.element.exclusions.shift.count <- function(i, functional.element.data) {
  
  acc.elements.gap.count <- functional.element.data$acc_elements_gap_counts[i]
  acc.elements.shift.count <- functional.element.data$acc_elements_shift_counts[i]
  
  gaps.count <- as.numeric(unlist(strsplit(acc.elements.gap.count, " ")))
  shifts.count <- as.numeric(unlist(strsplit(acc.elements.shift.count, " ")))
  
  index.no.gaps <- which(gaps.count == 0)
  print(index.no.gaps)
  
  if(length(index.no.gaps) > 0) {
    functional.element.shift.count.exclgt0gap  <- sum(shifts.count[index.no.gaps])
    acc.elements.count.exclgt0gap <- length(index.no.gaps)
  } else {
    functional.element.shift.count.exclgt0gap <- 0
    acc.elements.count.exclgt0gap <- 0
  }
  
  result <- list(functional.element.shift.count.exclgt0gap = functional.element.shift.count.exclgt0gap,
                 acc.elements.count.exclgt0gap = acc.elements.count.exclgt0gap)
  
  return(result)
  
}


filter.functional.element.one.acc.with.zero.gap <- function(functional.element.data){
  
  result <- functional.element.data
  
  functional.elements.count <- nrow(functional.element.data)
  
  correction.excl <- sapply(c(1:functional.elements.count), function(i) functional.element.exclusions.shift.count(i, functional.element.data))
  
  shifts.corrected <- unlist(correction.excl["functional.element.shift.count.exclgt0gap", ])
  
  count.corrected <- unlist(correction.excl["acc.elements.count.exclgt0gap", ])
  
  result$functional_element_shift_count_exclgt0gap <- shifts.corrected
  
  result$acc_elements_count_exclgt0gap <- count.corrected
 
  # conservo las regiones funcionales que tienen al menos un elemento acelerado sin gaps 
  result <- result[result$acc_elements_count_exclgt0gap > 0, ]
  
  return(result)
  
}

# functional.element.data.zerogap <- functional.element.data[functional.element.data$functional_element_gap_count == 0, ]

functional.element.data.oneacczerogap <- filter.functional.element.one.acc.with.zero.gap(functional.element.data)

#fitered.data <- functional.element.data.zerogap

fitered.data <- functional.element.data.oneacczerogap
  
##############################

fitered.data$functional.element.shift.count.exclgt0gap_rel <- fitered.data$functional_element_shift_count_exclgt0gap / fitered.data$width
fitered.data$acc.elements.count.exclgt0gap_rel <- fitered.data$acc_elements_count_exclgt0gap / fitered.data$width

plot(fitered.data$functional.element.shift.count.exclgt0gap_rel, fitered.data$acc.elements.count.exclgt0gap_rel)

hist(fitered.data$acc.elements.count.exclgt0gap_rel)

# guardo los datos de regiones funcionales (como resumen de los nonCod)
data.to.save <- merge(nonCod.data.functional, fitered.data, 
                      by.x=c("phastCons_id"), by.y = c("id"))

colnames(data.to.save) <- c( "phastCons_id", "acc_element_bed_chr", "acc_element_bed_start",
                             "acc_element_bed_end", "id", "acc_element_len",                          
                             "acc_element_size", "acc_element_gap_count", "acc_element_shift_count",
                             "hamming_ingroup_incl_gaps", "hamming_outgroup_incl_gaps", "hamming_in_out_incl_gaps",      
                             "hamming_ingroup_excl_gaps", "hamming_outgroup_excl_gaps", "hamming_in_out_excl_gaps",      
                             "acc_element_table_browser", "acc_element_shift_count_rel", "hamming_in_out_excl_gaps_scl",  
                             "shift_count_rel_scl", "norm2_scl", "coding",
                             "acc_elements_in_phastCons_count", "functional_element_shift_count", "functional_element_gap_count",  
                             "acc_elements_in_phastCons_shift_counts", "acc_elements_in_phastCons_gap_counts",
                             "phastCons_seqnames", "phastCons_start", "phastCons_end", "phastCons_width",                         
                             "phastCons_strand", "functional_elements_shift_rel",  "functional_elements_acc_elements_count_rel", 
                             "functional_element_shift_count_exclgt0gap", "acc_elements_count_exclgt0gap")


#out.file.name <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join/join_filtered_elements_norm_cod_nonCod_zerogap_functionalregions.csv"
out.file.name <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv"
write.table(data.to.save, out.file.name, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


colnames(data.to.save)
length(unique(data.to.save$phastCons_id))

# gr.functional.elements.df <- as.data.frame(gr.functional.elements)
# 
# elements.count.functional.element.data <- merge(elements.count.functional.element, gr.functional.elements.df, 
#                                                 by.x = c("phastCons_id"), by.y = c("id"), all.x = TRUE)
# 
# elements.count.functional.element.data$acc_elems_count_relative <- elements.count.functional.element.data$acc_elems_count / elements.count.functional.element.data$width

#pego la info del phastCons que contiene a cada elemento

#--------------------------------------------------------
#--------------------------------------------------------
# quiero ver si hay solapamiento parcial entre los elementos acelerados reportados: 
any.partial.overlap <- function(i, nonCod.gr) {
  gr.query <- nonCod.gr[i]
  overlaps.any <- findOverlaps(query = gr.query, subject = nonCod.gr, type = "any", minoverlap = 2)
  overlaps.eq  <- findOverlaps(query = gr.query, subject = nonCod.gr, type = "equal")
  overlaps.within <- findOverlaps(query = gr.query, subject = nonCod.gr, type = "within")
  overlaps.any.neq <- setdiff(overlaps.any, overlaps.eq)
  overlaps.any.neq.nw <- setdiff(overlaps.any.neq, overlaps.within)
  overlap.len <- length(overlaps.any.neq.nw)
  if(overlap.len > 0) {
    query.index <- queryHits(overlaps.any.neq.nw)
    subject.index <- subjectHits(overlaps.any.neq.nw)
    # reviso la inclusion inversa: subject en query
    overlap.subject.query <- findOverlaps(nonCod.gr[subject.index], gr.query[query.index], type = "within")
    is.inv.inclusion <- (length(overlap.subject.query) > 0)
    result <- ! is.inv.inclusion
  } else {
    result <- FALSE
  }
  return(result)
}

none.partial.overlps <- function() {
  nonCod.acc.elements.file.path <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join/join_filtered_elements_norm_nonCod.csv"
  nonCod.data <- read.table(nonCod.acc.elements.file.path, sep="\t", header = TRUE)
  nonCod.gr <- GRanges(seqnames = nonCod.data$bed_chr, IRanges(start = nonCod.data$bed_start, end = nonCod.data$bed_end))
  any.partial.overlap.item <- sapply(c(1:length(nonCod.gr)), function(i) any.partial.overlap(i, nonCod.gr))
  any.overlaps <- which(any.partial.overlap.item)
  result <- length(any.overlaps) == 0

  return(result)
}

# ningun overlap parcial:
no.partial.overlap <- none.partial.overlps()
no.partial.overlap

#--------------------------------------------------------
#--------------------------------------------------------


