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
    
  result.df <- data.frame(chr = character(),
                       start = numeric(),
                       end = numeric(),
                       stringsAsFactors = FALSE)
  
  for(chr in c(1:22)) {
  
    file.name <- paste0("chr", chr, "_intersection_sarcopterygii_mammals.bed")    
    file_path <- file.path(base.path, file.name)
    data <- read.table(file_path, sep=" ", header=FALSE, stringsAsFactors = FALSE)  
    colnames(data) <- c("chr", "start", "end")
    result.df <- rbind(result.df, data)
  
  }
  
  result <- GRanges(result.df$chr, IRanges(start = result.df$star, end = result.df$end))
  return(result)
  
}

build.functional.elements <- function(data.base.path) {

  #data.base.paths <- file.path(base_path, 'mammals/conservation/data/output/intersection_modneu_phastCons100way_sarcopterygii_mammals')
      #"/paula/2018/acc_regions/resultados/intersection_modneu_phastCons100way_sarcopterygii_mammals"
  gr <- load.phastCons(data.base.path)

  distance.threshold <- 20

  length.threshold <- 100

  gr.util <- collapse.neighbors(gr, distance.threshold, length.threshold)

  gr.util$id <- c(1:length(gr.util))
  
  gr.util.df <- as.data.frame(gr.util)
  
  out.file.path <- file.path(base_path, "mammals/data/results", 
  "intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed")
# "/paula/2018/acc_regions/resultados/intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed"
  write.table(gr.util.df, out.file.path, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)  
  
  return(gr.util)
  
}

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

get.acc.functional.element.data <- function(acc.elements.data.functional,
                                            gr.functional.elements.df) {
    
    # count hoy many acc elements are associated with each phastcons:
    elements.count.functional.element <- as.data.frame(
        table(acc.elements.data.functional$phastCons_id))
    colnames(elements.count.functional.element) <- c("phastCons_id", 
                                                     "acc_elems_count")
    
    # build a table with data associated with each functional element (phastCons)
    functional.element.count <- nrow(elements.count.functional.element)
    functional.element.data <- data.frame( id = numeric(functional.element.count),
                                           acc_elements_count = numeric(functional.element.count),
                                           functional_element_shift_count = numeric(functional.element.count),
                                           functional_element_gap_count = numeric(functional.element.count),
                                           acc_elements_shift_counts = character(functional.element.count), # lista de shifts por elementos
                                           acc_elements_gap_counts = character(functional.element.count), # lista de gaps por elementos
                                           stringsAsFactors = FALSE )
    
    functional.element.data$id <- elements.count.functional.element$phastCons_id
    functional.element.data$acc_elements_count <- 
        elements.count.functional.element$acc_elems_count
    
    count.data <- sapply(functional.element.data$id,
                         function(x) get.functional.element.counts(x, 
                                                                   acc.elements.data.functional))
    colnames(count.data) <- functional.element.data$id
    count.data <- t(count.data)
    functional.element.data$acc_elements_shift_counts <- 
        as.character(count.data[, "acc_elements_shift_counts"])
    functional.element.data$acc_elements_gap_counts <- 
        as.character(count.data[, "acc_elements_gap_counts"])
    functional.element.data$functional_element_shift_count <- 
        as.numeric(count.data[, "functional_element_shift_count"])
    functional.element.data$functional_element_gap_count <- 
        as.numeric(count.data[, "functional_element_gap_count"])
    # agrego los datos de localizacion del phastcons:
    functional.element.data <- merge(gr.functional.elements.df,
                                     functional.element.data, 
                                     by= c("id"), all.x = TRUE)    
    
    ##############elements.count.functional.element.data$acc_elems_count_relative <- elements.count.functional.element.data$acc_elems_count / elements.count.functional.element.data$width
    
    return(functional.element.data)
    
}

main <- function() {
    
    functional.elements.file.path <- file.path(base_path, "mammals/data/results", 
        "intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed")
    
    if(! file.exists(functional.elements.file.path)) {
        gr.functional.elements <- build.functional.elements(file.path(base_path, 
            "mammals/conservation/data/output/intersection_modneu_phastCons100way_sarcopterygii_mammals"))
    }
    
    gr.functional.elements.df <- read.table(functional.elements.file.path, 
                                            sep="\t", header = TRUE, 
                                            stringsAsFactors = FALSE)
    gr.functional.elements <- GRanges(gr.functional.elements.df$seqnames, 
                                      IRanges(gr.functional.elements.df$start, 
                                              gr.functional.elements.df$end))
    gr.functional.elements$id <- gr.functional.elements.df$id
    
    acc.elements.file.path <- file.path(base_path, 
        'mammals/scoring/data/join_filtered_elements_norm.csv')
    acc.elements.data <- read.table(acc.elements.file.path, sep="\t", 
                                    header = TRUE)
    acc.elements.gr <- GRanges(seqnames = acc.elements.data$bed_chr, 
                               IRanges(start = acc.elements.data$bed_start, 
                                       end = acc.elements.data$bed_end))
    
    overlaps <- findOverlaps(query = acc.elements.gr, 
                             subject = gr.functional.elements, type = "any")
    
    # length(overlaps)
    
    # using overlaps queryHits to index acc.elements.data
    acc.elements.data.functional <- acc.elements.data[queryHits(overlaps), ]
    phastcons_ids <- gr.functional.elements[subjectHits(overlaps)]@
        elementMetadata$id
    # add associated phastcons id to each acc element
    acc.elements.data.functional$phastCons_id <- phastcons_ids
    
    functional.element.data <- get.acc.functional.element.data(
        acc.elements.data.functional, gr.functional.elements.df)
    
    colnames(functional.element.data) <- c("phastCons_id", "phastCons_seqnames", 
                                           "phastCons_start", "phastCons_end",
                                           "phastCons_width", "phastCons_strand", 
                                           "acc_elements_count", 
                                           "functional_element_shift_count",
                                           "functional_element_gap_count",
                                           "acc_elements_shift_counts",
                                           "acc_elements_gap_counts" )
    
    ############ filtro
    
    # functional.element.data.zerogap <- functional.element.data[functional.element.data$functional_element_gap_count == 0, ]
    
    functional.element.data.oneacczerogap <- 
        filter.functional.element.one.acc.with.zero.gap(functional.element.data)
    
    #fitered.data <- functional.element.data.zerogap
    
    filtered.data <- functional.element.data.oneacczerogap
    
    ##############################
    
    filtered.data$functional.element.shift.count.exclgt0gap_rel <- 
        filtered.data$functional_element_shift_count_exclgt0gap / filtered.data$phastCons_width
    filtered.data$acc.elements.count.exclgt0gap_rel <- 
        filtered.data$acc_elements_count_exclgt0gap / filtered.data$phastCons_width
    
    #plot(fitered.data$functional.element.shift.count.exclgt0gap_rel, fitered.data$acc.elements.count.exclgt0gap_rel)
    
    #hist(fitered.data$acc.elements.count.exclgt0gap_rel)
    
    colnames(acc.elements.data.functional) <- c("id", "acc_element_len", 
                                                "acc_element_size",                        
                                                "acc_element_gap_count", 
                                                "acc_element_shift_count",
                                                "hamming_ingroup_incl_gaps",   
                                                "hamming_outgroup_incl_gaps",
                                                "hamming_in_out_incl_gaps",
                                                "hamming_ingroup_excl_gaps",
                                                "hamming_outgroup_excl_gaps",
                                                "hamming_in_out_excl_gaps",
                                                "acc_element_table_browser",
                                                "acc_element_bed_chr", 
                                                "acc_element_bed_start", 
                                                "acc_element_bed_end", 
                                                "acc_element_shift_count_rel", 
                                                "hamming_in_out_excl_gaps_scl",
                                                "shift_count_rel_scl",
                                                "norm2_scl", 
                                                "phastCons_id")
    
    
    data.to.save <- merge(acc.elements.data.functional, filtered.data, 
                          by=c("phastCons_id"))
    
    #out.file.name <- file.path(base_path, "/scoring/data/results/201901/filter/join/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv"
    out.file.name <- file.path(base_path, "mammals/data/output/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv")
    write.table(data.to.save, out.file.name, col.names = TRUE, row.names = FALSE, 
                quote = FALSE, sep = "\t")    
    
    
}




# test <- function() {
# 
# # test:
# 
# file_1 <- file.path(base_path, "mammals/data/output/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv")
# file_2 <- file.path(base_path, "mammals/data/output/mammals_acc_func_elements_all.txt")
# 
# data_1 <- read.table(file_1, sep = '\t', header = TRUE)
# data_2 <- read.table(file_2, sep = '\t', header = TRUE)
# 
# library(stringr)
# colnames(data_1) <- str_replace_all(colnames(data_1), pattern = "\\.", "_")
# colnames(data_2)
# 
# merge_data <- merge(data_1, data_2, by = c('id', 'phastCons_id',
#                            'acc_element_len', 'acc_element_size',
#                            'acc_element_gap_count', 'acc_element_shift_count',
#                            'phastCons_seqnames', 'phastCons_start',
#                            'phastCons_end', 'phastCons_width', 'phastCons_strand'),
#                     all.x = TRUE)
# 
# 
# }

