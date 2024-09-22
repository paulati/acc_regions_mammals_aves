base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'

parse_data <- function(file_path) {
  file_name <- basename(file_path)
  parts1 <- unlist(strsplit(file_name, split =".", fixed = TRUE))
  if(length(parts1) > 0) {
    parts2 <- unlist(strsplit(parts1[1], split ="_", fixed = TRUE))
    if(length(parts2) >= 3) {
      result <- list(chr= parts2[1], length=parts2[3])
    } else {
      result <- list(chr= NA, length=NA)
    }
  } else {
    result <- list(chr= NA, length=NA)
  }
  return(result)
}

source_base_path <- file.path(base_path, "mammals/acceleration/data/obs_phyloP/sarcopterygii_mammals")
out_base_path <- file.path(base_path, "mammals/acceleration/data/obs_phyloP/sarcopterygii_mammals/obsPhyloP_elements_bed")

files <- list.files(source_base_path, pattern = ".csv", full.names = TRUE)

file <- files[1]

for (file in files) {

  data <- read.table(file, sep="\t", header = TRUE )
  data_bed <- data[, c("chr", "start", "end")]
  
  atts <- parse_data(file)
  
  out_file_name <- paste0(atts$chr, "_obsPhyloP_", atts$length, ".bed")
  out_file_path <- file.path(out_base_path, out_file_name)
  
  write.table(data_bed, out_file_path, sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

}
