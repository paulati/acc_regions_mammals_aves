base_path <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join"
setwd(base_path)

file_name <- "join_filtered_elements.csv"

data <- read.table(file_name, sep="\t", stringsAsFactors = FALSE)
colnames(data) <- c("id", "len", "size", "gap_count", "shift_count", 
                    "hamming_ingroup_incl_gaps", "hamming_outgroup_incl_gaps", "hamming_in_out_incl_gaps", 
                    "hamming_ingroup_excl_gaps", "hamming_outgroup_excl_gaps", "hamming_in_out_excl_gaps",
                    "table_browser",  "bed_chr", "bed_start", "bed_end")

#cuantos registros hay en total:
nrow(data)

# size == len + gaps ?
# bad.indexes <- which(data$size != data$len + data$gap_count)
# length(bad.indexes)
# bad.data <- data[bad.indexes, ]
# ver que pasa con esats secuencias. bad.data deberia ser vacío

# diff.indexes <- which(data$hamming_ingroup_incl_gaps != data$hamming_ingroup_excl_gaps)

# diff.indexes <- which(data$hamming_outgroup_incl_gaps != data$hamming_outgroup_excl_gaps)

# # cuantos difieren en las distancia de hamming con y sin gaps
# which(data$hamming_in_out_incl_gaps < data$hamming_in_out_excl_gaps)
# tmp1 <- data[, ]
# nrow(tmp1)

# considero shifts solo cuando el cambio no es un gap
# y excluyo de las secuencias sobre las que calculo hamming unicamente las posiciones que son shift y gap
# esa es la razon por la que hamming con y sin gaps da lo mismo en las columnas del dataframe (demostrarlo)


# cuantos coinciden en las distancia de hamming con y sin gaps

# tmp2 <- data$hamming_in_out_excl_gaps == data$hamming_in_out_incl_gaps
# which(!tmp2)
# cuantos tienen distancia de hamming con gaps > distancia de hamming sin gaps

# cuantos tienen distancia de hamming con gaps < distancia de hamming sin gaps


data.mammals <- data[data$hamming_in_out_excl_gaps > 0, ]
data.vertebrates <- data[data$hamming_in_out_excl_gaps < 0, ]

# data.mammals <- data
########################################################################################

data.mammals$shift_count_rel <- data.mammals$shift_count / data.mammals$size

# escalo a 0 1
max.distance <- max(data.mammals$hamming_in_out_excl_gaps)
min.distance <- min(data.mammals$hamming_in_out_excl_gaps)
data.mammals$hamming_in_out_excl_gaps_scl <- (data.mammals$hamming_in_out_excl_gaps - rep(min.distance, nrow(data.mammals)))/ (max.distance - min.distance)

max.shift.count <- max(data.mammals$shift_count_rel)
min.shift.count <- min(data.mammals$shift_count_rel)
data.mammals$shift_count_rel_scl <- (data.mammals$shift_count_rel - rep(min.shift.count, nrow(data.mammals))) / (max.shift.count - min.shift.count)

plot(data.mammals$shift_count_rel_scl, data.mammals$hamming_in_out_scl)
hist(data.mammals$shift_count_rel_scl)
plot(data.mammals$shift_count_rel_scl)
hist(data.mammals$hamming_in_out_scl)

#max.distance <- data.mammals[(data.mammals$distance.rel > 0.8), ]
#max.distance.shift <- max.distance[max.distance$shift.rel > 0.6, ]

# calculo el modulo del vector formado por distance.rel y shift.rel
data.mammals$norm2_scl <- apply(data.mammals[, c("hamming_in_out_excl_gaps_scl", "shift_count_rel_scl")], 1, function(x) norm(x, type="2"))


data.mammals.sort <- data.mammals[order(-data.mammals$norm2_scl), ]

write.table(data.mammals.sort, "join_filtered_elements_norm.csv", sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#hist(data.mammals$norm2_scl)

# plot(data.mammals$norm2, data.mammals$V2)  


#agregar algun dato que tenga en cuenta el largo del cromosoma

table(data.mammals.sort$bed_chr)

########################################################################################

data.vertebrates$shift_count_rel <- data.vertebrates$shift_count / data.vertebrates$size

# escalo a 0 1
max.distance <- max(data.vertebrates$hamming_in_out_excl_gaps)
min.distance <- min(data.vertebrates$hamming_in_out_excl_gaps)
data.vertebrates$hamming_in_out_excl_gaps_scl <- 1 - (data.vertebrates$hamming_in_out_excl_gaps - rep(min.distance, nrow(data.vertebrates)))/ (max.distance - min.distance)

max.shift.count <- max(data.vertebrates$shift_count_rel)
min.shift.count <- min(data.vertebrates$shift_count_rel)
data.vertebrates$shift_count_rel_scl <- (data.vertebrates$shift_count_rel - rep(min.shift.count, nrow(data.vertebrates))) / (max.shift.count - min.shift.count)

plot(data.vertebrates$shift_count_rel_scl, data.vertebrates$hamming_in_out_scl)
hist(data.vertebrates$shift_count_rel_scl)
plot(data.vertebrates$shift_count_rel_scl)
hist(data.vertebrates$hamming_in_out_scl)

#max.distance <- data.mammals[(data.mammals$distance.rel > 0.8), ]
#max.distance.shift <- max.distance[max.distance$shift.rel > 0.6, ]

# calculo el modulo del vector formado por distance.rel y shift.rel
data.vertebrates$norm2_scl <- apply(data.vertebrates[, c("hamming_in_out_excl_gaps_scl", "shift_count_rel_scl")], 1, function(x) norm(x, type="2"))

data.vertebrates.sort <- data.vertebrates[order(-data.vertebrates$norm2_scl), ]
