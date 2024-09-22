cran_dependencies <- c('devtools', 'factoextra', 'ggplot2', 'aws.s3', 'R.utils',
                       'ape', 'BiocManager', 'rstudioapi')

bioc_dependencies <- c('GenomicRanges')


for(cran_dependence in cran_dependencies) {
    if (!require(cran_dependence, quietly = TRUE))
        install.packages(cran_dependence)
}

for(bioc_dependence in bioc_dependencies) {
    if (!require(bioc_dependence, quietly = TRUE))
        BiocManager::install(bioc_dependence)
}


if (!require(rphast, quietly = TRUE))
    devtools::install_github("CshlSiepelLab/RPHAST")



