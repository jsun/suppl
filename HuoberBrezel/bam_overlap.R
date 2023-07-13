library(gplots)

lib_names <- c('20181109.A-TaeRS_1_Tae_RS1_1', '20181109.A-TaeRS_1_Tae_RS1_2', '20181109.A-TaeRS_1_Tae_RS1_3',
               '20181109.A-TaeRS_1_Tae_RS1_16', '20181109.A-TaeRS_1_Tae_RS1_17', '20181109.A-TaeRS_1_Tae_RS1_18',
               'TaeRS2728_g086', 'TaeRS2728_g160', 'TaeRS2728_g167',
               'TaeRS2728_g091', 'TaeRS2728_g165', 'TaeRS2728_g172')
hisat_prefix <- '/Users/jsun/Desktop/HB/data/tagseq/bam/'
eagle_prefix <- '/Users/jsun/Desktop/HB/data/tagseq/eaglerc/'


statsdf <- matrix(NA, ncol = 9, nrow = length(lib_names))
rownames(statsdf) <- lib_names

for (lib_name in lib_names) {
    hisat_reads <- read.table(paste0(hisat_prefix, lib_name, '.chrA.read_name.txt'), header = FALSE)
    eagle_reads <- read.table(paste0(eagle_prefix, lib_name, '/', lib_name, '.chrA.read_name.txt'), header = FALSE)
    statsdf[lib_name, 1] <- length(setdiff(hisat_reads[, 1], eagle_reads[, 1]))
    statsdf[lib_name, 2] <- length(intersect(hisat_reads[, 1], eagle_reads[, 1]))
    statsdf[lib_name, 3] <- length(setdiff(eagle_reads[, 1], hisat_reads[, 1]))
        
    hisat_reads <- read.table(paste0(hisat_prefix, lib_name, '.chrB.read_name.txt'), header = FALSE)
    eagle_reads <- read.table(paste0(eagle_prefix, lib_name, '/', lib_name, '.chrB.read_name.txt'), header = FALSE)
    statsdf[lib_name, 4] <- length(setdiff(hisat_reads[, 1], eagle_reads[, 1]))
    statsdf[lib_name, 5] <- length(intersect(hisat_reads[, 1], eagle_reads[, 1]))
    statsdf[lib_name, 6] <- length(setdiff(eagle_reads[, 1], hisat_reads[, 1]))
    
    hisat_reads <- read.table(paste0(hisat_prefix, lib_name, '.chrD.read_name.txt'), header = FALSE)
    eagle_reads <- read.table(paste0(eagle_prefix, lib_name, '/', lib_name, '.chrD.read_name.txt'), header = FALSE)
    statsdf[lib_name, 7] <- length(setdiff(hisat_reads[, 1], eagle_reads[, 1]))
    statsdf[lib_name, 8] <- length(intersect(hisat_reads[, 1], eagle_reads[, 1]))
    statsdf[lib_name, 9] <- length(setdiff(eagle_reads[, 1], hisat_reads[, 1]))
}


statsdf



