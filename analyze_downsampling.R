library(tidyverse)
library(ggsci)
options(stringsAsFactors = FALSE)

DW_CODE <- c('1.00', '0.80', '0.60', '0.40', '0.20')
DW_RATE <- c(1.0, 0.8, 0.6, 0.4, 0.2)
LIBNAME <- c('control_1', 'control_2', 'control_3', 'cold_1', 'cold_2', 'cold_3')


# down-sampling was sampled from cleaned FASTQ, so I will estimate the size of un-cleaned FASTQ with downsampling rate
FQSIZE <- c(2604973 , 1596414, 3694206, 1764702, 4347049, 2604973)  # reads of the raw FASTQ


dfsz <- data.frame(n_read = NULL, n_gene = NULL, libname = NULL, dw_rate = NULL)

for (dw_code in DW_CODE) {
    fpath <- paste0('aaic_data/cs_cold_tagseq_dwsample/clean_tagseq_', dw_code, '/counts_1.0k/all.counts.gene.tsv.gz')
    x <- read.table(fpath, sep = '\t', header = TRUE)
    x <- x[, -1]
    x <- x[, c(1, 5, 6, 2, 3, 4)]
    colnames(x) <- LIBNAME
    
    dfsz <- rbind(dfsz, data.frame(n_read = round(FQSIZE * as.numeric(dw_code)),
                                   n_gene = colSums(x > 0), libname = LIBNAME, dw_rate = dw_code))
}



gp <- ggplot(dfsz, aes(x = n_read, y = n_gene, color = dw_rate)) +
        geom_point() + scale_color_npg() +
        xlab('number of reads') + ylab('number of genes')

gp


png(paste0('results/plots/downsampling.png'), 900, 600, res = 220)
print(gp)
dev.off()
 


