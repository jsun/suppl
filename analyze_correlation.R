library(tidyverse)
library(ggsci)
options(stringsAsFactors = FALSE)

LIBNAME <- c('control_1', 'control_2', 'control_3', 'cold_1', 'cold_2', 'cold_3')
LIBNAME <- c('replicate_1', 'replicate_2', 'replicate_3', 'replicate_4', 'replicate_5', 'replicate_6')


genemat2homeologmat <- function(mat) {
    abd <- read.table('aaic_data/homeolog.ABD.list', sep = '\t', header = FALSE)
    
    in_a <- abd[, 1] %in% rownames(mat)
    in_b <- abd[, 2] %in% rownames(mat)
    in_d <- abd[, 3] %in% rownames(mat)
    in_abd <- (in_a & in_b & in_d)
    
    abd <- abd[in_abd, ]
    
    mat.a <- mat[abd[, 1], ]
    mat.b <- mat[abd[, 2], ]
    mat.d <- mat[abd[, 3], ]
    
    list(A = mat.a, B = mat.b, D = mat.b)
}


calc_tpm <- function(mat) {
    tx2len <- read.table('aaic_data/cds.length.tsv', header = FALSE, sep = '\t')
    colnames(tx2len) <- c('tx', 'len')
    tx2len$gene <- sapply(strsplit(tx2len$tx, '\\.'), '[', 1)
    gene2lendf <- tx2len %>% select(gene, len) %>% group_by(gene) %>% summarise(len = max(len))
    gene2len <- gene2lendf$len
    names(gene2len) <- gene2lendf$gene
    
    matlen <- rep(NA, length = nrow(mat))
    names(matlen) <- rownames(mat)
    matlen[names(gene2len)] <- as.numeric(gene2len)
    
    matlen[is.na(matlen)] <- median(gene2len)
    
    cpk <- sweep(mat, 1, 1000 / matlen, '*')
    #tpm <- sweep(cpk, 2, 1e6 / colSums(cpk), '*')
    #tpm
    cpk
}



calc_corr <- function(matx, maty, labx = NULL, laby = NULL, cutx = 5, cuty = 10) {
    
    cormat <- matrix(NA, nrow = 6, ncol = 3)
    rownames(cormat) <- LIBNAME
    colnames(cormat) <- c('A', 'B', 'D')

    dfxy <- data.frame(x = NULL, y = NULL, libname = NULL, subgenome = NULL)
    for (i in LIBNAME) {
        for (g in c('A', 'B', 'D')) {
            dfxy_ig <- data.frame(x = matx[[g]][, i], y = maty[[g]][, i],
                                  libname = i, subgenome = g)
            dfxy_ig_high <- dfxy_ig[(dfxy_ig$x > cutx &  dfxy_ig$y > cuty), ]
            cormat[i, g] <- cor(dfxy_ig_high[, 1], dfxy_ig_high[, 2], method = 'pearson')
            dfxy <- rbind(dfxy, dfxy_ig)
        }
    }
    gp <- ggplot(dfxy, aes(x = x, y = y)) +
            geom_point() + facet_grid(libname ~ subgenome) +
            scale_x_log10() + scale_y_log10() + xlab(labx) + ylab(laby)
    
    dfxy <- dfxy  %>% filter(libname == LIBNAME[1])
    gp2 <- ggplot(dfxy, aes(x = x, y = y)) +
            geom_point() + facet_grid(. ~ subgenome) +
            scale_x_log10() + scale_y_log10() + xlab(labx) + ylab(laby)
    
    
    list(fig = gp, cor = cormat, fig_replicate = gp2)
}




## load read counts data



## CS 3'-end RNA-Seq / HISAT / IWGSC+1k
x <- read.table('aaic_data/cs_cold_tagseq_iwgscgff/counts_1.0k/all.counts.gene.tsv.gz', sep = '\t', header = TRUE)
tagseq_hisat <- as.matrix(x[, c(2, 6, 7, 3, 4, 5)])
colnames(tagseq_hisat) <- LIBNAME
rownames(tagseq_hisat) <- x[, 1]
tagseq_hisat_h <- genemat2homeologmat(tagseq_hisat)



## CS 3'-end RNA-Seq / HISAT / IWGSC (original)
x <- read.table('aaic_data/cs_cold_tagseq_iwgscgff/counts_iwgsc/all.counts.gene.tsv.gz', sep = '\t', header = TRUE)
tagseq_hisat0 <- as.matrix(x[, c(2, 6, 7, 3, 4, 5)])
colnames(tagseq_hisat0) <- LIBNAME
rownames(tagseq_hisat0) <- x[, 1]
tagseq_hisat_h0 <- genemat2homeologmat(tagseq_hisat0)



## CS paired-end RNA-Seq / HISAT / IWGSC
x <- read.table('aaic_data/cs_cold_fullseq/fullseq_hisat/counts_iwgsc/all.counts.gene.tsv.gz', sep = '\t', header = TRUE)
fullseq_hisat <- as.matrix(x[, -1])
colnames(fullseq_hisat) <- LIBNAME
rownames(fullseq_hisat) <- x[, 1]
fullseq_hisat_tpm <- calc_tpm(fullseq_hisat)
fullseq_hisat_h <- genemat2homeologmat(fullseq_hisat)
fullseq_hisat_tpm_h <- genemat2homeologmat(fullseq_hisat_tpm)



## CS paired-end RNA-Seq / EAGEL-RC / IWGSC
lbnms <- c('20181221.A-ZH_W2017_1_CS_2', '20181221.A-ZH_W2017_1_CS_3', '20181221.A-ZH_W2017_1_CS_4',
           '20181221.A-ZH_W2017_1_CS_cold_2', '20181221.A-ZH_W2017_1_CS_cold_3', '20181221.A-ZH_W2017_1_CS_cold_4')
fullseq_eagle <- NULL
for (lbnm in lbnms) {
    xa <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts/', lbnm, '.chrA.gene.counts.txt'), sep = '\t', header = TRUE)[, c(1, 7)]
    xb <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts/', lbnm, '.chrB.gene.counts.txt'), sep = '\t', header = TRUE)[, c(1, 7)]
    xd <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts/', lbnm, '.chrD.gene.counts.txt'), sep = '\t', header = TRUE)[, c(1, 7)]
    colnames(xa) <- colnames(xb) <- colnames(xd) <- c('gene', 'count')
    xabd <- rbind(xa, xb, xd)
    if (is.null(fullseq_eagle)) {
        fullseq_eagle <- as.matrix(xabd[, 2])
    } else {
        fullseq_eagle <- cbind(fullseq_eagle, as.matrix(xabd[, 2]))
    }
    rownames(fullseq_eagle) <- xabd[, 1]
}
colnames(fullseq_eagle) <- LIBNAME
fullseq_eagle_tpm <- calc_tpm(fullseq_eagle)
fullseq_eagle_h <- genemat2homeologmat(fullseq_eagle)
fullseq_eagle_tpm_h <- genemat2homeologmat(fullseq_eagle_tpm)





## check the order of genes in expression matrix are exactly same
table(rownames(tagseq_hisat_h0$A) == rownames(tagseq_hisat_h$A))
table(rownames(tagseq_hisat_h0$B) == rownames(tagseq_hisat_h$B))
table(rownames(tagseq_hisat_h0$D) == rownames(tagseq_hisat_h$D))
table(rownames(fullseq_eagle_tpm_h$A) == rownames(tagseq_hisat_h$A))
table(rownames(fullseq_eagle_tpm_h$B) == rownames(tagseq_hisat_h$B))
table(rownames(fullseq_eagle_tpm_h$D) == rownames(tagseq_hisat_h$D))
table(rownames(fullseq_hisat_tpm_h$A) == rownames(tagseq_hisat_h$A))
table(rownames(fullseq_hisat_tpm_h$B) == rownames(tagseq_hisat_h$B))
table(rownames(fullseq_hisat_tpm_h$D) == rownames(tagseq_hisat_h$D))




## calculate correlation


## 3'-end RNA-Seq (IWGSC+1k) vs paired-end RNA-Seq
tagseq_fullseqhisat <- calc_corr(tagseq_hisat_h, fullseq_hisat_tpm_h, "3'-end RNA-Seq (IWGSC+1k)", "paired-end RNA-Seq (HISAT)", 5, 10)
tagseq_fullseqeagle <- calc_corr(tagseq_hisat_h, fullseq_eagle_tpm_h, "3'-end RNA-Seq (IWGSC+1k)", "paired-end RNA-Seq (EAGLE-RC)", 5, 10)

## 3'-end RNA-Seq (IWGSC) vs paired-end RNA-Seq
tagseq0_fullseqhisat <- calc_corr(tagseq_hisat_h0, fullseq_hisat_tpm_h, "3'-end RNA-Seq (IWGSC)", "paired-end RNA-Seq (HISAT)", 5, 10)
tagseq0_fullseqeagle <- calc_corr(tagseq_hisat_h0, fullseq_eagle_tpm_h, "3'-end RNA-Seq (IWGSC)", "paired-end RNA-Seq (EAGLE-RC)", 5, 10)

## paired-end RNA-Seq HISAT vs paired-end RNA-Seq EAGLE-RC
fullseqhisat_fullseqeagle <- calc_corr(fullseq_hisat_tpm_h, fullseq_eagle_tpm_h, "paired-end RNA-Seq (HISAT)", "paired-end RNA-Seq (EAGLE-RC)", 10, 10)



## plot figures

png(paste0('results/plots/tagseqHISATIWGSC1.0k_fullseqHISAT_correlation.png'), 1500, 2300, res = 220)
print(tagseq_fullseqhisat$fig)
dev.off()
png(paste0('results/plots/tagseqHISATIWGSC1.0k_fullseqEAGLERC_correlation.png'), 1500, 2300, res = 220)
print(tagseq_fullseqeagle$fig)
dev.off()
png(paste0('results/plots/tagseqHISATIWGSC1.0k_fullseqEAGLERC_correlation_rep1.png'), 1500, 600, res = 220)
print(tagseq_fullseqeagle$fig_replicate)
dev.off()


png(paste0('results/plots/tagseqHISATIWGSC_fullseqHISAT_correlation.png'), 1500, 2300, res = 220)
print(tagseq0_fullseqhisat$fig)
dev.off()
png(paste0('results/plots/tagseqHISATIWGSC_fullseqEAGLERC_correlation.png'), 1500, 2300, res = 220)
print(tagseq0_fullseqeagle$fig)
dev.off()
png(paste0('results/plots/tagseqHISATIWGSC_fullseqEAGLERC_correlation_rep1.png'), 1500, 600, res = 220)
print(tagseq0_fullseqeagle$fig_replicate)
dev.off()


png(paste0('results/plots/fullseqHISAT_fullseqEAGLERC_correlation.png'), 1500, 2300, res = 220)
print(fullseqhisat_fullseqeagle$fig)
dev.off()


tagseq_fullseqhisat$cor
tagseq_fullseqeagle$cor
fullseqhisat_fullseqeagle$cor
tagseq0_fullseqhisat$cor
tagseq0_fullseqeagle$cor








