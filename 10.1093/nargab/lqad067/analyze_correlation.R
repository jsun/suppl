library(tidyverse)
library(ggsci)
library(ggExtra)
options(stringsAsFactors = FALSE)

# set TRUE or FALSE to switch the counts data of paired-end RNA-Seq
# TRUE: counts counted with IWGSC+1k; FLASE: with IWGSC
FULLSEQ_IWGSC1k <- TRUE


LIBNAME <- c('#1', '#2', '#3', '#4', '#5', '#6')

genemat2homeologmat <- function(mat) {
    abd <- read.table('data/homeolog.ABD.list', sep = '\t', header = FALSE)
    
    in_a <- abd[, 1] %in% rownames(mat)
    in_b <- abd[, 2] %in% rownames(mat)
    in_d <- abd[, 3] %in% rownames(mat)
    in_abd <- (in_a & in_b & in_d)
    
    abd <- abd[in_abd, ]
    
    mat.a <- mat[abd[, 1], ]
    mat.b <- mat[abd[, 2], ]
    mat.d <- mat[abd[, 3], ]
    
    list(A = mat.a, B = mat.b, D = mat.d)
}

get_gene_length <- function() {
    tx2len <- read.table('aaic_data/cds.length.tsv', header = FALSE, sep = '\t')
    colnames(tx2len) <- c('tx', 'len')
    tx2len$gene <- sapply(strsplit(tx2len$tx, '\\.'), '[', 1)
    gene2lendf <- tx2len %>% select(gene, len) %>% group_by(gene) %>% summarise(len = max(len))
    gene2len <- gene2lendf$len
    names(gene2len) <- gene2lendf$gene
    gene2len
}


calc_tpm <- function(mat) {
    gene2len <- get_gene_length()
    matlen <- rep(NA, length = nrow(mat))
    names(matlen) <- rownames(mat)
    matlen[names(gene2len)] <- as.numeric(gene2len)
    
    matlen[is.na(matlen)] <- median(gene2len)
    
    cpk <- sweep(mat, 1, 1000 / matlen, '*')
    tpm <- sweep(cpk, 2, 1e6 / colSums(cpk), '*')
    
    tpm
}

calc_cpm <- function(mat) {
    cpm <- sweep(mat, 2, 1e6 / colSums(mat), '*')
    cpm
}


calc_corr <- function(matx, maty, labx = NULL, laby = NULL, cutx = 5, cuty = 10) {
    outlier_upper <- 'TraesCS4B02G159100'
    outlier_lower <- 'TraesCS3B02G187700'
    
    cormat2 <- cormat <- rsqmat <- matrix(NA, nrow = 6, ncol = 3)
    rownames(cormat2) <- rownames(cormat) <- rownames(rsqmat) <- LIBNAME
    colnames(cormat2) <- colnames(cormat) <- colnames(rsqmat) <- c('A', 'B', 'D')

    dfxy <- data.frame(x = NULL, y = NULL, libname = NULL, subgenome = NULL)
    for (i in LIBNAME) {
        for (g in c('A', 'B', 'D')) {
            dfxy_ig <- data.frame(x = matx[[g]][, i], y = maty[[g]][, i],
                                  libname = i, subgenome = g)
            dfxy <- rbind(dfxy, dfxy_ig)
            dfxy_ig_high <- dfxy_ig[(dfxy_ig$x > cutx &  dfxy_ig$y > cuty), ]
            # corr (raw counts)
            cormat2[i, g] <- cor(dfxy_ig_high$x, dfxy_ig_high$y, method = 'pearson')
            
            dfxy_ig_high$x <- log10(dfxy_ig_high$x + 1)
            dfxy_ig_high$y <- log10(dfxy_ig_high$y + 1)
            # corr
            cormat[i, g] <- cor(dfxy_ig_high$x, dfxy_ig_high$y, method = 'pearson')
            # R^2
            xylm <- lm(y ~ x, data = dfxy_ig_high)
            rsqmat[i, g] <- summary(xylm)$r.squared
            
        }
    }
    dfxy$x <- log10(dfxy$x + 1)
    dfxy$y <- log10(dfxy$y + 1)
    bk <- c(0, 1, 2, 3, 4)
    lb <- c('0', '10', '100', '1000', '10000')
    lm <- c(-0.1, 4.1)
    gp <- ggplot(dfxy, aes(x = x, y = y)) +
            scale_x_continuous(breaks = bk, labels =lb, limits = lm) +
            scale_y_continuous(breaks = bk, labels =lb, limits = lm) +
            geom_point() + facet_grid(libname ~ subgenome) +
            xlab(labx) + ylab(laby) +
            coord_fixed()
    
    dfxy <- dfxy  %>% filter(libname == LIBNAME[1])
    gp2 <- ggplot(dfxy, aes(x = x, y = y)) +
            scale_x_continuous(breaks = bk, labels =lb, limits = lm) +
            scale_y_continuous(breaks = bk, labels =lb, limits = lm) +
            geom_point() + facet_grid(. ~ subgenome) +
            xlab(labx) + ylab(laby) + coord_fixed()
    
    
    list(fig = gp, cor = cormat, fig_replicate = gp2, rsq = rsqmat, cor2 = cormat)
}




## load read counts data
stop()


## CS 3'-end RNA-Seq / HISAT / IWGSC+1k
x <- read.table('data/tagseq/counts/cs.counts.gene.ext_1.0k.tsv.gz', sep = '\t', header = TRUE)
tagseq_hisat <- as.matrix(x[, -c(1:6)][, c(4, 5, 6, 1, 2, 3)])
colnames(tagseq_hisat) <- LIBNAME
rownames(tagseq_hisat) <- x[, 1]
tagseq_hisat <- calc_cpm(tagseq_hisat)
tagseq_hisat_h <- genemat2homeologmat(tagseq_hisat)



## CS 3'-end RNA-Seq / HISAT / IWGSC (original)
x <- read.table('data/tagseq/counts/cs.counts.gene.iwgsc.tsv.gz', sep = '\t', header = TRUE)
tagseq_hisat0 <- as.matrix(x[, -c(1:6)][, c(4, 5, 6, 1, 2, 3)])
colnames(tagseq_hisat0) <- LIBNAME
rownames(tagseq_hisat0) <- x[, 1]
tagseq_hisat0 <- calc_cpm(tagseq_hisat0)
tagseq_hisat_h0 <- genemat2homeologmat(tagseq_hisat0)



## CS paired-end RNA-Seq / HISAT / IWGSC
ftag <- ifelse (FULLSEQ_IWGSC1k, 'ext_1.0k', 'iwgsc')
x <- read.table(paste0('data/fullseq/counts/counts.gene.', ftag, '.tsv.gz'), sep = '\t', header = TRUE)
fullseq_hisat <- as.matrix(x[, -c(1:6)])
colnames(fullseq_hisat) <- LIBNAME
rownames(fullseq_hisat) <- x[, 1]
fullseq_hisat <- calc_tpm(fullseq_hisat)
fullseq_hisat_h <- genemat2homeologmat(fullseq_hisat)



## CS paired-end RNA-Seq / EAGEL-RC / IWGSC
## ftag <- ifelse (FULLSEQ_IWGSC1k, 'homeolog_counts_1.0k', 'homeolog_counts')
lbnms <- c('20181221.A-ZH_W2017_1_CS_2', '20181221.A-ZH_W2017_1_CS_3', '20181221.A-ZH_W2017_1_CS-4',
           '20181221.A-ZH_W2017_1_CS_cold_2', '20181221.A-ZH_W2017_1_CS_cold_3', '20181221.A-ZH_W2017_1_CS_cold_4')
fullseq_eagle <- NULL
for (lbnm in lbnms) {
    xa <- read.table(paste0('data/fullseq/countseaglrc/', lbnm, '.chrA.counts.homeolog.txt.gz'), sep = '\t', header = TRUE)[, c(1, 7)]
    xb <- read.table(paste0('data/fullseq/countseaglrc/', lbnm, '.chrB.counts.homeolog.txt.gz'), sep = '\t', header = TRUE)[, c(1, 7)]
    xd <- read.table(paste0('data/fullseq/countseaglrc/', lbnm, '.chrD.counts.homeolog.txt.gz'), sep = '\t', header = TRUE)[, c(1, 7)]
    colnames(xa) <- colnames(xb) <- colnames(xd) <- c('gene', 'count')
    xabd <- rbind(xa, xb, xd)
    fullseq_eagle <- cbind(fullseq_eagle, as.matrix(xabd[, 2]))
    rownames(fullseq_eagle) <- xabd[, 1]
}
colnames(fullseq_eagle) <- LIBNAME
fullseq_eagle <- calc_tpm(fullseq_eagle)
fullseq_eagle_h <- genemat2homeologmat(fullseq_eagle)





## check the order of genes in expression matrix are exactly same
chk <- c(
all(rownames(tagseq_hisat_h0$A) == rownames(tagseq_hisat_h$A)),
all(rownames(tagseq_hisat_h0$B) == rownames(tagseq_hisat_h$B)),
all(rownames(tagseq_hisat_h0$D) == rownames(tagseq_hisat_h$D)),
all(rownames(fullseq_eagle_h$A) == rownames(tagseq_hisat_h$A)),
all(rownames(fullseq_eagle_h$B) == rownames(tagseq_hisat_h$B)),
all(rownames(fullseq_eagle_h$D) == rownames(tagseq_hisat_h$D)),
all(rownames(fullseq_hisat_h$A) == rownames(tagseq_hisat_h$A)),
all(rownames(fullseq_hisat_h$B) == rownames(tagseq_hisat_h$B)),
all(rownames(fullseq_hisat_h$D) == rownames(tagseq_hisat_h$D))
)
if (!all(chk)) {stop('check gene orders!')}



## calculate correlation


## 3'-end RNA-Seq (IWGSC+1k) vs paired-end RNA-Seq
cutoff_tagseq <- 1
cutoff_fullseq <- 4
tagseq_fullseqhisat <- calc_corr(tagseq_hisat_h, fullseq_hisat_h, "Lasy-Seq (HISAT2)", "paired-end RNA-Seq (HISAT2)", cutoff_tagseq, cutoff_fullseq)
tagseq_fullseqeagle <- calc_corr(tagseq_hisat_h, fullseq_eagle_h, "Lasy-Seq (HISAT2)", "paired-end RNA-Seq (EAGLE-RC)", cutoff_tagseq, cutoff_fullseq)

## paired-end RNA-Seq HISAT vs paired-end RNA-Seq EAGLE-RC
fullseqhisat_fullseqeagle <- calc_corr(fullseq_hisat_h, fullseq_eagle_h, "paired-end RNA-Seq (HISAT)", "paired-end RNA-Seq (EAGLE-RC)", cutoff_fullseq, cutoff_fullseq)


## plot figures

png(paste0('results/plots/tagseqHISATIWGSC1.0k_fullseqHISAT_correlation.png'), 1500, 2500, res = 220)
print(tagseq_fullseqhisat$fig)
dev.off()
png(paste0('results/plots/tagseqHISATIWGSC1.0k_fullseqEAGLERC_correlation.png'), 1500, 2500, res = 220)
print(tagseq_fullseqeagle$fig)
dev.off()
png(paste0('results/plots/tagseqHISATIWGSC1.0k_fullseqEAGLERC_correlation_rep1.png'), 1800, 800, res = 220)
print(tagseq_fullseqeagle$fig_replicate)
dev.off()



png(paste0('results/plots/fullseqHISAT_fullseqEAGLERC_correlation.png'), 1500, 2500, res = 220)
print(fullseqhisat_fullseqeagle$fig)
dev.off()


tagseq_fullseqhisat$cor
tagseq_fullseqeagle$cor
fullseqhisat_fullseqeagle$cor




## normalized counts differences between 3'-end HISAT and paired-end RNA-Seq EAGLERC 
## are not correaltion with gene length

gene2len <- get_gene_length()
cormat <- rep(NA, length = 6)
dfxy <- data.frame(x = NULL, y = NULL, libname = NULL)
for (i in LIBNAME) {
    sn <- intersect(rownames(tagseq_hisat), rownames(fullseq_eagle))
    x <- tagseq_hisat[sn, i]
    y <- fullseq_eagle[sn, i]
    df_ <- data.frame(diff = y - x, gene_length = gene2len[sn], libname = i)
    df_ <- df_[(x > 0 & y > 0), ]
    dfxy <- rbind(dfxy, df_)
    
    png(paste0('results/plots/countdifference_genelength_', i, '.png'), 800, 800, res = 220)
    df_$diff <- log10(df_$diff)
    df_$gene_length <- log10(df_$gene_length)
    g <- ggplot(df_, aes(x = diff, y = gene_length)) + geom_point() +
         xlab('log10(difference)') + ylab('log10(gene length)')
    g <- ggMarginal(g, type = "histogram", margins = "both", size = 5)
    print(g)
    dev.off()
}


gp <- ggplot(dfxy, aes(x = diff, y = gene_length)) +
            geom_point() + facet_wrap(~libname, ncol = 3) +
            scale_x_log10() + scale_y_log10() + 
            xlab("paired-end RNA-Seq (EAGLE-RC) - Lasy-Seq (HISAT)") + ylab('gene length')
    
png(paste0('results/plots/countdifference_genelength.png'), 2300, 1500, res = 220)
print(gp)
dev.off()






