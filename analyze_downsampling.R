library(tidyverse)
library(ggsci)
options(stringsAsFactors = FALSE)

DW_CODE <- c('1.00', '0.80', '0.60', '0.40', '0.20')
DW_RATE <- c(1.0, 0.8, 0.6, 0.4, 0.2)
LIBNAME <- c('control_1', 'control_2', 'control_3', 'cold_1', 'cold_2', 'cold_3')
NGFF <- 107903

# down-sampling was sampled from cleaned FASTQ, so I will estimate the size of un-cleaned FASTQ with downsampling rate
FQSIZE <- c(2604973 , 1596414, 3694206, 1764702, 4347049, 2604973)  # reads of the raw FASTQ




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

    list(A = mat.a, B = mat.b, D = mat.d)
}




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
 
fit <- nls(n_gene ~ a * log10(n_read) + b, start = c(a = 1, b = 1), data = dfsz)
a <- summary(fit)$coefficients[2]
b <- summary(fit)$coefficients[1]
newx <- seq(min(dfsz$n_read), max(dfsz$n_read), 100)
dffit <- data.frame(n_read = newx, n_gene = predict(fit, newdata = data.frame(n_read = newx)))
print(predict(fit, newdata = data.frame(n_read = c(1e6, 2e6, 3e6, 4e6))))
print(predict(fit, newdata = data.frame(n_read = c(1e6, 2e6, 3e6, 4e6)))/NGFF)


gp <- ggplot() +
        geom_point(aes(x = n_read, y = n_gene, color = dw_rate), data = dfsz) +
        geom_line(aes(x = n_read, y = n_gene), data = dffit) +
        scale_color_npg() +
        xlab('number of reads') + ylab('number of genes') +
        theme(legend.position = 'top') +
        labs(color = 'rate')

png(paste0('results/plots/downsampling.png'), 900, 840, res = 220)
print(gp)
dev.off()








dfsz <- data.frame(n_read = NULL, n_gene = NULL, libname = NULL, dw_rate = NULL)
for (dw_code in DW_CODE) {
    fpath <- paste0('aaic_data/cs_cold_tagseq_dwsample/clean_tagseq_', dw_code, '/counts_1.0k/all.counts.gene.tsv.gz')
    .x <- read.table(fpath, sep = '\t', header = TRUE)
    x <- .x[, -1]
    x <- x[, c(1, 5, 6, 2, 3, 4)]
    rownames(x) <- .x[, 1]
    colnames(x) <- LIBNAME
    y <- genemat2homeologmat(x)
    z <- y$A + y$B + y$D
    dfsz <- rbind(dfsz, data.frame(n_read = round(FQSIZE * as.numeric(dw_code)),
                                   n_gene = colSums(z > 0), libname = LIBNAME, dw_rate = dw_code))
}

fit <- nls(n_gene ~ a * log10(n_read) + b, start = c(a = 1, b = 1), data = dfsz)
a <- summary(fit)$coefficients[2]
b <- summary(fit)$coefficients[1]
newx <- seq(min(dfsz$n_read), max(dfsz$n_read), 100)
dffit <- data.frame(n_read = newx, n_gene = predict(fit, newdata = data.frame(n_read = newx)))
print(predict(fit, newdata = data.frame(n_read = c(1e6, 2e6, 3e6, 4e6))))
print(predict(fit, newdata = data.frame(n_read = c(1e6, 2e6, 3e6, 4e6)))/NGFF)


gp <- ggplot() +
        geom_point(aes(x = n_read, y = n_gene, color = dw_rate), data = dfsz) +
        geom_line(aes(x = n_read, y = n_gene), data = dffit) +
        scale_color_npg() +
        xlab('number of reads') + ylab('number of homeolog triad') +
        theme(legend.position = 'top') +
        labs(color = 'rate')
png(paste0('results/plots/downsampling_homeolog.png'), 900, 840, res = 220)
print(gp)
dev.off()



 




# 
# number of genes/homeologs detected in paired-end RNA-Seq (EAGLE-RC, IWGSC+1k)
#
load_counts <- function() {
    tags <- c('20181221.A-ZH_W2017_1_CS_2', '20181221.A-ZH_W2017_1_CS_3', '20181221.A-ZH_W2017_1_CS_4',
              '20181221.A-ZH_W2017_1_CS_cold_2', '20181221.A-ZH_W2017_1_CS_cold_3', '20181221.A-ZH_W2017_1_CS_cold_4')
    cnt <- NULL
    for (tag in tags) {
        cg <- NULL
        for (chr in c('chrA', 'chrB', 'chrD')) {
            hc <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts_1.0k/', tag,
                                '.', chr, '.gene.counts.txt'), header = TRUE, sep = '\t', stringsAsFactors = F)
            uc <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/uniq_counts_1.0k/', tag,
                                '.', chr, '.uniq.gene.counts.txt'), header = TRUE, sep = '\t', stringsAsFactors = F)
            cc <- rep(0, nrow(hc))
            names(cc) <- hc[, 1]
            cc <- hc[, 7] + uc[, 7]
            if (is.null(cg)) { cg <- cc } else {cg <- c(cg, cc)}
        }
        if (is.null(cnt)) {
            cnt <- matrix(0, nrow = length(cg), ncol = length(tags))
            colnames(cnt) <- tags
        }
        cnt[, tag] <- cg
    }
    
    cnt
}

cnt <- load_counts()
colSums(cnt > 0)
mean(colSums(cnt > 0))



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

    list(A = mat.a, B = mat.b, D = mat.d)
}

lbnms <- c('20181221.A-ZH_W2017_1_CS_2', '20181221.A-ZH_W2017_1_CS_3', '20181221.A-ZH_W2017_1_CS_4',
           '20181221.A-ZH_W2017_1_CS_cold_2', '20181221.A-ZH_W2017_1_CS_cold_3', '20181221.A-ZH_W2017_1_CS_cold_4')
fullseq_eagle <- NULL
for (lbnm in lbnms) {
    xa <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts_1.0k/', lbnm, '.chrA.gene.counts.txt'), sep = '\t', header = TRUE)[, c(1, 7)]
    xb <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts_1.0k/', lbnm, '.chrB.gene.counts.txt'), sep = '\t', header = TRUE)[, c(1, 7)]
    xd <- read.table(paste0('aaic_data/cs_cold_fullseq/fullseq_eaglerc/homeolog_counts_1.0k/', lbnm, '.chrD.gene.counts.txt'), sep = '\t', header = TRUE)[, c(1, 7)]
    colnames(xa) <- colnames(xb) <- colnames(xd) <- c('gene', 'count')
    xabd <- rbind(xa, xb, xd)
    if (is.null(fullseq_eagle)) {
        fullseq_eagle <- as.matrix(xabd[, 2])
    } else {
        fullseq_eagle <- cbind(fullseq_eagle, as.matrix(xabd[, 2]))
    }
    rownames(fullseq_eagle) <- xabd[, 1]
}
colnames(fullseq_eagle) <- paste0('rep', 1:6)
fullseq_eagle_h <- genemat2homeologmat(fullseq_eagle)

h <- fullseq_eagle_h$A + fullseq_eagle_h$B + fullseq_eagle_h$D
colSums(h > 0)
mean(colSums(h > 0))






