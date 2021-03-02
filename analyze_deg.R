library(tidyverse)
library(ggsci)
library(ggExtra)
library(edgeR)
library(pROC)
options(stringsAsFactors = FALSE)



load_eagle_counts <- function() {
    load_count <- function(lib_name) {
        xa1 <- read.table(paste0('data/fullseq/countseaglrc/', lib_name, '.chrA.counts.homeolog.txt.gz'), header = TRUE)
        xa2 <- read.table(paste0('data/fullseq/countseaglrc/', lib_name, '.chrA.counts.specific.txt.gz'), header = TRUE)
        xb1 <- read.table(paste0('data/fullseq/countseaglrc/', lib_name, '.chrB.counts.homeolog.txt.gz'), header = TRUE)
        xb2 <- read.table(paste0('data/fullseq/countseaglrc/', lib_name, '.chrB.counts.specific.txt.gz'), header = TRUE)
        xd1 <- read.table(paste0('data/fullseq/countseaglrc/', lib_name, '.chrD.counts.homeolog.txt.gz'), header = TRUE)
        xd2 <- read.table(paste0('data/fullseq/countseaglrc/', lib_name, '.chrD.counts.specific.txt.gz'), header = TRUE)
        
        gid <- c(xa1$Geneid, xb1$Geneid, xd1$Geneid)
        gcounts <- as.matrix(c(xa1[, 7] + xa2[, 7], xb1[, 7] + xb2[, 7], xd1[, 7] + xd2[, 7]))
        rownames(gcounts) <- gid
        gcounts
    }
    
    x <- NULL
    lib_names <- c('20181221.A-ZH_W2017_1_CS_2', '20181221.A-ZH_W2017_1_CS_3', '20181221.A-ZH_W2017_1_CS-4',
                   '20181221.A-ZH_W2017_1_CS_cold_2', '20181221.A-ZH_W2017_1_CS_cold_3', '20181221.A-ZH_W2017_1_CS_cold_4')
    for (lib_name in lib_names) {
        x <- cbind(x, load_count(lib_name))
    }
    colnames(x) <- c('ctrl_1', 'ctrl_2', 'ctrl_3', 'cold_1', 'cold_2', 'cold_3')
    
    x
}




# fullseq IWGSC+1k counts (HISAT2)

x <- read.table(paste0('data/fullseq/counts/counts.gene.ext_1.0k.tsv.gz'), sep = '\t', header = TRUE)
fullseq_hisat <- as.matrix(x[, -c(1:6)])
rownames(fullseq_hisat) <- x$Geneid
colnames(fullseq_hisat) <- c('ctrl_1', 'ctrl_2', 'ctrl_3', 'cold_1', 'cold_2', 'cold_3')
fullseqhisat_all_zeros <- (rowSums(fullseq_hisat) == 0)



# fullseq IWGSC+1k counts (EAGLERC)

fullseq_eagle <- load_eagle_counts()
fullseqeagle_all_zeros <- (rowSums(fullseq_eagle) == 0)



# tagseq IWGSC+1k counts

x <- read.table('data/tagseq/counts/cs.counts.gene.ext_1.0k.tsv.gz', sep = '\t', header = TRUE)
tagseq_hisat <- as.matrix(x[, -c(1:6)][, c(4, 5, 6, 1, 2, 3)])
rownames(tagseq_hisat) <- x$Geneid
colnames(tagseq_hisat) <- c('ctrl_1', 'ctrl_2', 'ctrl_3', 'cold_1', 'cold_2', 'cold_3')
tagseq_all_zeros <- (rowSums(tagseq_hisat) == 0)





run_edger <- function(x, group) {
    
    min.count <- mean(3 / colSums(x) * 1e6)

    y <- DGEList(counts = x, group = group)
    
    keep <- filterByExpr(y, min.count = min.count, min.prop = 0.5)
    y <- y[keep,,keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    
    d <- model.matrix(~ group)
    y <- estimateDisp(y, d)

    fit <- glmQLFit(y, d)
    qlf <- glmQLFTest(fit, coef = 2)
    deg_table <- as.data.frame(topTags(qlf, n = nrow(y)))
    
    deg_table_full <- data.frame(logFC = rep(NA, length = nrow(x)),
                                 logCPM = NA, F = NA, PValue = NA, FDR = NA)
    rownames(deg_table_full) <- rownames(x)
    deg_table_full[rownames(deg_table), ] <- deg_table
    deg_table_full$PValue[is.na(deg_table_full$PValue)] <- 1.0
    deg_table_full$FDR[is.na(deg_table_full$FDR)] <- 1.0
    deg_table_full$DEG <- ifelse(deg_table_full$FDR < 0.1, TRUE, FALSE)

    deg_table_full
}



# DEG ROC analysis (tagseq vs fullseq HISAT)

group <- factor(rep(c('control', 'cold'), each = 3))
tagseq_deg <- run_edger(tagseq_hisat, group)
fullseq_deg <- run_edger(fullseq_hisat, group)

# remove all zero-counts in both approaches
tagseq_deg <- tagseq_deg[!(fullseqhisat_all_zeros | tagseq_all_zeros), ]
fullseq_deg <- fullseq_deg[!(fullseqhisat_all_zeros | tagseq_all_zeros), ]

dat <- data.frame(true = as.numeric(fullseq_deg$DEG), score = (1 - tagseq_deg$PValue))
roc_obj <- roc(true ~ score, data = dat, ci = FALSE)

roc_coordinates <- plot(roc_obj, identity = TRUE, print.thres = 'best',
                        print.thres.best.method = 'closest.topleft', legacy.axes = TRUE)
roc_coordinates <- data.frame(sensitivities = roc_coordinates$sensitivities,
                              specificities = roc_coordinates$specificities,
                              thresholds = roc_coordinates$thresholds)

roc_fig <- ggplot(roc_coordinates, aes(x = 1 - specificities, y = sensitivities)) +
                geom_line() + coord_fixed()

png(paste0('results/plots/DEG_ROC_fullseqhisat.png'), 1000, 900, res = 220)
print(roc_fig)
dev.off()




# DEG ROC analysis (tagseq vs fullseq EAGLERC) 

group <- factor(rep(c('control', 'cold'), each = 3))
tagseq_deg <- run_edger(tagseq_hisat, group)
fullseqeagle_deg <- run_edger(fullseq_eagle, group)

# remove all zero-counts in both approaches and get the common genes
tagseq_deg <- tagseq_deg[!tagseq_all_zeros, ]
fullseqeagle_deg <- fullseqeagle_deg[!fullseqeagle_all_zeros, ]
common_genes <- intersect(rownames(tagseq_deg), rownames(fullseqeagle_deg))
tagseq_deg <- tagseq_deg[common_genes, ]
fullseqeagle_deg <- fullseqeagle_deg[common_genes, ]

dat <- data.frame(true = as.numeric(fullseqeagle_deg$DEG), score = (1 - tagseq_deg$PValue))
roc_obj <- roc(true ~ score, data = dat, ci = FALSE)

roc_coordinates <- plot(roc_obj, identity = TRUE, print.thres = 'best',
                        print.thres.best.method = 'closest.topleft', legacy.axes = TRUE)
roc_coordinates <- data.frame(sensitivities = roc_coordinates$sensitivities,
                              specificities = roc_coordinates$specificities,
                              thresholds = roc_coordinates$thresholds)

roc_fig <- ggplot(roc_coordinates, aes(x = 1 - specificities, y = sensitivities)) +
                geom_line() + coord_fixed()

png(paste0('results/plots/DEG_ROC_fullseqeagle.png'), 1000, 900, res = 220)
print(roc_fig)
dev.off()





