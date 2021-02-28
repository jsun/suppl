library(tidyverse)
library(ggsci)
library(ggExtra)
library(edgeR)
library(pROC)
options(stringsAsFactors = FALSE)





# fullseq IWGSC+1k counts

x <- read.table(paste0('data/fullseq/counts/counts.gene.ext_1.0k.tsv.gz'), sep = '\t', header = TRUE)
fullseq_hisat <- as.matrix(x[, -c(1:6)])
rownames(fullseq_hisat) <- x$Geneid
colnames(fullseq_hisat) <- c('ctrl_1', 'ctrl_2', 'ctrl_3', 'cold_1', 'cold_2', 'cold_3')
fullseq_all_zeros <- (rowSums(fullseq_hisat) == 0)

# tagseq IWGSC+1k counts

x <- read.table('data/tagseq/counts/cs.counts.gene.ext_1.0k.tsv.gz', sep = '\t', header = TRUE)
tagseq_hisat <- as.matrix(x[, -c(1:6)][, c(4, 5, 6, 1, 2, 3)])
rownames(tagseq_hisat) <- x$Geneid
colnames(tagseq_hisat) <- c('ctrl_1', 'ctrl_2', 'ctrl_3', 'cold_1', 'cold_2', 'cold_3')
tagseq_all_zeros <- (rowSums(tagseq_hisat) == 0)




run_edger <- function(x, group, method = 'LRT') {
    
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



group <- factor(rep(c('control', 'cold'), each = 3))
tagseq_deg <- run_edger(tagseq_hisat, group)
fullseq_deg <- run_edger(fullseq_hisat, group)

# remove all zero-counts in both approaches
tagseq_deg <- tagseq_deg[!(fullseq_all_zeros | tagseq_all_zeros), ]
fullseq_deg <- fullseq_deg[!(fullseq_all_zeros | tagseq_all_zeros), ]

dat <- data.frame(true = as.numeric(fullseq_deg$DEG), score = (1 - tagseq_deg$PValue))
roc_obj <- roc(true ~ score, data = dat, ci = FALSE)



# ROC 
roc_coordinates <- plot(roc_obj, identity = TRUE, print.thres = 'best',
                        print.thres.best.method = 'closest.topleft', legacy.axes = TRUE)
roc_coordinates <- data.frame(sensitivities = roc_coordinates$sensitivities,
                              specificities = roc_coordinates$specificities,
                              thresholds = roc_coordinates$thresholds)

roc_fig <- ggplot(roc_coordinates, aes(x = 1 - specificities, y = sensitivities)) +
                geom_point() + coord_fixed()

png(paste0('results/plots/DEG_ROC.png'), 1000, 900, res = 220)
print(roc_fig)
dev.off()






