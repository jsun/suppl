#!/usr/bin/env Rscript
## Run: Rscript --vanilla --slave study_cdr3aaseq.R cm cm_output_prefix
##



## load package for estimating population sizes by capture-recapture study.
#library(Rcapture)
library(fossil)
library(VennDiagram)
library(ggplot2)
set.seed(2016)


argv <- commandArgs(TRUE)
cgene <- argv[1]
prefix <- argv[2]


## load data
all.dat <- read.table(paste0("fuguall.", cgene, ".capturestudy.dataset.tsv"),
                         stringsAsFactors = FALSE)

fugu1.dat <- all.dat[all.dat[, 2] == "fugu1", c(1, 3)]
fugu2.dat <- all.dat[all.dat[, 2] == "fugu2", c(1, 3)]
fugu3.dat <- all.dat[all.dat[, 2] == "fugu3", c(1, 3)]
fugu.dat <- list(`Fugu 1` = fugu1.dat, `Fugu 2` = fugu2.dat, `Fugu 3` = fugu3.dat)


venn.cluster <- list(`Fugu 1` = unique(fugu1.dat[, 1]),
                     `Fugu 2` = unique(fugu2.dat[, 1]),
                     `Fugu 3` = unique(fugu3.dat[, 1]))
venn.diagram(venn.cluster,
             filename = paste0(prefix, ".CDR3aaseq.", cgene, ".venn.tiff"),
             lwd = 0.5, margin = rep(.08, 4), height = 6, width = 6, units = "in",
             imagetype = "tiff", fontfamily = "sans", cex = 1,
             cat.fontfamily = "sans", cat.cex = 1,
             fill = c("#E41A1C", "#377EB8", "#4DAF4A"))
log.files <- list.files(dirname(prefix), pattern = paste0(basename(prefix), ".+log"))
file.remove(paste0(dirname(prefix), "/", log.files))






## sampling-resampling for estimating CDR3 AA population sizes
cdr3seq <- c(0, 0, 0, 0)
cdr3cls <- c(0, 0, 0, 0)
cdr3pop <- c(0, 0, 0, 0)

for (i in seq(fugu.dat)) {
    cdr3seq[i] <- nrow(fugu.dat[[i]])
    cdr3cls[i] <- length(table(fugu.dat[[i]][, 1]))
    cdr3pop[i] <- ACE(table(fugu.dat[[i]][, 1]))
}

pooled.fugu.dat <- c(fugu.dat[[1]][, 1], fugu.dat[[2]][, 1], fugu.dat[[3]][, 1])
cdr3seq[4] <- length(pooled.fugu.dat)
cdr3cls[4] <- length(table(pooled.fugu.dat))
cdr3pop[4] <- ACE(table(pooled.fugu.dat))

df <- cbind(cdr3seq, cdr3cls, round(cdr3pop))
colnames(df) <- c("#Seq", "#Cluster", "#Population")
print(df)



tb <- table(pooled.fugu.dat)
tb.sorted <- tb[order(tb, decreasing = TRUE)]
tb.sorted <- tb.sorted / sum(tb.sorted)

df2 <- data.frame(ab_rank = 1:length(tb.sorted),
                  rel_ab = as.numeric(tb.sorted))
g <- ggplot(df2, aes(x = ab_rank, y = rel_ab))
g <- g + geom_line()
g <- g + scale_y_log10()
g <- g + ylab('relative abundance (log10)')
g <- g + xlab('abundance rank (log10)')


#

pdf(paste0('abundance-rank-curve.', cgene, '.cdr3.pdf'), 4, 4)
plot(g)
dev.off()



## variable to record the estimated population sizes.
#estimated.population.sizes <- rep(0, length(fugu.dat))
## get all clusters
#cluster.names <- NULL
#for (fuguidx in seq(fugu.dat)) {
#    cluster.names <- union(cluster.names, fugu.dat[[fuguidx]][, 1])
#}
#capture.matrix <- matrix(0, nrow = length(cluster.names), ncol = length(fugu.dat))
#rownames(capture.matrix) <- cluster.names
#colnames(capture.matrix) <- names(fugu.dat)
# set captured clusters of each fugu
#for (fuguidx in seq(fugu.dat)) {
#    capture.matrix[fugu.dat[[fuguidx]][, 1], fuguidx] <- 1
#}
#est <- closedp(capture.matrix)
#est.pop.size <- est$parameters$MthC[1]
#pop.mat <- c(sapply(fugu.dat, function(x){length(unique(x[, 1]))}), est.pop.size)
#names(pop.mat) <- c(names(fugu.dat), "estimated size")
#write.table(pop.mat, file = paste0(prefix, ".CDR3aaseq.", cgene, ".estimated.population.sizes.tsv"),
#            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)







