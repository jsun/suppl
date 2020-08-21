#!/usr/bin/env Rscript
## Run: Rscript --vanilla --slave study_capturerecapture_cdr3seqaa.R cm
##
##
##





## load package for estimating population sizes by capture-recapture study.
require(Rcapture)
set.seed(2016)



argv <- commandArgs(TRUE)
cgene <- argv[1]


## load data
fugu.dat <- list(
    `Fugu 1` = read.table(paste0("fugu1.", cgene, ".capturestudy.dataset.tsv"), stringsAsFactors = FALSE),
    `Fugu 2` = read.table(paste0("fugu2.", cgene, ".capturestudy.dataset.tsv"), stringsAsFactors = FALSE),
    `Fugu 3` = read.table(paste0("fugu3.", cgene, ".capturestudy.dataset.tsv"), stringsAsFactors = FALSE)
)




## variable to record the estimated population sizes.
estimated.population.sizes <- rep(0, length(fugu.dat))


## perform 10 times capture-recapture.
ntry <- 10
for (fuguidx in seq(fugu.dat)) {
    cluster.labels <- fugu.dat[[fuguidx]][, 1]
    capture.matrix <- matrix(0, nrow = length(unique(cluster.labels)), ncol = ntry)
    rownames(capture.matrix) <- unique(cluster.labels)
    for (j in 1:ntry) {
        sampled.cluster.labels <- sample(cluster.labels,
                                     round(length(cluster.labels) * runif(1, 0.2, 1)))
        capture.matrix[sampled.cluster.labels, j] <- 1
    }
    est <- closedp(capture.matrix)
    estimated.population.sizes[fuguidx] <- est$parameters$MthC[1]
}

estimated.population.sizes <- matrix(estimated.population.sizes, ncol = 3)
colnames(estimated.population.sizes) <- names(fugu.dat)
rownames(estimated.population.sizes) <- "estimated size"

write.table(estimated.population.sizes, file = paste0("CDR3aaseq.", cgene, ".estimated.population.sizes.tsv"),
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

