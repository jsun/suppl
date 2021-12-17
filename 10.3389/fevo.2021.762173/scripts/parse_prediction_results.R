library(ggplot2)
library(ggsci)


calc.acc <- function(d, top.ranks = 1) {
    
    d.labels <- d[, 1]
    d.probs <- d[, 3:ncol(d)]
    d.ranks <- t(apply(-1 * d.probs, 1, rank))
    
    class.labels <- colnames(d.ranks)
    
    acc <- c(NA)
    for (i in 1:nrow(d.ranks)) {
        pred.labels <- class.labels[d.ranks[i, ] <= top.ranks]
        pred.labels <- gsub('(_M|_F)', '', pred.labels)
        acc <- c(acc, d.labels[i] %in% pred.labels)
    }
    acc <- acc[-1]
    
    acc
}


load.prediction.results <- function(f) {

    d <- read.table(f, sep = '\t', header = TRUE)
    d[, 1] <- gsub('(_M|_F)', '', sapply(strsplit(d[, 1], '/'),'[', 7))
    d[, 2] <- gsub('(_M|_F)', '', d[, 2])
    d

}


d <- load.prediction.results('results/dragonfly_whitebg.whitebg.tsv')
df <- data.frame(top1 = calc.acc(d, top.ranks = 1),
                 top3 = calc.acc(d, top.ranks = 3),
                 top5 = calc.acc(d, top.ranks = 5))
colSums(df) / nrow(df)
d <- load.prediction.results('results/dragonfly_whitebg.fieldphoto.tsv')
df <- data.frame(top1 = calc.acc(d, top.ranks = 1),
                 top3 = calc.acc(d, top.ranks = 3),
                 top5 = calc.acc(d, top.ranks = 5))
colSums(df) / nrow(df)





