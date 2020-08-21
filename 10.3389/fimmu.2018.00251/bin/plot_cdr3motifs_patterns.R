# R scripts for plotting CDR3 motifs.


## Paramater Settings
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
scale_colour_discrete <- function(...) scale_colour_brewer(palette = "Set1")
scale_fill_discrete <- function(...) scale_fill_brewer(palette = "Set1")


makeDataFrame <- function(motif.files) {
    y <- NULL
    for (f in motif.files) {
        d <- read.table(f, header  = TRUE)
        d <- d[- grep("\\*", d[, 1]), ]
        x <- d[, c(1, 5)]
        if (is.null(y)) y <- x
        else y <- merge(y, x, by.x = 1, by.y = 1, all = TRUE)
    }
    rownames(y) <- y[, 1]
    y <- y[, -1]
    y[is.na(y)] <- 0
    y <- y[rowSums(y) > 0, ]
    yp <- sweep(y, 2, 100 / colSums(y), "*")
    yp <- yp[apply(yp, 1, max) > 1, ]
     

}



fugu1.motif.f <- "fugu1.cm.dev.cdr3.fa.cdr3_motif"
fugu2.motif.f <- "fugu2.cm.dev.cdr3.fa.cdr3_motif"
fugu3.motif.f <- "fugu3.cm.dev.cdr3.fa.cdr3_motif"




x <- read.table('test', header = TRUE)
y <- x[apply(x[, -1], 1, max) > 1000, ]

m <- melt(y)

g <- ggplot(m, aes(x = motif, y = value))
g <- g + geom_bar(stat = 'identity', position = 'dodge')
g <- g + facet_grid(variable ~ .)
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g




