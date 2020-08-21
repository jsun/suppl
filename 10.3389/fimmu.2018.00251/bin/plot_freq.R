#!/usr/bin/env Rscript
## Run: Rscript --vanilla --slave plot_freq.R [prefix_sample1],[prefix_sample2],..., \
##                                            [sample name1],[sample name2],...,     \
##                                            prob, name, format, width, height, res







## Paramater Settings
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
scale_colour_discrete <- function(...) scale_colour_brewer(palette = "Set1")
scale_fill_discrete <- function(...) scale_fill_brewer(palette = "Set1")



## Argument Settings
argv <- commandArgs(TRUE)


dat.files    <- unlist(strsplit(argv[1], ","))
sample.names <- unlist(strsplit(argv[2], ","))
prob       <- ifelse(is.na(argv[3]), "FALSE", argv[3])
fig.name   <- ifelse(is.na(argv[4]), NA, argv[4])
fig.format <- ifelse(is.na(argv[5]), "tiff", argv[5])
fig.width  <- ifelse(is.na(argv[6]), 8, as.numeric(argv[6]))
fig.height <- ifelse(is.na(argv[7]), 4, as.numeric(argv[7]))
fig.dpi    <- ifelse(is.na(argv[8]), 200, as.numeric(argv[8]))

if (is.na(fig.name)) stop("Figure name should be given.")



# Functions

## Scientific labels of y-axis
fancy_scientific <- function(l) {
     l <- format(l, scientific = TRUE)
     l <- gsub("^(.*)e", "'\\1'e", l)
     l <- gsub("e\\+", "%*%10^", l)
     l[1] <- "0"
     parse(text = l)
}

read_data <- function(dat.files, sample.names, gene.cate = NULL, prob) {
    dat.matrix <- NULL
    for (dat.file in dat.files) {
        dat.f <- read.table(dat.file, skip = 1)
        dat.f[is.na(dat.f[, 1]), 1] <- "unidentifiable"
        if (is.null(dat.matrix)) {
            dat.matrix <- dat.f
        } else {
            dat.matrix <- merge(dat.matrix, dat.f, by.x = 1, by.y = 1, all = TRUE)
        }
    }
    dat.matrix[is.na(dat.matrix)] <- 0
    dat.matrix.na  <- dat.matrix[dat.matrix[, 1] == "unidentifiable", ]
    dat.matrix.val <- dat.matrix[dat.matrix[, 1] != "unidentifiable", ]
    dat.matrix <- rbind(dat.matrix.val[order(rowMeans(dat.matrix.val[, -1]), decreasing = TRUE), ], dat.matrix.na)
    dat.matrix.name <- dat.matrix[, 1]
    dat.matrix.val  <- dat.matrix[, -1]
    if (prob == "TRUE") dat.matrix.val <- sweep(dat.matrix.val, 2, 100 / colSums(dat.matrix.val), "*")
    dat.matrix <- cbind(dat.matrix.name, dat.matrix.val)
    colnames(dat.matrix) <- c("Gene", sample.names)
    dat.df <- melt(dat.matrix)
    dat.df$Gene <- factor(dat.df$Gene, levels = dat.matrix[, 1])
    dat.df$GeneType <- gene.cate
    dat.df
}

## Main Functions


## Read data
dat.df.v <- read_data(paste0(dat.files, ".v.freq.tsv"), sample.names, "V", prob)
dat.df.d <- read_data(paste0(dat.files, ".d.freq.tsv"), sample.names, "D", prob)
dat.df.j <- read_data(paste0(dat.files, ".j.freq.tsv"), sample.names, "J", prob)
dat.df <- rbind(rbind(dat.df.v, dat.df.d), dat.df.j)
dat.df$GeneType <- factor(dat.df$GeneType, levels = c("V", "D", "J"))

gene.lv <- c("v2.1", "v1.1", "v2.2", "v1.2", "v1.3",
             "v2.3", "v1.4", "v3.2", "v1.7", "v3.3",
             "v2.6", "v1.8", "v2.7", "v2.8", "v2.9", "v2.10",
             "v2.11", "v1.12", "v2.12", "v1.13", "v1.14", "v1.15",
             "v2.15", "v1.16", "v2.16", "v2.17", "v1.17", "v2.18",
             "v1.18", "v2.19", "v2.20", "v1.21")
if (fig.name == "Fig_Cm_VDJ_freq_hist.pdf") {
    gene.lv <- c(gene.lv, "Dm1", "Dm2", "Dm3", "Dm4", "Dm5", "Dm6", "ambiguous",
                 "Jm1", "Jm2", "Jm3", "Jm4", "Jm5")
} else {
    gene.lv <- c(gene.lv, "Dt", "ambiguous", "Jt")
}

dat.df$Gene <- factor(dat.df$Gene, levels = gene.lv)


g <- ggplot(dat.df, aes(x = Gene, y = value, fill = variable, na.rm = TRUE))
g <- g + geom_bar(stat = "identity", position = "dodge")
g <- g + theme_bw()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + scale_fill_brewer(palette = "Set1")
if (prob == "FALSE") g <- g + scale_y_continuous(labels = fancy_scientific)
if (prob == "FALSE") g <- g + xlab("") + ylab("Frequency")
if (prob == "TRUE") g <- g + xlab("") + ylab("Proportion (%)")
g <- g + guides(fill = guide_legend(title = "Sample"))
g <- g + facet_grid(. ~ GeneType,  scales = "free", space = "free")
g <- g + theme(plot.background = element_blank(),
               panel.background = element_blank(),
               axis.line = element_blank(), axis.ticks = element_blank(),
               strip.background = element_rect(fill = "white", colour = "white"),
               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



if (fig.format == "tiff") {
    tiff(fig.name, width = fig.width, height = fig.height, units = "in", res = fig.dpi)
} else if (fig.format == "pdf") {
    pdf(fig.name, width = fig.width, height = fig.height)
} else if (fig.format == "png") {
    png(fig.name, width = fig.width, height = fig.height, res = fig.dpi)
}

plot(g)
dev.off()




