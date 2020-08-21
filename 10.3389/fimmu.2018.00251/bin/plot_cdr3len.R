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



## Functions

## Scientific labels of y-axis
fancy_scientific <- function(l) {
     l <- format(l, scientific = TRUE)
     l <- gsub("^(.*)e", "'\\1'e", l)
     l <- gsub("e\\+", "%*%10^", l)
     l[1] <- "0"
     parse(text = l)
}

read_data <- function(dat.files, sample.names, prob) {
    dat.matrix <- NULL
    for (dat.file in dat.files) {
        dat.f <- read.table(dat.file, skip = 1)
        if (is.null(dat.matrix)) {
            dat.matrix <- dat.f
        } else {
            dat.matrix <- merge(dat.matrix, dat.f, by.x = 1, by.y = 1, all = TRUE)
        }
    }
    
    dat.matrix[is.na(dat.matrix)] <- 0
    dat.matrix <- dat.matrix[dat.matrix[, 1] > 0, ]
    dat.upper <- dat.matrix[dat.matrix[, 1] <= 20, ]
    dat.lower <- dat.matrix[dat.matrix[, 1] > 20, ]
    dat.lower.sum <- colSums(dat.lower)
    
    dat.matrix <- rbind(dat.upper, dat.lower.sum)
    
    dat.matrix.name <- dat.matrix[, 1]
    dat.matrix.val  <- dat.matrix[, -1]
    if (prob == "TRUE") dat.matrix.val <- sweep(dat.matrix.val, 2, 100 / colSums(dat.matrix.val), "*")
    
    colnames(dat.matrix.val) <- c(sample.names)
    rownames(dat.matrix.val) <- c(dat.matrix.name[-length(dat.matrix.name)], "20<")
    dat.df <- melt(as.matrix(dat.matrix.val))
    colnames(dat.df) <- c("Length", "variable", "value")
    lv <- c(as.character(sort(unique(as.integer(as.character(dat.df$Length[dat.df$Length != "20<"]))))), "20<")
    
    dat.df$Length <- factor(dat.df$Length, lv)
    dat.df
}

## Main Functions


## Read data
dat.cdr3len <- read_data(dat.files, sample.names, prob)

g <- ggplot(dat.cdr3len, aes(x = Length, y = value, fill = variable, na.rm = TRUE))
g <- g + geom_bar(stat = "identity", position = "dodge")
g <- g + theme_bw()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + scale_fill_brewer(palette = "Set1")
if (prob == "FALSE") g <- g + scale_y_continuous(labels = fancy_scientific)
if (prob == "FALSE") g <- g + xlab("") + ylab("Frequency")
if (prob == "TRUE") g <- g + xlab("") + ylab("Proportion (%)")
g <- g + guides(fill = guide_legend(title = "Sample"))
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




