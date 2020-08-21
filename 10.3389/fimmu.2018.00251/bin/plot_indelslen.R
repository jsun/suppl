#!/usr/bin/env Rscript





## Paramater Settings
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
scale_colour_discrete <- function(...) scale_colour_brewer(palette = "Set1")
scale_fill_discrete <- function(...) scale_fill_brewer(palette = "Set1")



## Argument Settings
argv <- commandArgs(TRUE)


dat.file   <- argv[1]
prob       <- ifelse(is.na(argv[2]), "FALSE", argv[2])
fig.name   <- ifelse(is.na(argv[3]), NA, argv[3])
fig.format <- ifelse(is.na(argv[4]), "tiff", argv[4])
fig.width  <- ifelse(is.na(argv[5]), 8, as.numeric(argv[5]))
fig.height <- ifelse(is.na(argv[6]), 4, as.numeric(argv[6]))
maxlen    <- ifelse(is.na(argv[7]), 20, as.numeric(argv[7]))

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

read_data <- function(dat.file, prob) {
    dat.matrix <- NULL
    dat.matrix <- read.table(dat.file)
    dat.matrix[is.na(dat.matrix)] <- 0
    dat.matrix <- cbind(as.integer(rownames(dat.matrix)), dat.matrix)
    
    dat.upper <- dat.matrix[dat.matrix[, 1] <= maxlen, ]
    dat.lower <- dat.matrix[dat.matrix[, 1] > maxlen, ]
    dat.lower.sum <- colSums(dat.lower)
    dat.matrix <- rbind(dat.upper, dat.lower.sum)
    
    dat.matrix.name <- dat.matrix[, 1]
    dat.matrix.val  <- dat.matrix[, -1]
    if (prob == "TRUE") dat.matrix.val <- sweep(dat.matrix.val, 2, 100 / colSums(dat.matrix.val), "*")
    
    colnames(dat.matrix.val) <- c("Fugu 1", "Fugu 2", "Fugu 3")
    rownames(dat.matrix.val) <- c(dat.matrix.name[-length(dat.matrix.name)], paste0(maxlen, "<"))
    dat.df <- melt(as.matrix(dat.matrix.val))
    colnames(dat.df) <- c("Length", "variable", "value")
    lv <- c(as.character(sort(unique(as.integer(as.character(dat.df$Length[dat.df$Length != paste0(maxlen, "<")]))))),
            paste0(maxlen, "<"))
    
    dat.df$Length <- factor(dat.df$Length, lv)
    dat.df
}

## Main Functions


## Read data
dat.df <- read_data(dat.file, prob)

g <- ggplot(dat.df, aes(x = Length, y = value, fill = variable, na.rm = TRUE))
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




