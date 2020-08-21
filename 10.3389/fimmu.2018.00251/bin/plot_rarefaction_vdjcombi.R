#!/usr/bin/env Rscript
## Run: Rscript --vanilla --slave plot_rarefaction_vdjcombi.R \
##                                            [prefix_sample1],[prefix_sample2],..., \
##                                            [sample name1],[sample name2],...,     \
##                                            name, format, width, height, res







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
fig.name   <- ifelse(is.na(argv[3]), NA, argv[3])
fig.format <- ifelse(is.na(argv[4]), "tiff", argv[4])
fig.width  <- ifelse(is.na(argv[5]), 8, as.numeric(argv[5]))
fig.height <- ifelse(is.na(argv[6]), 4, as.numeric(argv[6]))
fig.dpi    <- ifelse(is.na(argv[7]), 200, as.numeric(argv[7]))
log.file   <- paste0(fig.name, ".data.tsv")
if (is.na(fig.name)) stop("Figure name should be given.")




fancy_scientific <- function(l) {
     l <- format(l, scientific = TRUE)
     l <- gsub("^(.*)e", "'\\1'e", l)
     l <- gsub("e\\+", "%*%10^", l)
     l[1] <- "0"
     parse(text = l)
}



plotRarefactionCurve <- function(dat, plot.dat = NULL,
                                 xlab = "Sampling size", ylab = "VDJ combinations captured") {
    set.seed(2015)
    tmp.file = log.file
    e <- NULL

    sampling.sizes <- lapply(dat, function(x) {
        v1 <- c(0, 1e3, 2e3, 5e3, 1e4, 1.5e4, 2e4, 3e4, 4e4, 6e4, 8e4, 1e5, 1.5e5)
        if (2e5 < nrow(x)) {
            v2 <- c(seq(2e5, nrow(x), 2e5), nrow(x))
        } else {
            v2 <- nrow(x)
        }
        c(v1, v2)
    })

    libnames <- samplesize <- captured <- numeric()
    for (i in 1:length(sampling.sizes)) {
        libnames   <- c(libnames, rep(names(sampling.sizes)[i], length(sampling.sizes[[i]])))
        samplesize <- c(samplesize, sampling.sizes[[i]])
        captured   <- c(captured, rep(0, length(sampling.sizes[[i]])))
    }

    if (is.null(plot.dat)) {
        plot.dat <- data.frame(
            Sample = as.vector(libnames),
            samplesize = as.vector(samplesize),
            captured = as.vector(captured)
        )

        for (i in 1:nrow(plot.dat)) {
            idx <- (1:length(names(dat)))[names(dat) == plot.dat[i, 1]]
            resampling.times <- 1000
            resampled.combinations <- chao1 <- rep(0, resampling.times)
            for (j in 1:resampling.times) {
                sampled.index <- sample(1:nrow(dat[[idx]]), plot.dat[i, 2], replace = FALSE)
                sampled.data <- dat[[idx]][sampled.index, ]
                sampled.combinations <- apply(sampled.data[, 2:4], 1, paste0, collapse = "_")
                resampled.combinations[j] <- length(unique(sampled.combinations))
                if (length(sampled.combinations) > 0) 
                    e <- try(chao1[j] <- chao1(table(sampled.combinations)), silent = TRUE)
            }
            if (!is.null(e) && class(e) != "try-error")
                print(paste0("chao1 estimation (", plot.dat[i, 1], "):", mean(chao1)))
            plot.dat[i, 3] <- mean(resampled.combinations)
            write.table(plot.dat, file = tmp.file)
            print(paste0(round(i / nrow(plot.dat), 2) * 100, "% finished."))
        }
    }

    g <- ggplot(plot.dat, aes(x = samplesize, y = captured, colour = Sample, shape = Sample))
    g <- g + theme_bw()
    g <- g + geom_point() 
    g <- g + geom_line()
    g <- g + ylab(ylab)
    g <- g + xlab(xlab)
    g <- g + ylim(0, ceiling(max(plot.dat[, 3] / 1e2)) * 1e2)
    g <- g + scale_x_continuous(labels = fancy_scientific)
    g <- g + scale_colour_brewer(palette = "Set1")
    g

}

createDataFrame <- function(f) {
    dat <- read.table(f, skip = 1)
    dat <- dat[, -1]
    dat <- dat[!is.na(dat[, 2]), ]
    df1 <- df2 <- df3 <- df4 <- data.frame(id = NULL, v = NULL, d = NULL, j = NULL, stringsAsFactors = FALSE)
    for (i in 1:nrow(dat)) {
        if (i / nrow(dat) < 0.25) {
            df1 <- rbind(df1, data.frame(id = paste(i, 1:dat[i, 4]), v = dat[i, 1], d = dat[i, 2], j = dat[i, 3]))   
        } else if (i / nrow(dat) < 0.50){
            df2 <- rbind(df2, data.frame(id = paste(i, 1:dat[i, 4]), v = dat[i, 1], d = dat[i, 2], j = dat[i, 3]))   
        } else if (i / nrow(dat) < 0.75){
            df3 <- rbind(df3, data.frame(id = paste(i, 1:dat[i, 4]), v = dat[i, 1], d = dat[i, 2], j = dat[i, 3]))   
        } else {
            df4 <- rbind(df4, data.frame(id = paste(i, 1:dat[i, 4]), v = dat[i, 1], d = dat[i, 2], j = dat[i, 3]))   
        }
    }
    df <- rbind(rbind(df1, df2), rbind(df3, df4))
    df
}




fugu.dat <- vector("list", length(dat.files))
names(fugu.dat) <- sample.names
for (i in 1:length(dat.files)) fugu.dat[[i]] <- createDataFrame(dat.files[i])


if (!file.exists(log.file)) {
    fig <- plotRarefactionCurve(fugu.dat, NULL)
} else {
    plot.dat <- read.table(log.file)
    fig <- plotRarefactionCurve(fugu.dat, plot.dat)
}






if (fig.format == "tiff") {
    tiff(fig.name, width = fig.width, height = fig.height, units = "in", res = fig.dpi)
} else if (fig.format == "pdf") {
    pdf(fig.name, width = fig.width, height = fig.height)
} else if (fig.format == "png") {
    png(fig.name, width = fig.width, height = fig.height, res = fig.dpi)
}

plot(fig)
dev.off()


