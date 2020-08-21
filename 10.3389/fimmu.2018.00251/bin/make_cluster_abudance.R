#
library(ggplot2)
library(stringr)
x <- 'classabudance.ct.txt'
x <- 'classabudance.cm.txt'



d <- read.table(x, sep = "\t", stringsAsFactors = FALSE)

colClusterType = 2
colClusterSize = 3
colFuguId = 4
colCDR3AA = 5
colVdel = 6
colJdel = 7



df.C1.clusterSize <- data.frame(freq = table(d[is.C1, 1]), type = "shared with 1 Fugu")
df.C2.clusterSize <- data.frame(freq = table(d[is.C2, 1]), type = "shared with 2 Fugu")
df.C3.clusterSize <- data.frame(freq = table(d[is.C3, 1]), type = "shared with 3 Fugu")
df.all.clusterSize <- rbind(df.C1.clusterSize, df.C2.clusterSize, df.C3.clusterSize)
df.all.clusterSize[, 2] <- log10(df.all.clusterSize[, 2])

g.clusterSize <- ggplot(df.all.clusterSize, aes(x = freq.Freq, fill = type))
g.clusterSize <- g.clusterSize + geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.05)
g.clusterSize <- g.clusterSize + xlim(0, 3) + ylim(0, 30000)
g.clusterSize <- g.clusterSize + ylab("frequence") + xlab("log10(captured times)")
png("clusterabundance.caputredtimes.png")
plot(g.clusterSize)
dev.off()


df.C1.aaLen <- data.frame(len = str_length(d[is.C1, colCDR3AA]), type = "shared with 1 Fugu")
df.C2.aaLen <- data.frame(len = str_length(d[is.C2, colCDR3AA]), type = "shared with 2 Fugu")
df.C3.aaLen <- data.frame(len = str_length(d[is.C3, colCDR3AA]), type = "shared with 3 Fugu")
df.all.aaLen<- rbind(df.C1.aaLen, df.C2.aaLen, df.C3.aaLen)

g.aaLen <- ggplot(df.all.aaLen, aes(x = len, fill = type))
g.aaLen <- g.aaLen + geom_histogram(alpha = 0.5, position = "identity", binwidth = 1)
g.aaLen <- g.aaLen + xlim(0, 30)
png("clusterabundance.cdr3aalength.png")
plot(g.aaLen)
dev.off()




df.C1.vdelLen <- data.frame(len = str_length(d[is.C1, colVdel]), type = "shared with 1 Fugu")
df.C2.vdelLen <- data.frame(len = str_length(d[is.C2, colVdel]), type = "shared with 2 Fugu")
df.C3.vdelLen <- data.frame(len = str_length(d[is.C3, colVdel]), type = "shared with 3 Fugu")
df.all.vdelLen<- rbind(df.C1.vdelLen, df.C2.vdelLen, df.C3.vdelLen)

g.vdelLen <- ggplot(df.all.vdelLen, aes(x = len, fill = type))
g.vdelLen <- g.vdelLen + geom_histogram(alpha = 0.5, position = "identity", binwidth = 1)
g.vdelLen <- g.vdelLen + xlim(-1, 25)
png("clusterabundance.VDellength.png")
plot(g.vdelLen)
dev.off()




df.C1.jdelLen <- data.frame(len = str_length(d[is.C1, colJdel]), type = "shared with 1 fugu")
df.C2.jdelLen <- data.frame(len = str_length(d[is.C2, colJdel]), type = "shared with 2 fugu")
df.C3.jdelLen <- data.frame(len = str_length(d[is.C3, colJdel]), type = "shared with 3 fugu")
df.all.jdelLen<- rbind(df.C1.jdelLen, df.C2.jdelLen, df.C3.jdelLen)

g.jdelLen <- ggplot(df.all.jdelLen, aes(x = len, fill = type))
g.jdelLen <- g.jdelLen + geom_histogram(alpha = 0.5, position = "identity", binwidth = 1)
g.jdelLen <- g.jdelLen + xlim(-1, 25)
png("clusterabundance.JDellength.png")
plot(g.jdelLen)
dev.off()












ScientificNotation <- function(l) {
     l <- format(l, scientific = TRUE)
     l <- gsub("^(.*)e", "'\\1'e", l)
     l <- gsub("e\\+", "%*%10^", l)
     l[1] <- "0"
     parse(text = l)
}




x <- 'classabudance.cm.2.txt'
x <- 'classabudance.ct.2.txt'

d <- read.table(x, sep = "\t", stringsAsFactors = FALSE)
colnames(d) <- c("cluster", "type", "values")

d1 <- d[d$type == 1, ]
d2 <- d[d$type == 2, ]
d3 <- d[d$type == 3, ]

d$type[d$type == 1] <- "shared with 1 fugu"
d$type[d$type == 2] <- "shared with 2 fugu"
d$type[d$type == 3] <- "shared with 3 fugu"

d$values <- log10(d$values)

g <- ggplot(d, aes(x = values, group = type))
g <- g + geom_histogram(position = "identity", binwidth = 0.25)
g <- g + xlab("number of unique nt sequences per a CDR3 cluster") + ylab("frequences")
g <- g + theme_bw()
g <- g + theme(strip.background = element_rect(fill = "white", colour = "white"))
g <- g + scale_x_continuous(breaks = c(0:3), labels = c(1, 10, 100, 1000))
g <- g + facet_grid(type ~ ., scales = "free_y")
pdf("number_of_shared_nt_per_fugu.ct.pdf", 4, 6)
plot(g)
dev.off()




library(ggtern)
x <- 'classabudance.cm.2.txt'
#x <- 'classabudance.ct.2.txt'

xp <- paste0(x, '.public')

d <- read.table(xp, sep = "\t", stringsAsFactors = FALSE)
colnames(d) <- c("cluster", "fugu1", "fugu2", "fugu3", "total")
rownames(d) <- d$cluster
e <- read.table(x, sep = "\t", stringsAsFactors = FALSE)
colnames(e) <- c("cluster", "type", "values")
e <- e[e$type == 3, ]

m <- median(e$values[e$type == 3])


.d <- d[e[e$values == m, 'cluster'], ]
.d$fugu1 <- .d$fugu1 + rnorm(nrow(.d), 0.1, 0.2)
.d$fugu2 <- .d$fugu2 + rnorm(nrow(.d), 0.1, 0.2)
.d$fugu3 <- .d$fugu3 + rnorm(nrow(.d), 0.1, 0.2)
g <- ggtern(.d, aes(fugu1, fugu2, fugu3))
g <- g + geom_point(alpha = 0.3, size = 1.4)
g <- g + theme_showarrows()
g

pdf("ternplot.ct.pdf", 5, 5)
plot(g)
dev.off()







