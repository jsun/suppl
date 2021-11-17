
x <- read.table('meshmatrix.tsv.gz', header = TRUE)
y <- x[, -c(1, 2)]


s <- sort(colSums(y), decreasing  = TRUE)
head(sort(s, decreasing = TRUE), 10)

s2 <- s[s > 500]
barplot(s2, las=2)


