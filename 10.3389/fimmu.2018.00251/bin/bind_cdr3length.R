




f1 <- read.table('fugustatsNS_cm.sample.Fugu1.cdr3_prot_length.freq.tsv', stringsAsFactors = F)
f2 <- read.table('fugustatsNS_cm.sample.Fugu2.cdr3_prot_length.freq.tsv', stringsAsFactors = F)
f3 <- read.table('fugustatsNS_cm.sample.Fugu3.cdr3_prot_length.freq.tsv', stringsAsFactors = F)

lenname <- sort(unique(as.integer(c(rownames(f1), rownames(f2), rownames(f3)))))

f1tmp <- intersect(lenname, rownames(f1))
f2tmp <- intersect(lenname, rownames(f2))
f3tmp <- intersect(lenname, rownames(f3))


lenmat <- matrix(0, ncol = 3, nrow = max(lenname) + 1)
rownames(lenmat) <- as.character(0:max(lenname))


lenmat[f1tmp, 1] <- f1[f1tmp, 1]
lenmat[f2tmp, 2] <- f2[f2tmp, 1]
lenmat[f3tmp, 3] <- f3[f3tmp, 1]





