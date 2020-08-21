




f1 <- read.table('fugustats_ct.Fugu1.vdj.freq.tsv', stringsAsFactors = F)
f2 <- read.table('fugustats_ct.Fugu2.vdj.freq.tsv', stringsAsFactors = F)
f3 <- read.table('fugustatS_ct.Fugu3.vdj.freq.tsv', stringsAsFactors = F)
f1 <- f1[!is.na(f1$d), ]
f2 <- f2[!is.na(f2$d), ]
f3 <- f3[!is.na(f3$d), ]
rownames(f1) <- apply(f1[, 1:3], 1, paste0, collapse = "__")
rownames(f2) <- apply(f2[, 1:3], 1, paste0, collapse = "__")
rownames(f3) <- apply(f3[, 1:3], 1, paste0, collapse = "__")
allname <- unique(c(rownames(f1), rownames(f2), rownames(f3)))
alldat <- matrix(0, ncol = 3, nrow = length(allname))
rownames(alldat) <- allname
f1tmp <- intersect(allname, rownames(f1))
f2tmp <- intersect(allname, rownames(f2))
f3tmp <- intersect(allname, rownames(f3))
alldat[f1tmp, 1] <- f1[f1tmp, 4]
alldat[f2tmp, 2] <- f2[f2tmp, 4]
alldat[f3tmp, 3] <- f3[f3tmp, 4]
alldat2 <- data.frame(t(as.data.frame(strsplit( allname, "__"))), alldat)
rownames(alldat2) <- NULL



