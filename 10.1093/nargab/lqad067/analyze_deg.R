library(tidyverse)
library(ggsci)
library(ggExtra)
library(edgeR)
library(topGO)
library(pROC)
options(stringsAsFactors = FALSE)


LIBNAME <- c('a_ctrl_1', 'a_ctrl_2', 'a_ctrl_3', 'b_cold_1', 'b_cold_2', 'b_cold_3')


load_eagle_counts <- function() {
    load_count <- function(lib_name) {
        xa1 <- read.table(paste0('data/fullseq/countseaglerc/',
                                 lib_name, '.chrA.counts.homeolog.txt.gz'), header = TRUE)
        xa2 <- read.table(paste0('data/fullseq/countseaglerc/',
                                 lib_name, '.chrA.counts.specific.txt.gz'), header = TRUE)
        xb1 <- read.table(paste0('data/fullseq/countseaglerc/',
                                 lib_name, '.chrB.counts.homeolog.txt.gz'), header = TRUE)
        xb2 <- read.table(paste0('data/fullseq/countseaglerc/',
                                 lib_name, '.chrB.counts.specific.txt.gz'), header = TRUE)
        xd1 <- read.table(paste0('data/fullseq/countseaglerc/',
                                 lib_name, '.chrD.counts.homeolog.txt.gz'), header = TRUE)
        xd2 <- read.table(paste0('data/fullseq/countseaglerc/',
                                  lib_name, '.chrD.counts.specific.txt.gz'), header = TRUE)
        
        gid <- c(xa1$Geneid, xb1$Geneid, xd1$Geneid)
        gcounts <- as.matrix(c(xa1[, 7] + xa2[, 7], xb1[, 7] + xb2[, 7], xd1[, 7] + xd2[, 7]))
        rownames(gcounts) <- gid
        gcounts
    }
    
    x <- NULL
    lib_names <- c('20181221.A-ZH_W2017_1_CS_2', '20181221.A-ZH_W2017_1_CS_3',
                   '20181221.A-ZH_W2017_1_CS_4', '20181221.A-ZH_W2017_1_CS_cold_2',
                   '20181221.A-ZH_W2017_1_CS_cold_3', '20181221.A-ZH_W2017_1_CS_cold_4')
    for (lib_name in lib_names) {
        x <- cbind(x, load_count(lib_name))
    }
    colnames(x) <- LIBNAME
    
    x
}


bind_gene_symbols <- function(x) {
    d1 <- read.table('data/gene2symbol.tsv.gz', header = TRUE, sep = '\t')
    d2 <- read.table('data/gene2desc.tsv.gz', header = TRUE, sep = '\t', quote = '')
    d3 <- read.table('data/gene2go.tsv.gz', header = TRUE, sep = '\t', quote = '')
    dict1 <- d1[, 2]
    names(dict1) <- d1[, 1]
    dict2 <- d2[, 2]
    names(dict2) <- d2[, 1]
    dict3 <- d3[, 2]
    names(dict3) <- d3[, 1]
    y <- data.frame(x,
                    symbol = dict1[rownames(x)],
                    desc = dict2[rownames(x)],
                    GOBP = dict3[rownames(x)])
    y
}


run_edger <- function(x, group, plantid = NULL) {
    min.count <- mean(3 / colSums(x) * 1e6)
    y <- DGEList(counts = x, group = group)
    keep <- filterByExpr(y, min.count = min.count, min.prop = 0.5)
    y <- y[keep,,keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    if (is.null(plantid)) {
        d <- model.matrix(~ group)
        y <- estimateDisp(y, d)
        fit <- glmQLFit(y, d)
        qlf <- glmLRT(fit, coef = 2)
    } else {
        d <- model.matrix(~ plantid + group)
        y <- estimateDisp(y, d)
        fit <- glmQLFit(y, d)
        qlf <- glmLRT(fit, coef = 4)
    }
    deg_table <- as.data.frame(topTags(qlf, n = nrow(y)))
    
    deg_table_full <- data.frame(logFC = rep(NA, length = nrow(x)),
                                 logCPM = NA, F = NA, PValue = NA, FDR = NA)
    rownames(deg_table_full) <- rownames(x)
    deg_table_full[rownames(deg_table), ] <- deg_table
    deg_table_full$PValue[is.na(deg_table_full$PValue)] <- 1.0
    deg_table_full$FDR[is.na(deg_table_full$FDR)] <- 1.0
    deg_table_full$DEG <- ifelse(deg_table_full$FDR < 0.1, TRUE, FALSE)
    deg_table_full
}


run_topgo <- function(allgene, siggene) {
    geneID2GO <- readMappings('data/gene2go.tsv.gz')

    allgene <- as.character(allgene)
    siggene <- as.character(siggene)
    allgene.f <- as.numeric(allgene %in% siggene)
    names(allgene.f) <- allgene
    topgo.obj <- new('topGOdata', ontology = 'BP', allGenes = allgene.f,
                     geneSel = function(x) x > 0,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    elim_fisher    <- runTest(topgo.obj, algorithm = 'elim', statistic = 'fisher')
    topgotable  <- GenTable(topgo.obj, numChar = 200,
                            elimFisher = elim_fisher,
                            orderBy = 'elimFisher', ranksOf = 'classicFisher', topNodes = 600)
    #topgotable <- topgotable[topgotable[, 3] > 10, ]
    #topgotable <- topgotable[topgotable[, 3] < 500, ]
    topgotable <- data.frame(topgotable)
    topgotable
}

deg2go <- function(x, cutoff = 500, output_prefix = NULL) {
    x <- x[!is.na(x$F),]
    x <- x[order(x$PValue), ]
    x <- bind_gene_symbols(x)
    
    g_background <- rownames(x)
    g_deg <- NA
    if (cutoff <= 1) {
        g_deg <- rownames(x)[x$FDR < cutoff & abs(x$logFC) > 1]
    } else {
        g_deg <- rownames(x)[1:cutoff]
    }
    
    x_go <- run_topgo(g_background, g_deg)
    
    if (!is.null(output_prefix)) {
        write.table(x, file = paste0(output_prefix, '_DEG.txt'), sep = '\t',
                    row.names = TRUE, col.names = TRUE)
        write.table(x_go, file = paste0(output_prefix, '_DEGGO.txt'), sep = '\t',
                    row.names = TRUE, col.names = TRUE)
    }
}



# fullseq IWGSC+1k counts (HISAT2)
x <- read.table('data/fullseq/counts/counts.gene.ext_1.0k.tsv.gz',
                sep = '\t', header = TRUE)
fullseq_hisat <- as.matrix(x[, -c(1:6)])
rownames(fullseq_hisat) <- x$Geneid
colnames(fullseq_hisat) <- LIBNAME
fullseqhisat_all_zeros <- (rowSums(fullseq_hisat) == 0)



# downsampled fullseq IWGSC+1k counts (HISAT2)
x <- read.table('data/dwfullseq/counts/counts.gene.ext_1.0k.tsv.gz',
                sep = '\t', header = TRUE)
dwfullseq_hisat <- as.matrix(x[, -c(1:6)])
rownames(dwfullseq_hisat) <- x$Geneid
colnames(dwfullseq_hisat) <- LIBNAME
dwfullseqhisat_all_zeros <- (rowSums(dwfullseq_hisat) == 0)



# fullseq IWGSC+1k counts (EAGLERC)
fullseq_eagle <- load_eagle_counts()
fullseqeagle_all_zeros <- (rowSums(fullseq_eagle) == 0)



# tagseq IWGSC+1k counts
x <- read.table('data/tagseq/counts/cs.counts.gene.ext_1.0k.tsv.gz',
                 sep = '\t', header = TRUE)
tagseq_hisat <- as.matrix(x[, -c(1:6)][, c(4, 5, 6, 1, 2, 3)])
rownames(tagseq_hisat) <- x$Geneid
colnames(tagseq_hisat) <- LIBNAME
tagseq_all_zeros <- (rowSums(tagseq_hisat) == 0)




# DEG analysis (tagseq vs fullseq), HISAT2
group <- factor(rep(c('a_control', 'b_cold'), each = 3))
plantid <- rep(c('1', '2', '3'), times = 2)
tagseq_deg <- run_edger(tagseq_hisat, group, plantid)
fullseq_deg <- run_edger(fullseq_hisat, group, plantid)

td <- tagseq_deg[order(tagseq_deg$PValue), ]
fd <- fullseq_deg[order(fullseq_deg$PValue), ]

td_degid <- rownames(td)[td$FDR < 0.05 & abs(td$logFC) > 1]
fd_degid <- rownames(fd)[fd$FDR < 0.05 & abs(fd$logFC) > 1]

td <- bind_gene_symbols(td)
fd <- bind_gene_symbols(fd)

deg2go(tagseq_deg, 0.1, 'results/tagseq_hisat2_')
deg2go(fullseq_deg, 0.1, 'results/fullseq_hisat2_')

common_deg <- td[td_degid[td_degid %in% fd_degid], ]
write.table(common_deg, 'results/tagseq_fullseq_shared_DEG.tsv',
            sep = '\t', col.names = TRUE, row.names = TRUE)

print(length(td_degid))
print(length(fd_degid))
print(sum(td_degid %in% fd_degid))



go_tagseq <- read.table('results/tagseq_hisat2__DEGGO.txt', quote = '"', sep = '\t')
go_fullseq <- read.table('results/fullseq_hisat2__DEGGO.txt', quote = '"', sep = '\t')

go_tagseq_sig <- go_tagseq$Term[go_tagseq$elimFisher < 0.05]
go_fullseq_sig <- go_fullseq$Term[go_fullseq$elimFisher < 0.05]





# DEG analysis (tagseq vs dwfullseq), HISAT2
group <- factor(rep(c('a_control', 'b_cold'), each = 3))
plantid <- rep(c('1', '2', '3'), times = 2)
tagseq_deg <- run_edger(tagseq_hisat, group, plantid)
dwfullseq_deg <- run_edger(dwfullseq_hisat, group, plantid)

td <- tagseq_deg[order(tagseq_deg$PValue), ]
fd <- dwfullseq_deg[order(dwfullseq_deg$PValue), ]

td_degid <- rownames(td)[td$FDR < 0.05 & abs(td$logFC) > 1]
fd_degid <- rownames(fd)[fd$FDR < 0.05 & abs(fd$logFC) > 1]

td <- bind_gene_symbols(td)
fd <- bind_gene_symbols(fd)

deg2go(tagseq_deg, 0.1, 'results/tagseq_hisat2_')
deg2go(dwfullseq_deg, 0.1, 'results/dwfullseq_hisat2_')

common_deg <- td[td_degid[td_degid %in% fd_degid], ]
write.table(common_deg, 'results/tagseq_dwfullseq_shared_DEG.tsv',
            sep = '\t', col.names = TRUE, row.names = TRUE)

print(length(td_degid))
print(length(fd_degid))
print(sum(td_degid %in% fd_degid))




data4plots <- function(x) {
    group <- rep(c('control', 'cold'), each = 3)
    design  <- model.matrix(~ group)

    min.count <- mean(3 / colSums(x) * 1e6)
    y <- DGEList(counts = x, group = group)
    keep <- filterByExpr(y, min.count = min.count, min.prop = 0.5)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    
    y
}

d1 <- data4plots(dwfullseq_hisat)
d2 <- data4plots(tagseq_hisat)

png('results/maplot_rep1_dwfullseq.png', 1000, 1000, res = 220)
plotSmear(d1, ylim = c(-10, 10))
dev.off()

png('results/maplot_rep1_tagseq.png', 1000, 1000, res = 220)
plotSmear(d2, ylim = c(-10, 10))
dev.off()


png('results/BCAplot_rep1_dwfullseq.png', 1000, 1000, res = 220)
plotBCV(d1, ylim = c(0, 1.7))
dev.off()

png('results/BCAplot_rep1_tagseq.png', 1000, 1000, res = 220)
plotBCV(d2, ylim = c(0, 1.7))
dev.off()





