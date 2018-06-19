#
# main analhysis code

library(MASS)
library(edgeR)
library(org.At.tair.db)
library(topGO)
library(ggplot2)
library(ggsci)
library(reshape2)
library(openxlsx)
library(gplots)
library(RColorBrewer)
library(eulerr)

FIG_PATH <- '~/research/cardamine3/fig'



COLSFUNC <- function(n = 100) {
    #.colfunc <- colorRampPalette(c("#FFFFFF", "#D6604D", "#67000D"))
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu")))
    .colfunc(n)
}
COLSFUNCSP <- function(n = 100) {
    #.colfunc <- colorRampPalette(c("#FFFFFF", "#D6604D", "#67000D"))
    .colfunc <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
    .colfunc(n)
}
COLSFUNC2 <- function(n = 100) {
    .colfunc <- colorRampPalette(c("#053061", "#4393C3", "#FFFFFF", "#D6604D", "#67001F"))
    .colfunc(n)
}
COLSFUNCC2 <- function(n = 100) {
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu"))[9:5])
    .colfunc(n)
}



WriteExcel <- function(data.list, file.name = NULL) {
    dat <- list()
    nam <- NULL
    for (i in names(data.list)) {
        if (class(data.list[[i]]) == "list") {
            for (j in names(data.list[[i]])) {
                nam <- c(nam, paste0(i, ".", j))
                dat <- c(dat, list(data.frame(Rowname = rownames(data.list[[i]][[j]]), data.list[[i]][[j]])))
            }
        } else {
            nam <- c(nam, i)
            dat <- c(dat, list(data.frame(Rowname = rownames(data.list[[i]]), data.list[[i]])))
        }
    }
    names(dat) <- nam
    if (is.null(file.name)) stop("Need file name.")
    if (is.null(names(dat))) stop("Need list names.")
    if (file.exists(file.name)) file.remove(file.name)
    gc()
    write.xlsx(data.list, file = file.name)
    gc()
}




# read count and FPKM data
load.expression <- function(D) {
    aaic.data.path <- '~/research/cardamine3/data/hir_std'
    count.h.files <- list.files(aaic.data.path, pattern = 'orig.txt', full.names = TRUE)
    count.a.files <- list.files(aaic.data.path, pattern = 'other.txt', full.names = TRUE)
    count.u.files <- list.files(aaic.data.path, pattern = 'common.txt', full.names = TRUE)
    
    .get.counts <- function(count.files) {
        gene.names <- NULL
        read.count <- NULL
        for (f in sort(count.files)) {
            x <- read.table(f, sep = '\t', header = TRUE)
            gene.names <- x[, 1]
            read.count <- cbind(read.count, x[, 7])
        }
        # make column name
        v <- gsub('__on__hir.orig.txt|__on__hir.other.txt|__on__hir.common.txt', '', sort(count.files))
        v <- sapply(strsplit(v, '/'), function(x) {x[[8]]})
        v <- gsub('BG', 'BG_', v)
        v <- gsub('FRbr', 'KT5_', v)
        v <- gsub('FRsl', 'KT2_', v)
        v <- gsub('FRsh', 'SH_', v)
        v <- gsub('G', 'GX_', v)
        v <- gsub('IR', 'IR1_', v)
        v <- gsub('_H', '_H_', v)
        v <- gsub('_A', '_A_', v)
        v <- gsub('_F', '_F_', v)
        rownames(read.count) <- gene.names
        colnames(read.count) <- v
        as.matrix(read.count)
    }

    raw.count.h <- .get.counts(count.h.files)
    raw.count.a <- .get.counts(count.a.files)
    raw.count.u <- .get.counts(count.u.files)
    
    h.origin.ratio <- raw.count.h / (raw.count.h + raw.count.a)
    h.origin.ratio[, grep('_H_', colnames(h.origin.ratio))] <- 1
    h.origin.ratio[, grep('_A_', colnames(h.origin.ratio))] <- 0
    
    count.h <- raw.count.h
    count.a <- raw.count.a
    for (j in 1:ncol(raw.count.h)) {
        count.h[, j] <- apply(cbind(raw.count.h[, j], raw.count.u[, j] * h.origin.ratio[, j]), 1,
                              function(x) {sum(x, na.rm = TRUE)})
        count.a[, j] <- apply(cbind(raw.count.a[, j], raw.count.u[, j] * (1-h.origin.ratio[, j])), 1,
                              function(x) {sum(x, na.rm = TRUE)})
    }

    is.hir.lib <- grep('_H_', colnames(count.h))
    is.ama.lib <- grep('_A_', colnames(count.h))
    count.h[, is.hir.lib] <- raw.count.h[, is.hir.lib] + raw.count.u[, is.hir.lib] + raw.count.a[, is.hir.lib]
    count.a[, is.ama.lib] <- raw.count.h[, is.ama.lib] + raw.count.u[, is.ama.lib] + raw.count.a[, is.ama.lib]
    count.h[, is.ama.lib] <- 0
    count.a[, is.hir.lib] <- 0
    
    
    D$h.origin.ratio <- h.origin.ratio
    D$raw.count <- list(H = raw.count.h, U = raw.count.u, A = raw.count.a)
    D$count <- list(H = count.h, A = count.a)
    D
}


load.expression.zurich <- function(FL) {
    x <- read.table('~/research/cardamine3/zurich/CountQC_Cflex_all_on_Chir_2016-06-27--14-12-11/Count_QC/Count_QC-raw-count.txt', header = TRUE, sep = '\t')
    y <- x[, -1]
    rownames(y) <- x[, 1]
    
    raw.count.h <- y[, grep('orig', colnames(y))]
    raw.count.a <- y[, grep('other', colnames(y))]
    raw.count.u <- y[, grep('common', colnames(y))]
    
    v <- colnames(raw.count.h)
    v <- gsub('_Chir_other|_Chir_common|_Chir_orig', '', v)
    v <- gsub('FRbr', 'br_F_', v)
    v <- gsub('FRsl', 'sl_F_', v)
    v <- gsub('FRsh', 'sh_F_', v)
    v <- gsub('IR', 'ir_F_', v)
    
    colnames(raw.count.h) <- colnames(raw.count.a) <- colnames(raw.count.u) <- v
    
    h.origin.ratio <- raw.count.h / (raw.count.h + raw.count.a)
    h.origin.ratio[, grep('_H_', colnames(h.origin.ratio))] <- 1
    h.origin.ratio[, grep('_A_', colnames(h.origin.ratio))] <- 0
    
    count.h <- raw.count.h
    count.a <- raw.count.a
    for (j in 1:ncol(raw.count.h)) {
        count.h[, j] <- apply(cbind(raw.count.h[, j], raw.count.u[, j] * h.origin.ratio[, j]), 1,
                              function(x) {sum(x, na.rm = TRUE)})
        count.a[, j] <- apply(cbind(raw.count.a[, j], raw.count.u[, j] * (1-h.origin.ratio[, j])), 1,
                              function(x) {sum(x, na.rm = TRUE)})
    }

    FL$h.origin.ratio <- h.origin.ratio
    FL$raw.count <- list(H = raw.count.h, U = raw.count.u, A = raw.count.a)
    FL$count <- list(H = count.h, A = count.a)
    FL
    
}



calc.fpkm <- function(D) {
    .calc.fpkm <- function(y, len) {
        cpm  <- sweep(y, 2, 1e6 / colSums(y), '*')
        fpkm <- sweep(cpm, 1, 1e4 / len, '*')
        fpkm[is.na(fpkm)] <- 0
        fpkm
    }
    
    l <- read.table('~/research/cardamine3/data/genome/chir_cdna_length.tsv', header = T)
    len <- l[, 2]
    names(len) <- l[, 1]
    len <- len[rownames(D$count$H)]
    fpkm.h <- .calc.fpkm(D$count$H, len)
    fpkm.a <- .calc.fpkm(D$count$A, len)
    fpkm.total <- .calc.fpkm(D$raw.count$H + D$raw.count$A + D$raw.count$U, len)
    D$fpkm <- list(H = fpkm.h, A = fpkm.a)
    D$total.fpkm <- fpkm.total
    D
}



calc.exp.genes <- function(D) {
    D$exp.homeolog <- list(
        H = D$fpkm$H > 1,
        A = D$fpkm$A > 1
    )
    D
}




BindTairName <- function(dat) {
    if (!is.matrix(dat) && !is.data.frame(dat)) {
        dat <- as.matrix(dat)
        colnames(dat) <- "CARHR"
        rownames(dat) <- dat
    }
    rnm <- colnames(dat)
    tairid <- as.character((CARHR2TAIR[as.character(rownames(dat))]))
    genename <- as.character((CARHR2NAME[as.character(rownames(dat))]))
    genedesc <- as.character((CARHR2DESC[as.character(rownames(dat))]))
    tairid[tairid == "NULL"] <- ""
    genename[genename == "NULL"] <- ""
    genedesc[genedesc == "NULL"] <- ""
    dat <- data.frame(dat, TAIR = tairid, name = genename, description = genedesc, stringsAsFactors = F)
    dat[is.na(dat)] <- ""
    colnames(dat) <- c(rnm, "TAIR", "name", "description")
    dat
}




doEdgeR <- function(x, g, fdr.cutoff, figname) {
    x <- x[rowSums(cpm(x) > 1) > 0, ]
    dds <- DGEList(counts = x, group = g)
    dds <- calcNormFactors(dds)
    dsn <- model.matrix(~ g)
    dds <- estimateDisp(dds, dsn)
    fit <- glmFit(dds, dsn)
    lrt <- glmLRT(fit, coef = 2)
    df <- as.data.frame(topTags(lrt, n = nrow(x)))
    png(paste0('fig/', figname, '.png'), 500, 400)
    plotSmear(dds, de.tags  = rownames(df)[df$FDR < fdr.cutoff])
    dev.off()
    df <- data.frame(gene = as.character(rownames(df)), BindTairName(df))
    df <- data.frame(CarID = as.character(rownames(df)), ATID = df$TAIR, Symbols = df$name, Description = df$description, log2Ratio = df$logFC, FDR = df$FDR)
    df
}







.f.id <- 'data/genome/chi_m25.txt'
carhrf <- read.table(.f.id, header = TRUE, sep = "\t", fill = TRUE)
carhrf[, 6] <- gsub(" \\[.+\\]$", "", carhrf[,6])
carhrid <- unique(as.character(carhrf[, 1]))
CARHR2TAIR <- vector("list", length(carhrid))
CARHR2NAME <- vector("list", length(carhrid))
CARHR2DESC <- vector("list", length(carhrid))
names(CARHR2TAIR) <- names(CARHR2NAME) <- names(CARHR2DESC) <- carhrid
for (i in 1:length(CARHR2TAIR)) {
  CARHR2TAIR[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 2])
  CARHR2NAME[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 5])
  CARHR2DESC[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 6])
}
carhrf <- read.table(.f.id, header = TRUE, sep = "\t", fill = TRUE)
carhrf[, 6] <- gsub(" \\[.+\\]$", "", carhrf[,6])
tairid <- unique(as.character(carhrf[, 2]))
TAIR2CARHR <- vector("list", length(tairid))
names(TAIR2CARHR) <- tairid
for (i in 1:length(TAIR2CARHR)) {
  TAIR2CARHR[[i]] <- as.character(carhrf[carhrf[, 2] == tairid[i], 1])
}





.gettairs <- read.table('./zurich/car2tair_rbh_obh_20180423/chir2atid_rbh_obh_20160711/results/rbh_chir2tair/carhr2tair.txt', sep = '\t', stringsAsFactors = F)
.getcarhrs <- read.table('./zurich/car2tair_rbh_obh_20180423/chir2atid_rbh_obh_20160711/results/rbh_chir2tair/tair2carhr.txt', sep = '\t', stringsAsFactors = F)
.carhr2tarivec <- .gettairs[, 2]
names(.carhr2tarivec) <- .gettairs[, 1]
.tair2carhrvec <- .getcarhrs[, 2]
names(.tair2carhrvec) <- .getcarhrs[, 1]

getTAIRs <- function(x) {
    y <- .carhr2tarivec[x]
    y <- y[!is.na(y)]
    y <- paste(y, collapse = ',')
    y <- unlist(strsplit(y, ','))
    y
}

getCARHRs <- function(x) {
    y <- .tair2carhrvec[x]
    y <- y[!is.na(y)]
    y <- paste(y, collapse = ',')
    y <- unlist(strsplit(y, ','))
    y
}


CARHR2TAIR <- .gettairs[, 2]
names(CARHR2TAIR) <- .gettairs[, 1]
TAIR2CARHR <- .getcarhrs[, 2]
names(TAIR2CARHR) <- .getcarhrs[, 1]


.x <- read.table('./data/gene_aliases_20140331.txt.adjusted', sep = '\t')
TAIR2NAME <- .x[, 2]
TAIR2DESC <- .x[, 3]
names(TAIR2NAME) <- names(TAIR2DESC) <- .x[, 1]

CARHR2SYMBOL <- TAIR2NAME[CARHR2TAIR]
CARHR2DESC   <- TAIR2DESC[CARHR2TAIR]
names(CARHR2SYMBOL) <- names(CARHR2DESC) <- names(CARHR2TAIR)

for (.i in names(CARHR2DESC)) {
    if (!is.na(CARHR2DESC[[.i]]) && CARHR2DESC[[.i]] == 'NA') {
        CARHR2DESC[[.i]] <- NA
    }
}
for (.i in names(CARHR2SYMBOL)) {
    if (!is.na(CARHR2SYMBOL[[.i]]) && CARHR2SYMBOL[[.i]] == 'NA') {
        CARHR2SYMBOL[[.i]] <- NA
    }
}
CARHR2NAME <- CARHR2SYMBOL



#carhrf <- read.table(.f.id, header = TRUE, sep = "\t", fill = TRUE)
#carhrf[, 6] <- gsub(" \\[.+\\]$", "", carhrf[,6])
#carhrid <- unique(as.character(carhrf[, 1]))
#carhrID2tairID <- vector("list", length(carhrid))
#carhrID2carhrName <- vector("list", length(carhrid))
#carhrID2carhrDesc <- vector("list", length(carhrid))
#names(carhrID2tairID) <- names(carhrID2carhrName) <- names(carhrID2carhrDesc) <- carhrid
#for (i in 1:length(carhrID2tairID)) {
#  carhrID2tairID[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 2])
#  carhrID2carhrName[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 5])
#  carhrID2carhrDesc[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 6])
#}
#carhrf <- read.table(.f.id, header = TRUE, sep = "\t", fill = TRUE)
#carhrf[, 6] <- gsub(" \\[.+\\]$", "", carhrf[,6])
#tairid <- unique(as.character(carhrf[, 2]))
#tairID2carhrID <- vector("list", length(tairid))
#names(tairID2carhrID) <- tairid
#for (i in 1:length(tairID2carhrID)) {
#  tairID2carhrID[[i]] <- as.character(carhrf[carhrf[, 2] == tairid[i], 1])
#}


geneID2GO <- readMappings('data/AT_GO_all.txt')
#geneID2GO <- readMappings('topgo/AT_GO_all_subset.txt2')
#geneID2GO <- readMappings('topgo/AT_GO_all_subset.txt.20180605')
#.x <- read.table('topgo/chi_m25.txt.adjusted')
#carhrID2tairID <- as.character(.x[, 2])
#names(carhrID2tairID) <- as.character(.x[, 1])
#tairID2carhrID <- as.character(.x[, 1])
#names(tairID2carhrID) <- as.character(.x[, 2])
#rm(.x)

.annotatedTairs <- read.table('data/go.num.txt', sep = '\t')
annotatedTairs <- .annotatedTairs[, 2]
names(annotatedTairs) <- .annotatedTairs[, 1]

doTopGO <- function(allgene, siggene) {
    allgene <- as.character(allgene)
    siggene <- as.character(siggene)
    allgene.tair <- getTAIRs(allgene)
    siggene.tair <- getTAIRs(siggene)
    allgene.f <- as.numeric(allgene.tair %in% siggene.tair)
    names(allgene.f) <- allgene.tair
    topgo.obj <- new('topGOdata', ontology = 'BP', allGenes = allgene.f, geneSel = function(x) x > 0,
                     #annot = annFUN.org, mapping = 'org.At.tair.db')
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    #classic_fisher <- runTest(topgo.obj, algorithm = "classic", statistic = "fisher")
    elim_fisher    <- runTest(topgo.obj, algorithm = "elim", statistic = "fisher")
    #weight01_fisher  <- runTest(topgo.obj, algorithm = "weight01", statistic = "fisher")
    topgotable  <- GenTable(topgo.obj, numChar = 200,
                            # classicFisher = classic_fisher,
                            elimFisher = elim_fisher,
                            #weight01 = weight01_fisher,
                            orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = 600)
    topgotable <- topgotable[topgotable[, 3] > 10, ]
    topgotable <- topgotable[topgotable[, 3] < 500, ]
    topgotable <- data.frame(topgotable, TAIR10annotated = annotatedTairs[topgotable[, 1]])
    topgotable
}



barplotExp <- function(interesth, fig.prefix = 'test') {
    dir.create(paste0('./fig/barplots/', fig.prefix), recursive = TRUE)
    cname <- colnames(D$fpkm$H)
    tag <- paste0(fig.prefix, '/KT5_IR1')

    ir.hir <- D$fpkm$H[, grep(paste0('IR1_H_0502'), cname)][interesth, ]
    ir.ama <- D$fpkm$A[, grep(paste0('IR1_A_0502'), cname)][interesth, ]
    br.hir <- D$fpkm$H[, grep(paste0('KT5_H_0502'), cname)][interesth, ]
    br.ama <- D$fpkm$A[, grep(paste0('KT5_A_0502'), cname)][interesth, ]
    ir.flxH <- D$fpkm$H[, grep(paste0('IR1_F_0502'), cname)][interesth, ]
    ir.flxA <- D$fpkm$A[, grep(paste0('IR1_F_0502'), cname)][interesth, ]
    br.flxH <- D$fpkm$H[, grep(paste0('KT5_F_0502'), cname)][interesth, ]
    br.flxA <- D$fpkm$A[, grep(paste0('KT5_F_0502'), cname)][interesth, ]

    .meltmat <- function(y, species, site) {
        df <- NULL
        if (ncol(y) > 0) {
            y <- log10(y + 1)
            df <- data.frame(gene = rownames(y), mean = apply(y, 1, mean), sd = apply(y, 1, sd), species = species, site = site)
        }
        df
    }
 
    df <- rbind(data.frame(.meltmat(ir.hir, species = 'C. hirsuta', 'IR1')), 
                data.frame(.meltmat(ir.ama, species = 'C. amara', 'IR1')),
                data.frame(.meltmat(br.hir, species = 'C. hirsuta', 'KT5')),
                data.frame(.meltmat(br.ama, species = 'C. amara', 'KT5')),
                data.frame(.meltmat(ir.flxH, species = 'C. flexuosa (H)', 'IR1')),
                data.frame(.meltmat(ir.flxA, species = 'C. flexuosa (A)', 'IR1')),
                data.frame(.meltmat(br.flxH, species = 'C. flexuosa (H)', 'KT5')),
                data.frame(.meltmat(br.flxA, species = 'C. flexuosa (A)', 'KT5')))
    df$species <- factor(df$species, levels = c('C. hirsuta', 'C. flexuosa (H)', 'C. flexuosa (A)', 'C. amara'))
    df$site <- factor(df$site, levels = c('IR1', 'KT5'))
    for (gn in unique(df$gene)) {
            gn.tair <- CARHR2TAIR[[gn]]
            gn.name <- CARHR2NAME[[gn]]
            gn.fullname <- gn
            if (!is.null(gn.tair) && !is.null(gn.name)) {
                gn.fullname <- paste0(gn, ' (', gn.tair, '/', gn.name, ')')
            } else if (!is.null(gn.tair) && is.null(gn.name)) {
                gn.fullname <- paste0(gn, ' (', gn.tair, ')')
            }
            dfg <- df[df$gene == gn, ]
            g <- ggplot(dfg, aes(x = species, y = mean))
            g <- g + geom_bar(stat = 'identity') + ggtitle(gn.fullname)
            g <- g + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3))
            g <- g + theme(axis.text.x=element_text(angle = 30, vjust = 0.5, hjust = 0.5))
            g <- g + facet_grid(. ~ site, scales = "free_x") + ylab('log10(FPKM)') + xlab('') 
            
            png(paste0('./fig/barplots/', fig.prefix, '/', gn, '.png'), 400, 300)
            print(g)
            dev.off()
    }
 
}


ak.deg <- function(D) {
    fDEG <- aDEG <- hDEG <- NULL
    fDEGat <- aDEGat <- hDEGat <- NULL
    fdr.cutoff <- 0.05
    
    # DEG for amara
    d <- D$count$H + D$count$A
    a <- d[, grep('_A_', colnames(d))]
    f <- d[, grep('_F_', colnames(d))]
    h <- d[, grep('_H_', colnames(d))]
    fH <- D$count$H[, grep('_F_', colnames(D$count$H))]
    fA <- D$count$A[, grep('_F_', colnames(D$count$A))]

    # DEG for flexuosa
    for (dt in c('0418', '0502', '0516')) {
        fdeg <- fdegat <- NULL
        .f <- f[, grep(dt, colnames(f))]
        # ir x sl
        .f1 <- .f[, grep('IR1|KT2', colnames(.f))]
        .l1 <- sapply(strsplit(colnames(.f1), '_'), function(x) {x[[1]]})
        fdeg[['KT2_IR1']] <- doEdgeR(.f1, .l1, fdr.cutoff, paste0('flexall_KT2_vs_IR1_', dt))
        fdegat[['KT2_IR1']] <- fdeg[['KT2_IR1']][(fdeg[['KT2_IR1']][, 1] %in% names(CARHR2TAIR)), ]
        # ir x br
        .f2 <- .f[, grep('IR1|KT5', colnames(.f))]
        .l2 <- sapply(strsplit(colnames(.f2), '_'), function(x) {x[[1]]})
        fdeg[['KT5_IR1']] <- doEdgeR(.f2, .l2, fdr.cutoff, paste0('flexall_KT5_vs_IR1_', dt))
        fdegat[['KT5_IR1']] <- fdeg[['KT5_IR1']][(fdeg[['KT5_IR1']][, 1] %in% names(CARHR2TAIR)), ]
        # sl x br
        .f3 <- .f[, grep('KT2|KT5', colnames(.f))]
        .l3 <- sapply(strsplit(colnames(.f3), '_'), function(x) {x[[1]]})
        fdeg[['KT5_KT2']] <- doEdgeR(.f3, .l3, fdr.cutoff, paste0('flexall_KT5_KT2_', dt))
        fdegat[['KT5_KT2']] <- fdeg[['KT5_KT2']][(fdeg[['KT5_KT2']][, 1] %in% names(CARHR2TAIR)), ]
        fDEG[[dt]] <- fdeg
        fDEGat[[dt]] <- fdegat
    }
    fDEG2 <- fDEG2at <- NULL
    for (nm1 in names(fDEG)) {
        for (nm2 in names(fDEG[[nm1]])) {
            print(paste0(nm1, ' ', nm2, ' ', sum(fDEG[[nm1]][[nm2]]$FDR < fdr.cutoff),
                         '  AT:', sum(fDEGat[[nm1]][[nm2]]$FDR < fdr.cutoff)))
            fDEG2[[paste0(nm1, '_', nm2)]] <- fDEG[[nm1]][[nm2]]
            fDEG2at[[paste0(nm1, '_', nm2)]] <- fDEGat[[nm1]][[nm2]]
        }
    }
    WriteExcel(fDEG2, file = 'result/DEG.xlsx')
    WriteExcel(fDEG2at, file = 'result/DEG.AT.xlsx')
    
    
    # topGO
    fDEGGO <- fDEG
    fDEGGO2 <- NULL
    for (dt in names(fDEG)) {
        for (cp in names(fDEG[[dt]])) {
            degtable <- fDEG[[dt]][[cp]]
            allgene <- degtable$CarID
            siggene <- degtable$CarID[degtable$FDR < fdr.cutoff]
            fDEGGO[[dt]][[cp]] <- doTopGO(allgene, siggene)
            fDEGGO2[[paste0(dt, '_', cp)]] <- fDEGGO[[dt]][[cp]]
        }
    }
    WriteExcel(fDEGGO2, file = 'result/DEGGO.xlsx')
    
 
    used.go.id <- unique(c(unlist(sapply(fDEGGO[[1]], function(x){x[, 1]})),
                           unlist(sapply(fDEGGO[[2]], function(x){x[, 1]})),
                           unlist(sapply(fDEGGO[[3]], function(x){x[, 1]}))))
    gotable   <- matrix(NA, ncol = 9, nrow = length(used.go.id))
    gotable.b <- matrix(0, ncol = 9, nrow = length(used.go.id))
    rownames(gotable) <- rownames(gotable.b) <- used.go.id
    colnames(gotable) <- colnames(gotable.b) <- paste0(rep(names(fDEG), each = 3), '_', rep(names(fDEG[[1]]), times = 3))
    term2desc <- rep('', length = length(used.go.id))
    names(term2desc) <- used.go.id
    for (dt in names(fDEGGO)) {
        for (cp in names(fDEGGO[[dt]])) {
            toptable <- fDEGGO[[dt]][[cp]]
            gotable[toptable[, 1], paste0(dt, '_', cp)] <- as.numeric(toptable[, 6])
            gotable.b[toptable[, 1], paste0(dt, '_', cp)] <- as.numeric(as.numeric(toptable[, 6]) < 0.05)
            term2desc[toptable[, 1]] <- toptable[, 2]
        }
    }
    gotable.b <- data.frame(go = rownames(gotable.b), desc = term2desc[rownames(gotable.b)],
                            gotable.b, num.of.case = rowSums(gotable.b))
    gotable.b <- gotable.b[order(- gotable.b[['num.of.case']]), ]
    write.table(gotable.b, file = 'result/DEG_GOTable.summary.xls', sep = '\t')
    
 

   
    # GO:0009414
    .getGOGene <- function(gonum) {
        x.tair <- read.table(paste0('topgo/GO', gonum, '.txt'), sep = '\t')
        x.tair <- as.character(x.tair[, 1])
        x.carhr <- getCARHRs(x.tair)
        x.carhr <- unique(as.character(x.carhr[!is.na(x.carhr)]))
        x.carhr
    }
    watdep.carhr <- .getGOGene('0009414')
    acqres.carhr <- .getGOGene('0009627')
    rscold.carhr <- .getGOGene('0009409')
    rsheat.carhr <- .getGOGene('0009408')
    
    vk <- NULL
    for (dt in names(fDEG)) {
        v <- list(KT2_IR1 = intersect(fDEG[[dt]][['KT2_IR1']]$CarID[fDEG[[dt]][['KT2_IR1']]$FDR < fdr.cutoff], watdep.carhr),
                  KT5_IR1 = intersect(fDEG[[dt]][['KT5_IR1']]$CarID[fDEG[[dt]][['KT5_IR1']]$FDR < fdr.cutoff], watdep.carhr),
                  KT5_KT2 = intersect(fDEG[[dt]][['KT5_KT2']]$CarID[fDEG[[dt]][['KT5_KT2']]$FDR < fdr.cutoff], watdep.carhr))
		vuq <- unique(unlist(v))
		vdf <- matrix(FALSE, ncol = 3, nrow = length(vuq))
		rownames(vdf) <- vuq
		colnames(vdf) <- c('KT2_IR1', 'KT5_IR1', 'KT5_KT2')
		vdf[v[['KT2_IR1']], 'KT2_IR1'] <- TRUE
		vdf[v[['KT5_IR1']], 'KT5_IR1'] <- TRUE
		vdf[v[['KT5_KT2']], 'KT5_KT2'] <- TRUE
		fit <- euler(vdf)
        png(paste0('fig/vennDEG-date-', dt, '.png'), 500, 500)
		print(plot(fit, quantities = F, lty = 1, shape = 'ellipse', labels = F,
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
        png(paste0('fig/vennDEG-date-', dt, '.labeled.png'), 500, 500)
		print(plot(fit, quantities = list(cex = 1.9), lty = 1, shape = 'ellipse',
			 labels = list(cex = 1.4), 
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
    }
	
	
    vs <- NULL
    for (cp in names(fDEG[[1]])) {
        v <- list(d0418 = intersect(fDEG[['0418']][[cp]]$CarID[fDEG[['0418']][[cp]]$FDR < fdr.cutoff], watdep.carhr),
                  d0502 = intersect(fDEG[['0502']][[cp]]$CarID[fDEG[['0502']][[cp]]$FDR < fdr.cutoff], watdep.carhr),
                  d0516 = intersect(fDEG[['0516']][[cp]]$CarID[fDEG[['0516']][[cp]]$FDR < fdr.cutoff], watdep.carhr))
		vuq <- unique(unlist(v))
		vdf <- matrix(FALSE, ncol = 3, nrow = length(vuq))
		rownames(vdf) <- vuq
		colnames(vdf) <- c('0418', '0502', '0516')
		vdf[v[['d0418']], '0418'] <- TRUE
		vdf[v[['d0502']], '0502'] <- TRUE
		vdf[v[['d0516']], '0516'] <- TRUE
		fit <- euler(vdf)
        png(paste0('fig/vennDEG-site-', cp, '.png'), 500, 500)
		print(plot(fit, quantities = F, lty = 1, shape = 'ellipse', labels = F,
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
        png(paste0('fig/vennDEG-site-', cp, '.labeled.png'), 500, 500)
		print(plot(fit, quantities = list(cex = 1.9), lty = 1, shape = 'ellipse',
			 labels = list(cex = 1.4), 
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
    }
    
    interestd <- NULL
    for (nm1 in names(fDEG)) {
        for (nm2 in names(fDEG[[nm1]])) {
            interestd <- c(interestd, intersect(fDEG[[nm1]][[nm2]]$CarID[fDEG[[nm1]][[nm2]]$FDR < fdr.cutoff], watdep.carhr))
        }
    }
    interestd <- unique(interestd)
    barplotExp(interestd, 'flexuosa_DEG_waterdeprivation')
    WriteExcel(list(Sheet1 = BindTairName(interestd)), './result/flexuosa_DEG_waterdeprivation.xlsx')


    # 0502
    .ir1.hir   <- h[, grep('IR1_H_0502', colnames(h))]
    .ir1.flex  <- f[, grep('IR1_F_0502', colnames(f))]
    .ir1.flexa <- fA[, grep('IR1_F_0502', colnames(fA))]
    .ir1.flexh  <- fH[, grep('IR1_F_0502', colnames(fH))]
    .kt5.flex  <- f[, grep('KT5_F_0502', colnames(f))]
    .kt5.flexa <- fA[, grep('KT5_F_0502', colnames(fA))]
    .kt5.flexh <- fH[, grep('KT5_F_0502', colnames(fH))]
    .kt5.ama   <- a[, grep('KT5_A_0502', colnames(a))]
    
    # amara vs hirsuta // KT5 vs IR1, 0502
    zpDEG <- doEdgeR(cbind(.ir1.hir, .kt5.ama), c(rep('IRI-H', 3), rep('KT5-A', 3)), fdr.cutoff, 'IR1-H_KT5-A')
    zpdeg.allgene <- zpDEG$CarID
    zpdeg.siggene <- zpDEG$CarID[zpDEG$FDR < fdr.cutoff]
    zpdeg.siggeneup <- zpDEG$CarID[zpDEG$FDR < fdr.cutoff & zpDEG$log2Ratio > 0]
    zpdeg.siggenedw <- zpDEG$CarID[zpDEG$FDR < fdr.cutoff & zpDEG$log2Ratio < 0]
    zpDEGGO <- doTopGO(zpdeg.allgene, zpdeg.siggene)
    zpDEGGOup <- doTopGO(zpdeg.allgene, zpdeg.siggeneup)
    zpDEGGOdw <- doTopGO(zpdeg.allgene, zpdeg.siggenedw)
    WriteExcel(list(DEG = zpDEG, DEGAT = zpDEG[(zpDEG$CarID %in% names(CARHR2TAIR)), ],
                    GO = zpDEGGO, GO_A_highexp = zpDEGGOup, GO_H_highexp = zpDEGGOdw),
                    file = 'result/amara0502_hirsuta0502_DEG.xlsx')
    zpDEG.watdep <- intersect(as.character(zpdeg.siggene), as.character(watdep.carhr))
    barplotExp(zpDEG.watdep, 'amara0502_vs_hirsuta0502_waterdeprivation')
    zpDEG.rscold <- intersect(as.character(zpdeg.siggene), as.character(rscold.carhr))
    barplotExp(zpDEG.rscold, 'amara0502_vs_hirsuta0502_responsetocold')
    zpDEG.rsheat <- intersect(as.character(zpdeg.siggene), as.character(rsheat.carhr))
    barplotExp(zpDEG.rsheat, 'amara0502_vs_hirsuta0502_responsetoheat')

    # flexuoasa vs hirsuta // IR1, 0502
    zhDEG <- doEdgeR(cbind(.ir1.flex, .ir1.hir), c(rep('IR1-F', 3), rep('IR1-H', 3)), fdr.cutoff, 'IR1-F_IR1-H')
    zhdeg.allgene <- zhDEG$CarID
    zhdeg.siggene <- zhDEG$CarID[zhDEG$FDR < fdr.cutoff]
    zhdeg.siggeneup <- zhDEG$CarID[zhDEG$FDR < fdr.cutoff & zhDEG$log2Ratio > 0]
    zhdeg.siggenedw <- zhDEG$CarID[zhDEG$FDR < fdr.cutoff & zhDEG$log2Ratio < 0]
    zhDEGGO <- doTopGO(zhdeg.allgene, zhdeg.siggene)
    zhDEGGOup <- doTopGO(zhdeg.allgene, zhdeg.siggeneup)
    zhDEGGOdw <- doTopGO(zhdeg.allgene, zhdeg.siggenedw)
    WriteExcel(list(DEG = zhDEG, DEGAT = zhDEG[(zhDEG$CarID %in% names(CARHR2TAIR)), ],
                    GO = zhDEGGO, GO_H_highexp = zhDEGGOup, GO_F_highexp = zhDEGGOdw),
                    file = 'result/flexuosa0502IR1_hirsuta0502IR1_DEG.xlsx')


    # flexuoasa vs amara // KT5, 0502
    zaDEG <- doEdgeR(cbind(.kt5.flex, .kt5.ama), c(rep('KT5-F', 3), rep('KT5-A', 3)), fdr.cutoff, 'KT5-F_KT5-A')
    zadeg.allgene <- zaDEG$CarID
    zadeg.siggene <- zaDEG$CarID[zaDEG$FDR < fdr.cutoff]
    zadeg.siggeneup <- zaDEG$CarID[zaDEG$FDR < fdr.cutoff & zaDEG$log2Ratio > 0]
    zadeg.siggenedw <- zaDEG$CarID[zaDEG$FDR < fdr.cutoff & zaDEG$log2Ratio < 0]
    zaDEGGO <- doTopGO(zadeg.allgene, zadeg.siggene)
    zaDEGGOup <- doTopGO(zadeg.allgene, zadeg.siggeneup)
    zaDEGGOdw <- doTopGO(zadeg.allgene, zadeg.siggenedw)
    WriteExcel(list(DEG = zaDEG, DEGAT = zaDEG[(zaDEG$CarID %in% names(CARHR2TAIR)), ],
                    GO = zaDEGGO, GO_H_highexp = zaDEGGOup, GO_A_highexp = zaDEGGOdw),
                    file = 'result/flexuosa0502KT5_amara0502KT5_DEG.xlsx')

    
    ziDEG <- fDEG[['0502']][['KT5_IR1']]
    zideg.allgene <- ziDEG$CarID
    zideg.siggene <- ziDEG$CarID[ziDEG$FDR < fdr.cutoff]
    zideg.siggeneup <- ziDEG$CarID[ziDEG$FDR < fdr.cutoff & ziDEG$log2Ratio > 0]
    zideg.siggenedw <- ziDEG$CarID[ziDEG$FDR < fdr.cutoff & ziDEG$log2Ratio < 0]
    ziDEGGO <- doTopGO(zideg.allgene, zideg.siggene)
    ziDEGGOup <- doTopGO(zideg.allgene, zideg.siggeneup)
    ziDEGGOdw <- doTopGO(zideg.allgene, zideg.siggenedw)
    WriteExcel(list(DEG = ziDEG, DEGAT = ziDEG[(ziDEG$CarID %in% names(CARHR2TAIR)), ],
                    GO = ziDEGGO, GO_KT5_highexp = ziDEGGOup, GO_IR1_highexp = ziDEGGOdw),
                    file = 'result/flexuosa0502KT5_flexuosa0502IR1_DEG.xlsx')
    
    zideg.siggene.watdep <- intersect(as.character(zideg.siggene), as.character(watdep.carhr))
    zpdeg.siggene.watdep <- intersect(as.character(zpdeg.siggene), as.character(watdep.carhr))
    barplotExp(intersect(zideg.siggene.watdep, zpdeg.siggene.watdep), 'parent0502DEG_flexuosa0502DEG_shared_waterdeprivation')
    
    
 
    # flexuoasa H vs hirsuta // IR1, 0502
    fhDEG <- doEdgeR(cbind(.ir1.flexh, .ir1.hir), c(rep('IR1-Fh', 3), rep('IR1-H', 3)), fdr.cutoff, 'IR1-Fh_IR1-H')
    fhdeg.allgene <- fhDEG$CarID
    fhdeg.siggene <- fhDEG$CarID[fhDEG$FDR < fdr.cutoff]
    fhdeg.siggeneup <- fhDEG$CarID[fhDEG$FDR < fdr.cutoff & fhDEG$log2Ratio > 0]
    fhdeg.siggenedw <- fhDEG$CarID[fhDEG$FDR < fdr.cutoff & fhDEG$log2Ratio < 0]
    fhDEGGO <- doTopGO(fhdeg.allgene, fhdeg.siggene)
    fhDEGGOup <- doTopGO(fhdeg.allgene, fhdeg.siggeneup)
    fhDEGGOdw <- doTopGO(fhdeg.allgene, fhdeg.siggenedw)
    WriteExcel(list(DEG = fhDEG, DEGAT = fhDEG[(fhDEG$CarID %in% names(CARHR2TAIR)), ],
                    GO = fhDEGGO, GO_H_highexp = fhDEGGOup, GO_Fh_highexp = fhDEGGOdw),
                    file = 'result/flexuosaH0502IR1_hirsuta0502IR1_DEG.xlsx')


    # flexuoasa A vs amara // KT5, 0502
    faDEG <- doEdgeR(cbind(.kt5.flexa, .kt5.ama), c(rep('KT5-Fa', 3), rep('KT5-A', 3)), fdr.cutoff, 'KT5-Fa_KT5-A')
    fadeg.allgene <- faDEG$CarID
    fadeg.siggene <- faDEG$CarID[faDEG$FDR < fdr.cutoff]
    fadeg.siggeneup <- faDEG$CarID[faDEG$FDR < fdr.cutoff & faDEG$log2Ratio > 0]
    fadeg.siggenedw <- faDEG$CarID[faDEG$FDR < fdr.cutoff & faDEG$log2Ratio < 0]
    faDEGGO <- doTopGO(fadeg.allgene, fadeg.siggene)
    faDEGGOup <- doTopGO(fadeg.allgene, fadeg.siggeneup)
    faDEGGOdw <- doTopGO(fadeg.allgene, fadeg.siggenedw)
    WriteExcel(list(DEG = faDEG, DEGAT = faDEG[(faDEG$CarID %in% names(CARHR2TAIR)), ],
                    GO = faDEGGO, GO_A_highexp = faDEGGOup, GO_Fa_highexp = faDEGGOdw),
                    file = 'result/flexuosaA0502KT5_amara0502KT5_DEG.xlsx')

   
    # ratio-changes
    fRT <- NULL
    dh <- D$fpkm$H
    da <- D$fpkm$A
    fh <- dh[, grep('_F_', colnames(dh))]
    fa <- da[, grep('_F_', colnames(da))]
    for (dt in c('0418', '0502', '0516')) {
        frt <- NULL
        .fh <- fh[, grep(dt, colnames(fh))]
        .fa <- fa[, grep(dt, colnames(fa))]
        
        dat <- cbind(.fh[, grep('KT2', colnames(.fh))], .fa[, grep('KT2', colnames(.fa))],
                     .fh[, grep('IR1', colnames(.fh))], .fa[, grep('IR1', colnames(.fa))])
        tmpdir <- paste0('homeoroq_tmp/', dt, '_KT2_IR1')
        if (!dir.exists(tmpdir)) {
            dir.create(tmpdir, recursive = TRUE)
        }
        dat <- dat[rowSums(dat) > 0, ]
        write.table(dat, col.names = T, row.names = T, sep = '\t', quote = F, file = paste0(tmpdir, '/data.xls'))
        
        dat <- cbind(.fh[, grep('KT5', colnames(.fh))], .fa[, grep('KT5', colnames(.fa))],
                     .fh[, grep('IR1', colnames(.fh))], .fa[, grep('IR1', colnames(.fa))])
        tmpdir <- paste0('homeoroq_tmp/', dt, '_KT5_IR1')
        if (!dir.exists(tmpdir)) {
            dir.create(tmpdir, recursive = TRUE)
        }
        dat <- dat[rowSums(dat) > 0, ]
        write.table(dat, col.names = T, row.names = T, sep = '\t', quote = F, file = paste0(tmpdir, '/data.xls'))
        
        dat <- cbind(.fh[, grep('KT5', colnames(.fh))], .fa[, grep('KT5', colnames(.fa))],
                     .fh[, grep('KT2', colnames(.fh))], .fa[, grep('KT2', colnames(.fa))])
        tmpdir <- paste0('homeoroq_tmp/', dt, '_KT5_KT2')
        if (!dir.exists(tmpdir)) {
            dir.create(tmpdir, recursive = TRUE)
        }
        dat <- dat[rowSums(dat) > 0, ]
        write.table(dat, col.names = T, row.names = T, sep = '\t', quote = F, file = paste0(tmpdir, '/data.xls'))
    }
    # goto homeoroq_tmp directory and run the shell.
    # run.calcpval.1.sh
    # run.calcpval.2.sh
    fRT <- fRTtb <- fRTat <- NULL
    for (dname in list.dirs('homeoroq_tmp')) {
        if (dname == 'homeoroq_tmp') {next}
        x <- read.table(paste0(dname, '/result.txt'), header = TRUE)
        x <- data.frame(x, sigGene = as.numeric(x$ratioSD < 0.3 & x$padj < fdr.cutoff))
        y <- x[, c('gene', 'padj', 'ratioSD', 'sigGene')]
        rownames(x) <- x[, 1]
        rownames(y) <- y[, 1]
        x <- BindTairName(x)
        y <- BindTairName(y)
        print(paste0(dname, ' ', sum(x$sigGene),
                     '  AT:', sum(x[x$gene %in% names(CARHR2TAIR), ]$sigGene)))
        x <- x[order(-x$sigGene), ]
        y <- y[order(-y$sigGene), ]
        fRT[[strsplit(dname, '/')[[1]][2]]] <- x
        fRTtb[[strsplit(dname, '/')[[1]][2]]] <- y
        fRTat[[strsplit(dname, '/')[[1]][2]]] <- y[y$gene %in% names(CARHR2TAIR), ]
    }
    WriteExcel(fRTtb, file = 'result/RatioDiffGene.xlsx')
    WriteExcel(fRTat, file = 'result/RatioDiffGeneAT.xlsx')

    # topGO
    fRTGO <- fRT
    for (dn in names(fRT)) {
        degtable <- fRT[[dn]]
        allgene <- as.character(degtable$gene)
        siggene <- as.character(degtable$gene[degtable$sigGene > 0])
        fRTGO[[dn]] <- doTopGO(allgene, siggene)
    }
    WriteExcel(fRTGO, file = 'result/RatioDiffGeneGO.xlsx')
    
    used.go.id <- unique(c(unlist(sapply(fRTGO, function(x){x[, 1]}))))
    gotable   <- matrix(NA, ncol = 9, nrow = length(used.go.id))
    gotable.b <- matrix(0, ncol = 9, nrow = length(used.go.id))
    rownames(gotable) <- rownames(gotable.b) <- used.go.id
    colnames(gotable) <- colnames(gotable.b) <- paste0(rep(names(fDEG), each = 3), '_', rep(names(fDEG[[1]]), times = 3))
    term2desc <- rep('', length = length(used.go.id))
    names(term2desc) <- used.go.id
    for (dt in names(fRTGO)) {
        toptable <- fRTGO[[dt]]
        gotable[toptable[, 1], dt] <- as.numeric(toptable[, 6])
        gotable.b[toptable[, 1], dt] <- as.numeric(as.numeric(toptable[, 6]) < 0.05)
        term2desc[toptable[, 1]] <- toptable[, 2]
    }
    gotable.b <- data.frame(go = rownames(gotable.b), desc = term2desc[rownames(gotable.b)],
                            gotable.b, num.of.case = rowSums(gotable.b))
    gotable.b <- gotable.b[order(- gotable.b[['num.of.case']]), ]
    write.table(gotable.b, file = 'result/RatioDiffGene_GOTable.summary.xls', sep = '\t')
    
    
    for (dt in c('0418', '0502', '0516')) {
        frt <- fRT[grep(dt, names(fRT))]
        frt2 <- NULL
        for (dm in names(frt)) {
            frt2[[gsub(paste0(dt, '_'), '', dm)]] <- intersect(frt[[dm]]$gene[frt[[dm]]$sigGene == 1], watdep.carhr)
        }
		vuq <- unique(unlist(frt2))
		vdf <- matrix(FALSE, ncol = 3, nrow = length(vuq))
		rownames(vdf) <- vuq
		colnames(vdf) <- c('KT2_IR1', 'KT5_IR1', 'KT5_KT2')
		vdf[frt2[['KT2_IR1']], 'KT2_IR1'] <- TRUE
		vdf[frt2[['KT5_IR1']], 'KT5_IR1'] <- TRUE
		vdf[frt2[['KT5_KT2']], 'KT5_KT2'] <- TRUE
		fit <- euler(vdf)
        png(paste0('fig/vennDiffRatio-date-', dt, '.png'), 500, 500)
		print(plot(fit, quantities = F, lty = 1, shape = 'ellipse', labels = F,
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
        png(paste0('fig/vennDiffRatio-date-', dt, '.labeled.png'), 500, 500)
		print(plot(fit, quantities = list(cex = 1.9), lty = 1, shape = 'ellipse',
			 labels = list(cex = 1.4), 
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
 	}


    for (dt in c('KT2_IR1', 'KT5_KT2', 'KT5_IR1')) {
        frt <- fRT[grep(dt, names(fRT))]
        frt2 <- NULL
        for (dm in names(frt)) {
            frt2[[gsub(paste0('_', dt), '', dm)]] <- intersect(frt[[dm]]$gene[frt[[dm]]$sigGene == 1], watdep.carhr)
        }
		vuq <- unique(unlist(frt2))
		vdf <- matrix(FALSE, ncol = 3, nrow = length(vuq))
		rownames(vdf) <- vuq
		colnames(vdf) <- c('0418', '0502', '0516')
		vdf[frt2[['0418']], '0418'] <- TRUE
		vdf[frt2[['0502']], '0502'] <- TRUE
		vdf[frt2[['0516']], '0516'] <- TRUE
		fit <- euler(vdf)
        png(paste0('fig/vennDiffRatio-site-', dt, '.png'), 500, 500)
		print(plot(fit, quantities = F, lty = 1, shape = 'ellipse', labels = F,
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
        png(paste0('fig/vennDiffRatio-site-', dt, '.labeled.png'), 500, 500)
		print(plot(fit, quantities = list(cex = 1.9), lty = 1, shape = 'ecllipse',
			 labels = list(cex = 1.4), 
			 fill = list(fill = c("purple", "darkorange", "forestgreen"), alpha = 0.3)))
		dev.off()
	}


    # barplots
    interesth <- interesth0502 <- NULL
    for (dt in c('KT2_IR1', 'KT5_KT2', 'KT2_IR1')) {
        frt <- fRT[grep(dt, names(fRT))]
        for (dm in names(frt)) {
            interesth <- c(interesth, intersect(frt[[dm]]$gene[frt[[dm]]$sigGene == 1], watdep.carhr))
        }
    }
    interesth <- unique(interesth)
    WriteExcel(list(Sheet1 = BindTairName(interesth)), './result/flexuosa_sigDiffRatio_waterdeprivation.xlsx')
    barplotExp(interesth, 'flexuosa_sigDiffRatio_waterdeprivation')
    
    venn(list(interesth = interesth, interestd = interestd))
    barplotExp(intersect(interesth, interestd), 'flexuosa_DEGandSigDiffRatio_waterdeprivation')
    WriteExcel(list(Sheet1 = BindTairName(intersect(interesth, interestd))),
                     './result/flexuosa_DEGandSigDiffRatio_waterdeprivation.xlsx')

    
    venn(list(flexDEG = zideg.siggene.watdep, parentDEG = zpdeg.siggene.watdep, sigDiffRatio = interesth))
    barplotExp(intersect(interesth, intersect(zideg.siggene.watdep, zpdeg.siggene.watdep)), 'target_waterdeprivation')
    WriteExcel(list(Sheet1 = BindTairName(intersect(interesth, intersect(zideg.siggene.watdep, zpdeg.siggene.watdep)))),
                     './result/target_waterdeprivation.xlsx')
    
}



stop()


D <- NULL
D <- load.expression(D)
D <- calc.fpkm(D)
D <- calc.exp.genes(D)
ak.deg(D)



# expressed genes
png(paste0(FIG_PATH, '/exp.H.homeologs.png'), 800, 500)
barplot(colSums(D$exp.homeolog$H), las = 2)
dev.off()
png(paste0(FIG_PATH, '/exp.A.homeologs.png'), 800, 500)
barplot(colSums(D$exp.homeolog$A), las = 2)
dev.off()


# PCA
fpkm.mat <- cbind(D$fpkm$H, D$fpkm$A)
fpkm.col <- paste0(c(rep('H', ncol(D$fpkm$H)), rep('A', ncol(D$fpkm$A))), '_', colnames(fpkm.mat))
fpkm.mat <- fpkm.mat[apply(D$exp.homeolog$H | D$exp.homeolog$A, 1, any), ]
is.zero.lib <- colSums(fpkm.mat) == 0
fpkm.mat <- fpkm.mat[, !is.zero.lib]
fpkm.col <- fpkm.col[!is.zero.lib]
colnames(fpkm.mat) <- fpkm.col

pcaobj <- prcomp(t(log10(as.matrix(fpkm.mat + 1))), scale = FALSE)

df <- data.frame(
    homeolog = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[1]]}),
    site = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[2]]}),
    strain = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[3]]}),
    date = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[4]]}),
    id = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[5]]}),
    pcaobj$x
)
df$strain <- as.character(df$strain)
df$date <- as.character(df$date)
df$strain[df$strain == 'A'] <- 'C. amara'
df$strain[df$strain == 'F'] <- 'C. flexuosa'
df$strain[df$strain == 'H'] <- 'C. hirsuta'
df$date[df$date == '0418'] <- 'April 18'
df$date[df$date == '0502'] <- 'May 2'
df$date[df$date == '0516'] <- 'May 16'
df$strain <- factor(df$strain, levels = c('C. flexuosa', 'C. amara', 'C. hirsuta'))
df$date   <- factor(df$date,   levels = c('April 18', 'May 2', 'May 16'))


png(paste0(FIG_PATH, '/pca.all.libs.png'), 800, 300)
v <- pcaobj$sdev ^ 2
names(v) <- colnames(pcaobj$x)
barplot(v / sum(v) * 100, las = 2)
dev.off()

g <- ggplot(df, aes(x = PC1, y = PC2, color = strain, shape = homeolog))
g <- g + geom_point() + scale_color_jama()
png(paste0(FIG_PATH, '/pca.flexuosahomeolog.libs.PC12.png'), 500, 500)
print(g)
dev.off()
g <- ggplot(df, aes(x = PC2, y = PC3, color = strain, shape = homeolog))
g <- g + geom_point() + scale_color_jama()
png(paste0(FIG_PATH, '/pca.flexuosahomeolog.libs.PC23.png'), 500, 500)
print(g)
dev.off()

g <- ggplot(df, aes(x = PC1, y = PC2, color = site, shape = homeolog))
g <- g + geom_point() + scale_color_uchicago()
png(paste0(FIG_PATH, '/pca.flexuosahomeolog.libs.PC12.site.png'), 500, 500)
print(g)
dev.off()
g <- ggplot(df, aes(x = PC2, y = PC3, color = site, shape = homeolog))
g <- g + geom_point() + scale_color_uchicago()
png(paste0(FIG_PATH, '/pca.flexuosahomeolog.libs.PC23.site.png'), 500, 500)
print(g)
dev.off()

g <- ggplot(df, aes(x = PC1, y = PC2, color = date, shape = homeolog))
g <- g + geom_point() + scale_color_nejm()
png(paste0(FIG_PATH, '/pca.flexuosahomeolog.libs.PC12.date.png'), 500, 500)
print(g)
dev.off()
g <- ggplot(df, aes(x = PC2, y = PC3, color = date, shape = homeolog))
g <- g + geom_point() + scale_color_nejm()
png(paste0(FIG_PATH, '/pca.flexuosahomeolog.libs.PC23.date.png'), 500, 500)
print(g)
dev.off()




# PCA for total PFKM
fpkm.mat <- D$total.fpkm
fpkm.mat <- fpkm.mat[apply(fpkm.mat > 1, 1, any), ]
pcaobj <- prcomp(t(log10(as.matrix(fpkm.mat + 1))), scale = FALSE)

df <- data.frame(
    site = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[1]]}),
    strain = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[2]]}),
    date = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[3]]}),
    id = sapply(strsplit(rownames(pcaobj$x), '_'), function(x) {x[[4]]}),
    pcaobj$x
)
df$strain <- as.character(df$strain)
df$date <- as.character(df$date)
df$strain[df$strain == 'A'] <- 'C. amara'
df$strain[df$strain == 'F'] <- 'C. flexuosa'
df$strain[df$strain == 'H'] <- 'C. hirsuta'
df$date[df$date == '0418'] <- 'April 18'
df$date[df$date == '0502'] <- 'May 2'
df$date[df$date == '0516'] <- 'May 16'
df$strain <- factor(df$strain, levels = c('C. flexuosa', 'C. amara', 'C. hirsuta'))
df$date   <- factor(df$date,   levels = c('April 18', 'May 2', 'May 16'))


png(paste0(FIG_PATH, '/pca.flexuosaall.libs.png'), 800, 300)
v <- pcaobj$sdev ^ 2
names(v) <- colnames(pcaobj$x)
barplot(v / sum(v) * 100, las = 2)
dev.off()

g <- ggplot(df, aes(x = PC1, y = PC2, color = strain))
g <- g + geom_point() + scale_color_jama()
png(paste0(FIG_PATH, '/pca.flexuosaall.libs.PC12.png'), 500, 500)
print(g)
dev.off()
g <- ggplot(df, aes(x = PC2, y = PC3, color = strain))
g <- g + geom_point() + scale_color_jama()
png(paste0(FIG_PATH, '/pca.flexuosaall.libs.PC23.png'), 500, 500)
print(g)
dev.off()

g <- ggplot(df, aes(x = PC1, y = PC2, color = site))
g <- g + geom_point() + scale_color_uchicago()
png(paste0(FIG_PATH, '/pca.flexuosaall.libs.PC12.site.png'), 500, 500)
print(g)
dev.off()
g <- ggplot(df, aes(x = PC2, y = PC3, color = site))
g <- g + geom_point() + scale_color_uchicago()
png(paste0(FIG_PATH, '/pca.flexuosaall.libs.PC23.site.png'), 500, 500)
print(g)
dev.off()

g <- ggplot(df, aes(x = PC1, y = PC2, color = date))
g <- g + geom_point() + scale_color_nejm()
png(paste0(FIG_PATH, '/pca.flexuosaall.libs.PC12.date.png'), 500, 500)
print(g)
dev.off()
g <- ggplot(df, aes(x = PC2, y = PC3, color = date))
g <- g + geom_point() + scale_color_nejm()
png(paste0(FIG_PATH, '/pca.flexuosaall.libs.PC23.date.png'), 500, 500)
print(g)
dev.off()





# FPKM correlation
total.fpkm <- D$total.fpkm
total.fpkm <- total.fpkm[apply(total.fpkm > 1, 1, any), ]
p.total.fpkm <- log10(total.fpkm[, - grep('_F_', colnames(total.fpkm))] + 1)
f.total.fpkm <- log10(total.fpkm[, grep('_F_', colnames(total.fpkm))] + 1)
cc <- matrix(0, nrow = ncol(p.total.fpkm), ncol = ncol(f.total.fpkm))
colnames(cc) <- colnames(f.total.fpkm)
rownames(cc) <- colnames(p.total.fpkm)
for (i in 1:nrow(cc)) {
    for (j in 1:ncol(cc)) {
        cc[i, j] <- cor(p.total.fpkm[, i], f.total.fpkm[, j], method = 'pearson')
    }
}
png(paste0(FIG_PATH, '/heatmap_totallogfpkm_alldate.png'), 800, 300)
heatmap.2(as.matrix(1 - cc), trace = "none", scale = "none", col = COLSFUNCC2(),
          hclustfun = function(x) hclust(x, method = 'ward'), key.title = '',
          density.info = "none", key.xlab = "1 - rho", margins = c(8, 10),
          cexCol = 0.8, cexRow = 0.8)
dev.off()


total.fpkm <- D$total.fpkm
total.fpkm <- total.fpkm[apply(total.fpkm > 1, 1, any), ]
total.fpkm <- total.fpkm[, grep('0502', colnames(total.fpkm))]
p.total.fpkm <- log10(total.fpkm[, - grep('_F_', colnames(total.fpkm))] + 1)
f.total.fpkm <- log10(total.fpkm[, grep('_F_', colnames(total.fpkm))] + 1)
cc <- matrix(0, nrow = ncol(p.total.fpkm), ncol = ncol(f.total.fpkm))
colnames(cc) <- colnames(f.total.fpkm)
rownames(cc) <- colnames(p.total.fpkm)
for (i in 1:nrow(cc)) {
    for (j in 1:ncol(cc)) {
        cc[i, j] <- cor(p.total.fpkm[, i], f.total.fpkm[, j], method = 'pearson')
    }
}
png(paste0(FIG_PATH, '/heatmap_totallogfpkm_0502.png'), 480, 300)
heatmap.2(as.matrix(1 - cc), trace = "none", scale = "none", col = COLSFUNCC2(),
          hclustfun = function(x) hclust(x, method = 'ward'), key.title = '',
          density.info = "none", key.xlab = "1 - rho", margins = c(8, 10),
          cexCol = 0.8, cexRow = 0.8)
dev.off()


# FPKM correlation
total.fpkm <- cbind(D$fpkm$H, D$fpkm$A)
colnames(total.fpkm) <- paste0(colnames(total.fpkm), rep(c('_H', '_A'), each = ncol(total.fpkm)/2))
total.fpkm <- total.fpkm[apply(total.fpkm > 1, 1, any), ]
total.fpkm <- total.fpkm[, colSums(total.fpkm) > 0]
p.total.fpkm <- log10(total.fpkm[, - grep('_F_', colnames(total.fpkm))] + 1)
f.total.fpkm <- log10(total.fpkm[, grep('_F_', colnames(total.fpkm))] + 1)
cc <- matrix(0, nrow = ncol(p.total.fpkm), ncol = ncol(f.total.fpkm))
colnames(cc) <- colnames(f.total.fpkm)
rownames(cc) <- colnames(p.total.fpkm)
for (i in 1:nrow(cc)) {
    for (j in 1:ncol(cc)) {
        cc[i, j] <- cor(p.total.fpkm[, i], f.total.fpkm[, j], method = 'pearson')
    }
}
png(paste0(FIG_PATH, '/heatmap_homeologfpkm_alldate.png'), 1500, 300)
heatmap.2(as.matrix(1 - cc), trace = "none", scale = "none", col = COLSFUNCC2(),
          hclustfun = function(x) hclust(x, method = 'ward'), key.title = '',
          density.info = "none", key.xlab = "1 - rho", margins = c(8, 14),
          cexCol = 0.8, cexRow = 0.8)
dev.off()



# FPKM correlation
total.fpkm <- cbind(D$fpkm$H, D$fpkm$A)
total.fpkm <- total.fpkm[, grep('0502', colnames(total.fpkm))]
colnames(total.fpkm) <- paste0(colnames(total.fpkm), rep(c('_H', '_A'), each = ncol(total.fpkm)/2))
total.fpkm <- total.fpkm[apply(total.fpkm > 1, 1, any), ]
total.fpkm <- total.fpkm[, colSums(total.fpkm) > 0]
p.total.fpkm <- log10(total.fpkm[, - grep('_F_', colnames(total.fpkm))] + 1)
f.total.fpkm <- log10(total.fpkm[, grep('_F_', colnames(total.fpkm))] + 1)
cc <- matrix(0, nrow = ncol(p.total.fpkm), ncol = ncol(f.total.fpkm))
colnames(cc) <- colnames(f.total.fpkm)
rownames(cc) <- colnames(p.total.fpkm)
for (i in 1:nrow(cc)) {
    for (j in 1:ncol(cc)) {
        cc[i, j] <- cor(p.total.fpkm[, i], f.total.fpkm[, j], method = 'pearson')
    }
}
png(paste0(FIG_PATH, '/heatmap_homeologfpkm_0502.png'), 900, 300)
heatmap.2(as.matrix(1 - cc), trace = "none", scale = "none", col = COLSFUNCC2(),
          hclustfun = function(x) hclust(x, method = 'ward'), key.title = '',
          density.info = "none", key.xlab = "1 - rho", margins = c(8, 14),
          cexCol = 0.8, cexRow = 0.8)
dev.off()













