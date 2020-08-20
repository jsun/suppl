## R code for homeolog expression analysis

set.seed(2016)
library(MASS)
library(XLConnect)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggExtra)
library(gplots)
library(VennDiagram)
library(GO.db)
library(org.At.tair.db)
library(clusterProfiler)


options(stringsAsFactors = FALSE)


## GLOBAL PARAMTERS
FPKM_CUTOFF <- 1.0


## INIT DIRECTORIES
WORK_DIR <- dirname('./')
PATH <- list(
    res = paste0(WORK_DIR, "/result_files"),
    lib = paste0(WORK_DIR, "/lib"),
    src = paste0(WORK_DIR, "/lib/src"),
    dat = paste0(WORK_DIR, "/data"),
    tmp = paste0(WORK_DIR, "/tmp")
)
for (.p in PATH) if(!file.exists(.p)) dir.create(.p)
path.lib      <- paste0(PATH$res, "/library_features")
path.de       <- paste0(PATH$res, "/differential_expression")
path.geneprof <- paste0(PATH$res, "/expression_profile")
path.expbias  <- paste0(PATH$res, "/expression_bias")
path.expbias.path      <- paste0(path.expbias, "/path")
for (.path.i in c(path.lib, path.de, path.geneprof, path.expbias, path.expbias.path)) {
    if (!file.exists(.path.i)) dir.create(.path.i)
}







## INIT COLOR SETTINGS
.brewer.set1    <- brewer.pal(9, "Set1")
.brewer.pastel1 <- brewer.pal(9, "Pastel1")
.brewer.dark2   <- brewer.pal(8, "Dark2")
.brewer.greens  <- brewer.pal(9, "Greens")
COLS <- list(
    ama  = .brewer.set1[2],     # C. amara
    riv  = .brewer.set1[1],     # C. rivularis
    hir  = .brewer.set1[1],     # C. hirsuta
    ins  = .brewer.set1[3],     # C. insueta
    insA = .brewer.set1[4],     # C. insueta (A)
    insR = .brewer.set1[5],     # C. insueta (R)
    DE   = .brewer.set1[5],     # Differentially expressed genes
    NDE  = .brewer.set1[9]      # Non-differentially expressed genes
)

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
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu"))[5:9])
    .colfunc(n)
}



SPLABELS  <- c("AA", "IA", "IR", "RR")
SPLABELS2 <- c("AA", "IA", "II", "IR", "RR")
TIMELABELS  <- c("0 hr", "2 hr", "4 hr", "8 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr")
TIMELABELS2 <-  c("0 hr vs. 2 hr",   "0 hr vs. 4 hr",  "0 hr vs. 8 hr", "0 hr vs. 12 hr",
                  "0 hr vs. 24 hr", "0 hr vs. 48 hr", "0 hr vs. 72 hr", "0 hr vs. 96 hr")
INCH2CM <- 1 / 2.54





#' Split string of TAIR ID and convert them to CARHR ID.
#'
#'
TAIRSTR2CARHRVEC <- function(x, sep = '/') {
    paste0(unlist(TAIR2CARHR[unlist(strsplit(x, '/'))]), collapse = sep)
}



#' Plot multiple ggplot objects in a figure
#'
#'
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#' Save list object into Excel file.
#'
#'
SaveExcel <- function(data.list, file.name = NULL) {
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
    workbook <- loadWorkbook(file.name, create = TRUE)
    createSheet(workbook, names(dat))
    writeWorksheet(workbook, dat, names(dat), header = TRUE)
    saveWorkbook(workbook)
    gc()
}





GetDesign <- function(tissue = "leaf", species = NULL) {
    design <- read.table(paste0(PATH$lib, "/src/lib_design.txt"), header = TRUE)
    keeped.libs <- (design$tissue == tissue)
    design      <- design[keeped.libs, ]
    libs.order  <- order(design$tissue, design$species, design$time)
    design      <- design[libs.order, ]
    design.sp   <- NULL
    if (!is.null(species)) {
        for (sp in species) {
            design.sp <- rbind(design.sp, design[grep(sp, design$species), ])
        }
        design <- design.sp
        design$species <- factor(design$species, levels = species)
    }
    design
}




GetCounts <- function(ca, design) {
    counts.aorigin <- counts.rorigin <- counts.common <- 
        vector("list", length = length(design$library))
    names(counts.aorigin) <- names(counts.rorigin) <- names(counts.common) <-
        design$library
    
    gene.names <- NULL
    
    for (libname in design$library) {
        .a <- read.table(paste0(PATH$dat, "/counts/", libname,
                                "__on__camara.filtered.rc_orig.sorted.txt"),
                         sep = "\t", fill = TRUE)
        .r <- read.table(paste0(PATH$dat, "/counts/", libname,
                                "__on__camara.filtered.rc_other.sorted.txt"),
                         sep = "\t", fill = TRUE)
        .c <- read.table(paste0(PATH$dat, "/counts/", libname,
                                "__on__camara.filtered.rc_common.sorted.txt"),
                         sep = "\t", fill = TRUE)
        gene.names <- .a[, 1]
        counts.aorigin[[libname]] <- .a[, 2]
        counts.rorigin[[libname]] <- .r[, 2]
        counts.common[[libname]] <- .c[, 2]
    }
    
    counts.aorigin <- as.data.frame(counts.aorigin)[-c(grep("^__", gene.names)), ]
    counts.rorigin <- as.data.frame(counts.rorigin)[-c(grep("^__", gene.names)), ]
    counts.common  <- as.data.frame(counts.common)[-c(grep("^__", gene.names)), ]
    
    print(data.frame(A      = colSums(counts.aorigin),
                     common = colSums(counts.common),
                     R      = colSums(counts.rorigin)))
    
    gene.names  <- gene.names[-c(grep("^__", gene.names))]
    
    counts.AA <- counts.aorigin[, 1:9] + counts.rorigin[, 1:9] + counts.common[, 1:9]
    counts.RR <- counts.aorigin[, 19:27] + counts.rorigin[, 19:27] + counts.common[, 19:27]
    
    aorigin.ratio <- (counts.aorigin / (counts.aorigin + counts.rorigin))[, 10:18]
    aorigin.ratio[(counts.aorigin[, 10:18] + counts.rorigin[, 10:18]) == 0 & counts.common[, 10:18] != 0] <- 1/3
    counts.IA <- counts.aorigin[, 10:18] + counts.common[, 10:18] * aorigin.ratio
    counts.IR <- counts.rorigin[, 10:18] + counts.common[, 10:18] * (1 - aorigin.ratio)
    counts.IA[is.na(counts.IA)] <- 0
    counts.IR[is.na(counts.IR)] <- 0
    
    ## add gene names and time course labels
    colnames(counts.AA) <- colnames(counts.RR) <- colnames(counts.IA) <- colnames(counts.IR) <-
            colnames(aorigin.ratio) <- TIMELABELS
    rownames(counts.AA) <- rownames(counts.RR) <- rownames(counts.IA) <- rownames(counts.IR) <-
            rownames(aorigin.ratio) <- gene.names
    
    ## remove all zero rows
    is.allzero <- as.logical(rowSums(counts.aorigin + counts.rorigin + counts.common) == 0)
    gene.names <- gene.names[!is.allzero]
    counts.AA <- counts.AA[!is.allzero, ]
    counts.RR <- counts.RR[!is.allzero, ]
    counts.IA <- counts.IA[!is.allzero, ]
    counts.IR <- counts.IR[!is.allzero, ]
    aorigin.ratio <- aorigin.ratio[!is.allzero, ]
    
    ca$counts <- list(AA = as.matrix(counts.AA), IA = as.matrix(counts.IA),
                      IR = as.matrix(counts.IR), RR = as.matrix(counts.RR))
    ca$aorigin.ratio <- as.matrix(aorigin.ratio)
    ca$gene.names <- gene.names
    ca
}

CalcFPKM <- function(ca) {
    gene.length <- read.table(paste0(PATH$lib, "/src/gene_length.tsv"), header = FALSE)
    row.names(gene.length) <- gene.length[, 1]
    gene.length <- gene.length[ca$gene.names, -1]
    
    fpkm <- ca$counts
    
    ## AA
    cpm.A <- sweep(ca$counts$AA, 2, 1e6 / colSums(ca$counts$AA), "*")
    fpkm$AA <- sweep(cpm.A, 1, 1e3 / gene.length[, 1], "*")
    
    ## RR
    cpm.R <- sweep(ca$counts$RR, 2, 1e6 / colSums(ca$counts$RR), "*")
    fpkm$RR <- sweep(cpm.R, 1, 1e3 / gene.length[, 1], "*")
    
    ## IA, IR
    print("IA counts ---> IA FPKM")
    print("IR counts ---> IR FPKM")
    print("IA counts + IR counts---> II FPKM")
    cpm.IA <- sweep(ca$counts$IA, 2, 1e6 / colSums(ca$counts$IA), "*")
    fpkm$IA <- sweep(cpm.IA, 1, 1e3 / gene.length[, 1], "*")
    cpm.IR <- sweep(ca$counts$IR, 2, 1e6 / colSums(ca$counts$IR), "*")
    fpkm$IR <- sweep(cpm.IR, 1, 1e3 / gene.length[, 1], "*")
    .counts.I <- ca$counts$IA + ca$counts$IR
    cpm.II <- sweep(.counts.I, 2, 1e6 / colSums(.counts.I), "*")
    fpkm$II <- sweep(cpm.II, 1, 1e3 / gene.length[, 1], "*")
   
    ca$gene.length <- gene.length
    ca$fpkm <- fpkm
    ca
}





DoGO2 <- function(sig.genes = NULL, all.genes = NULL, ontology = "BP",
                  p.cutoff = 1, q.cutoff = 1) {
    all.genes <- as.character(all.genes)
    sig.genes <- as.character(sig.genes)
    all.genes.at <- as.character(unlist(CARHR2TAIR[all.genes]))
    sig.genes.at <- as.character(unlist(CARHR2TAIR[sig.genes]))
    all.genes.at <- all.genes.at[all.genes.at != ""]
    sig.genes.at <- sig.genes.at[sig.genes.at != ""]
    ego <- enrichGO(gene = sig.genes.at, universe = all.genes.at, OrgDb = org.At.tair.db,
                    ont = ontology, pAdjustMethod = "BH", keyType = "TAIR",
                    pvalueCutoff = p.cutoff, qvalueCutoff = q.cutoff, readable = FALSE)
    egos <- dropGO(ego, level = 1)
    egos <- dropGO(egos, level = 2)
    egos <- dropGO(egos, level = 3)
    egostable <- as.data.frame(egos)
    numgenes <- sapply(strsplit(egostable$BgRatio, "/"), function(x) as.integer(x[1]))
    egostable <- egostable[(10 < numgenes) & (numgenes < 500), ]
    egostable <- egostable[, colnames(egostable) != "p.adjust"]
    egostable$qvalue <- p.adjust(egostable$pvalue, method = "BH")
    
    .carhr <- rep("", length = nrow(egostable))
    for (.i in 1:nrow(egostable)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(egostable$geneID[.i])
    }
    egostable <- data.frame(egostable, carID = .carhr)
    
    list(GO = ego, GOTABLE = as.data.frame(ego), GOSMPLTABLE = egostable)
}









LargeFPKM_GO <- function(ca) {
    A <- ca$fpkm$AA
    R <- ca$fpkm$RR
    I <- ca$fpkm$IA + ca$fpkm$IR
    A <- A[rowSums(A > FPKM_CUTOFF) > 0, ]
    R <- R[rowSums(R > FPKM_CUTOFF) > 0, ]
    I <- I[rowSums(I > FPKM_CUTOFF) > 0, ]
    A <- log10(A+1)
    R <- log10(R+1)
    I <- log10(I+1)
    sd.A <- apply(A, 1, sd, na.rm = T)
    sd.R <- apply(R, 1, sd, na.rm = T)
    sd.I <- apply(I, 1, sd, na.rm = T)
    mu.A <- apply(A, 1, mean, na.rm = T)
    mu.R <- apply(R, 1, mean, na.rm = T)
    mu.I <- apply(I, 1, mean, na.rm = T)
    cv.A <- sd.A / mu.A
    cv.R <- sd.R / mu.R
    cv.I <- sd.I / mu.I
    
    # camara
    kamara <- cv.A > 0.20 & mu.A > 1
    goobj <- DoGO(sig.genes = names(kamara[kamara]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableama<- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableama))
    for (.i in 1:nrow(gotableama)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableama$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableama, .carhr)),
              file.name = paste0(path.expbias, "/GO_camara_CV0.20_and_mu1.0.xlsx"))
    gc()
    # crivularis
    krivularis <- cv.R > 0.20 & mu.R > 1
    goobj <- DoGO(sig.genes = names(krivularis[krivularis]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableriv <- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableriv))
    for (.i in 1:nrow(gotableriv)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableriv$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableriv, .carhr)),
              file.name = paste0(path.expbias, "/GO_rivularis_CV0.20_and_mu1.0.xlsx"))
    gc()
    # cinsueta
    kinsueta <- cv.I > 0.20 & mu.I > 1
    goobj <- DoGO(sig.genes = names(kinsueta[kinsueta]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableins <- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableins))
    for (.i in 1:nrow(gotableins)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableins$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableins,.carhr)),
              file.name = paste0(path.expbias, "/GO_cinuseta_CV0.20_and_mu1.0.xlsx"))
    gc()
    
    # cinsueta (var AOriginRatio)
    theta.ar <- ca$aorigin.ratio[names(mu.I), ]
    var.theta.ar <- apply(theta.ar, 1, var, na.rm = T)
    keep2c <- var.theta.ar > 0.01  & mu.I > 1
    goobj <- DoGO(sig.genes = names(keep2c[keep2c]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableins2 <- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableins2))
    for (.i in 1:nrow(gotableins2)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableins2$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableins2, .carhr)),
              file.name = paste0(path.expbias, "/GO_cinsueta_VarARatio0.01_and_mu1.0.v2.xlsx"))

    
    goterms <- unique(c(gotableama$ID[gotableama$qvalue < 0.1],
                        gotableriv$ID[gotableriv$qvalue < 0.1],
                        gotableins$ID[gotableins$qvalue < 0.1],
                        gotableins2$ID[gotableins2$qvalue < 0.1]))
    
    go2desc <- rep("", length=length(goterms))
    gomat <- matrix(NA, ncol = 4, nrow = length(goterms))
    rownames(gomat) <- names(go2desc) <- goterms
    colnames(gomat) <- c("AA_CvMean", "II_CvMean", "II_VRatioMean", "RR_CvMean")
    go2desc[gotableama$ID[gotableama$qvalue < 0.1]] <- gotableama$Description[gotableama$qvalue < 0.1]
    go2desc[gotableriv$ID[gotableriv$qvalue < 0.1]] <- gotableriv$Description[gotableriv$qvalue < 0.1]
    go2desc[gotableins$ID[gotableins$qvalue < 0.1]] <- gotableins$Description[gotableins$qvalue < 0.1]
    go2desc[gotableins2$ID[gotableins2$qvalue < 0.1]] <- gotableins2$Description[gotableins2$qvalue < 0.1]
    
    gomat[gotableama$ID[gotableama$qvalue < 0.1], 1] <- -log10(gotableama$qvalue)[gotableama$qvalue < 0.1]
    gomat[gotableins$ID[gotableins$qvalue < 0.1], 2] <- -log10(gotableins$qvalue)[gotableins$qvalue < 0.1]
    gomat[gotableins2$ID[gotableins2$qvalue < 0.1], 3] <- -log10(gotableins2$qvalue)[gotableins2$qvalue < 0.1]
    gomat[gotableriv$ID[gotableriv$qvalue < 0.1], 4] <- -log10(gotableriv$qvalue)[gotableriv$qvalue < 0.1]
    
    
    gomat2 <- gomat
    gomat2[is.na(gomat2)] <- 1
    h <- heatmap.2(gomat2, trace = "none", scale = "none", dendrogram = "row", Colv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNCC2(), key.xlab = "log2(FPKM + 1)")
    df <- data.frame(gomat, desc = go2desc)
    df <- df[h$rowInd, ]
    SaveExcel(list(Sheet1=df), file.name = paste0(path.lib, "/CVLogFPKM0.2_meanLogFPKM1.0.xlsx"))
}






PlotRqsMDS <- function(ca) {
    rsqmat <- matrix(NA, ncol = 4 * 9, nrow = 4 * 9)
    f <- cbind(ca$fpkm$AA, ca$fpkm$IA, ca$fpkm$IR, ca$fpkm$RR)
    f <- log10(f[ca$EXP$gene$ALL, ] + 1)
    
    sd.log10.f <- apply(f, 1, sd, na.rm = T)
    mu.log10.f <- apply(f, 1, mean, na.rm = T)
    cv.log10.f <- sd.log10.f / mu.log10.f
    
    #f <- f[(cv.log10.f > 1.5), ]
    
    for (i in 2:(4*9)) {
        for (j in 1:(i-1)) {
            rsqmat[i, j] <- summary(lm(f[, i] ~ f[, j]))$r.squared
        }
    }
    
    colnames(rsqmat) <- rownames(rsqmat) <- 
        paste0(rep(SPLABELS, each = 9), "__", rep(TIMELABELS, times = 4))
    
    mds <- data.frame(cmdscale(as.dist(1 - rsqmat)))
    mdsdf <- data.frame(
        x = mds$X1, y = mds$X2,
        species = sapply(strsplit(rownames(mds), "__"), function(x) x[1]),
        time = gsub(" hr", "", sapply(strsplit(rownames(mds), "__"), function(x) x[2]))
    )
    g <- ggplot(mdsdf, aes(x = x, y = y, color = species, label = time))
    g <- g + geom_text(size=6) + coord_fixed()
    g <- g + xlab("X1") + ylab("X2")
    g <- g + scale_color_manual(values = c(COLS$ama, COLS$insA, COLS$insR, COLS$riv))
    g <- g + theme_bw()
    g <- g + xlim(-0.5, 0.5) + ylim(-0.5, 0.5)
    g
    
    png(paste0(path.lib, "/logFPKM.MDS.png"), 340, 340)
    plot(g)
    dev.off()
}


PlotPCA <- function(ca) {

    f <- cbind(ca$fpkm$AA, ca$fpkm$IA, ca$fpkm$IR, ca$fpkm$RR)
    f <- log10(f[ca$EXP$gene$ALL, ] + 1)
    
    fpkm.class <- rep(c('C. amara', 'IA', 'IR', 'C. rivularis'), each = 9)
    time <- gsub(' hr', '', colnames(f))
    colnames(f) <- paste0(fpkm.class, ' ', time)

    pcaobj <- prcomp(t(f), scale = FALSE)
    # plot proportion of variances of PCA
    pca.prop.vars <- data.frame(Var = summary(pcaobj)$importance[2, ],
                                PC = factor(colnames(summary(pcaobj)$importance),
                                            levels = colnames(summary(pcaobj)$importance)))
    pca.prop.vars.ggplot <- ggplot(pca.prop.vars, aes(x = PC, y = Var))
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + geom_bar(stat = 'identity')
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + theme_bw()
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + xlab('principal component')
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + ylab('proportion of variance')
    png(paste0(path.lib, '/pca.prop.vars.png'), 800, 400)
    plot(pca.prop.vars.ggplot)
    dev.off()
    
    pca.pc.coords <- data.frame(pcaobj$x,
                                fpkm = fpkm.class,
                                time = time)
    pca.pc.coords$time <- factor(pca.pc.coords$time,
                                 levels = c('0', '2', '4', '8', '12', '24', '48', '72', '96'))

    pca.pc.coords.ggplot <- ggplot(pca.pc.coords, aes(x = PC1, y = PC2, shape = time, color = fpkm))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + geom_point() + scale_color_nejm() + coord_fixed() 
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + scale_shape_manual(values=c(16, 17, 15, 1:6))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + theme_bw() + guides(color = guide_legend(title = ''))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlab(paste0('PC1 (', round(summary(pcaobj)$importance[2, 1] * 100, 1), '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + ylab(paste0('PC2 (', round(summary(pcaobj)$importance[2, 2] * 100, 1), '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlim(-30, 30) + ylim(-30, 30)
    png(paste0(path.lib, '/pca.pc1pc2.png'), 920, 1400, res=250)
    plot(pca.pc.coords.ggplot)
    dev.off()

    pca.pc.coords.ggplot <- ggplot(pca.pc.coords, aes(x = PC2, y = PC3, color = fpkm, shape = time))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + geom_point() + scale_color_nejm() + coord_fixed()
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + scale_shape_manual(values=c(16, 17, 15, 1:6))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + theme_bw() + guides(color = guide_legend(title = ''))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlab(paste0('PC2 (', round(summary(pcaobj)$importance[2, 1] * 100, 1), '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + ylab(paste0('PC3 (', round(summary(pcaobj)$importance[2, 3] * 100, 1), '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlim(-30, 30) + ylim(-30, 30)
    png(paste0(path.lib, '/pca.pc2pc3.png'), 920, 1400, res=250)
    plot(pca.pc.coords.ggplot)
    dev.off()

}





PlotAoriginRatioDist <- function(ca) {
    r <- ca$aorigin.ratio
    dfr <- NULL
    for (i in 1:9) {
        target.homeologs <- (ca$fpkm$IA[, i] > FPKM_CUTOFF) | (ca$fpkm$IR[, i] > FPKM_CUTOFF)
        ri <- r[target.homeologs, i]
        dfi <- data.frame(value = ri, gene = names(ri), time = paste0(TIMELABELS[i], " (", length(ri), " homeologs)"))
        dfr <- rbind(dfr, dfi)
    }
    
    dfr$time <- factor(dfr$time, levels = unique(dfr$time))
    
    g <- ggplot(dfr, aes(x = value))
    g <- g + geom_histogram(binwidth = 0.025)
    g <- g + theme_bw()
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none",
                   strip.background = element_rect(fill = "white", colour = "white"))
    g <- g + scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), minor_breaks = NULL)
    g <- g + geom_vline(xintercept = 1/3, color = '#D55E00')
    g <- g + facet_wrap(~ time, ncol = 3)
    g <- g + xlab("A-origin ratio") + ylab("frequency")
    ca$BIAS$aorigin.ratio.dist <- g
    
    pdf(paste0(path.expbias, "/AOriginRatio.hist.lib.pdf"), 6, 6)
    options(digits=2)
    plot(g)
    dev.off()
    options(digits=6)
    
     
    R <- (ca$fpkm$AA/2) / (ca$fpkm$AA/2 + ca$fpkm$RR)
    dfr <- NULL
    for (i in 1:9) {
        target.homeologs <- (ca$fpkm$AA[, i] > FPKM_CUTOFF) | (ca$fpkm$RR[, i] > FPKM_CUTOFF)
        #target.homeologs <- (ca$fpkm$AA[, i] / 2 + ca$fpkm$RR[, i] > FPKM_CUTOFF)
        ri <- R[target.homeologs, i]
        dfi <- data.frame(value = ri, gene = names(ri), time = paste0(TIMELABELS[i], " (", length(ri), " homeologs)"))
        dfr <- rbind(dfr, dfi)
    }

    dfr$time <- factor(dfr$time, levels = paste0(TIMELABELS, " (", apply(target.homeologs, 2, sum), " homeologs)"))

    g <- ggplot(dfr, aes(x = value))
    g <- g + geom_histogram(binwidth = 0.025)
    g <- g + theme_bw()
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none",
                   strip.background = element_rect(fill = "white", colour = "white"))
    g <- g + scale_x_continuous(breaks = c(0, 0.2, 1/3, 0.4, 0.6, 0.8, 1), minor_breaks = NULL)
    g <- g + facet_grid(~ time)#, ncol = 5)
    g <- g + xlab("A-origin ratio") + ylab("frequency")

    png(paste0(path.expbias, "/AOriginRatio.parent.hist.lib.png"), 1000, 240)
    options(digits=2)
    plot(g)
    dev.off()
    options(digits=6)

    
    
    dfrR <- NULL
    crrR <- NULL
    for (i in 1:9) {
        .r <- r[, i]
        .R <- R[, i]
        .k <- rowSums(cbind(ca$fpkm$IA[, i], ca$fpkm$AA[, i], ca$fpkm$RR[, i], ca$fpkm$IR[, i]) > FPKM_CUTOFF) > 0
        .k.isval <- !is.na(rowSums(cbind(.r, .R)))
        dfrR <- rbind(dfrR, data.frame(r = .r[.k & .k.isval], R = .R[.k & .k.isval], timepoint = TIMELABELS[i]))
        crrR <- c(crrR, cor(.r[.k & .k.isval], .R[.k & .k.isval]))
    }
    dfrR$timepoint <- factor(dfrR$timepoint, levels = TIMELABELS)
    g <- ggplot(dfrR, aes(x = r, y = R))
    g <- g + geom_point(alpha = 0.4, size = 0.6) + theme_bw()
    g <- g + xlab('A-origin ratio of C. insueta') + ylab('A-origin ratio of parental species')
    g <- g + facet_wrap(~ timepoint, ncol = 3) + coord_fixed()
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                   strip.background = element_rect(fill = "white", colour = "white"))
    print(crrR)
    png(paste0(path.expbias, "/AOriginRatio.parent.and.child.png"), 800, 800)
    options(digits=2)
    plot(g)
    dev.off()
    options(digits=6)

    
    
    dfr <- NULL
    for (i in 1:9) {
        target.homeologs <- (ca$fpkm$IA[, i] > FPKM_CUTOFF) | (ca$fpkm$IR[, i] > FPKM_CUTOFF)
        ri <- r[target.homeologs, i]
        dfi <- data.frame(value = ri, gene = names(ri), time = TIMELABELS[i])
        dfr <- rbind(dfr, dfi)
    }
    
    gene2chr_txt <- read.table(paste0(PATH$lib, "/src/gene_chr.txt"), header = F, sep = "\t")
    gene2chr <- gene2chr_txt[, 2]
    names(gene2chr) <- gene2chr_txt[, 1]
    dfr <- data.frame(dfr, chr = gene2chr[dfr$gene])
    dfr <- dfr[grep("Chr", dfr$chr), ]
    dfr <- dfr[dfr$chr != 'ChrC', ]
    dfr <- dfr[dfr$chr != 'ChrM', ]
    g <- ggplot(dfr, aes(x = value))
    g <- g + geom_histogram(binwidth = 0.025)
    g <- g + theme_bw()
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none",
                   strip.background = element_rect(fill = "white", colour = "white"))
    g <- g + scale_x_continuous(breaks = c(0, 0.2, 1/3, 0.4, 0.6, 0.8, 1), minor_breaks = NULL)
    g <- g + facet_grid(chr ~ time, scale="free_y")
    g <- g + xlab("A-origin ratio") + ylab("frequency")
    png(paste0(path.expbias, "/AOriginRatio.hist.chr.png"), 1500, 1200)
    options(digits=2)
    plot(g)
    dev.off()
    options(digits=6)

    invisible(ca)
}





PlotProfile <- function(g, ca, cls = TRUE, m = c(10, 10), col = NULL) {
    g <- intersect(g, rownames(ca$fpkm$AA))
    g.fpkm.A <- ca$fpkm$AA[g, ]
    g.fpkm.a <- ca$fpkm$IA[g, ]
    g.fpkm.r <- ca$fpkm$IR[g, ]
    g.fpkm.R <- ca$fpkm$RR[g, ]
    g.fpkm <- cbind(cbind(g.fpkm.A, g.fpkm.a), cbind(g.fpkm.r, g.fpkm.R))
    df <- BindTairName(g)[, -4]
    df[df[, 3] == "NA", 3] <- df[df[, 3] == "NA", 2]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 2]
    df[, 3] <- gsub("\\|\\|", ",", df[, 3])
    df[df[, 2] == "", 2] <- df[df[, 2] == "", 1]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 1]
    if (!is.null(col)) df[, 3] <- paste0("[", col, "] ", df[, 3])
    colnames(g.fpkm) <- paste0(rep(SPLABELS, each = 9), " (", TIMELABELS, ")")
    rownames(g.fpkm) <- df[, 3]
    g.fpkm <- log10(g.fpkm + 1)
    if (cls) {
        heatmap.2(g.fpkm, trace = "none", scale = "row", dendrogram = "row", Colv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "z-score", margins = m,
              colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    } else {
        heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "none", Colv = F, Rowv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "z-score", margins = m,
              colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    }
}


PlotProfile.FC <- function(g, ca, cls = TRUE, m = c(10, 10), col = NULL) {
    g <- intersect(g, rownames(ca$fpkm$AA))
    g.fpkm.A <- log2(ca$fpkm$AA[g, ] + 1)
    g.fpkm.a <- log2(ca$fpkm$IA[g, ] + 1)
    g.fpkm.r <- log2(ca$fpkm$IR[g, ] + 1)
    g.fpkm.R <- log2(ca$fpkm$RR[g, ] + 1)
    g.fpkm.A <- (g.fpkm.A - g.fpkm.A[, 1])#[, -1]
    g.fpkm.a <- (g.fpkm.a - g.fpkm.a[, 1])#[, -1]
    g.fpkm.r <- (g.fpkm.r - g.fpkm.r[, 1])#[, -1]
    g.fpkm.R <- (g.fpkm.R - g.fpkm.R[, 1])#[, -1]

    g.fpkm <- cbind(cbind(g.fpkm.A, g.fpkm.a), cbind(g.fpkm.r, g.fpkm.R))
    df <- BindTairName(g)[, -4]
    df[df[, 3] == "NA", 3] <- df[df[, 3] == "NA", 2]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 2]
    df[, 3] <- gsub("\\|\\|", ",", df[, 3])
    df[df[, 2] == "", 2] <- df[df[, 2] == "", 1]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 1]
    if (!is.null(col)) df[, 3] <- paste0("[", col, "] ", df[, 3])

    splabels <- c('C. amara', 'C. insueta (A-sub)', 'C. insueta (R-sub)', 'C. rivularis')
    colnames(g.fpkm) <- paste0(rep(splabels, each = 9), " (", TIMELABELS, ")")
    
    rownames(g.fpkm) <- df[, 3]
    g.fpkm <- g.fpkm[, -c(1, 9+1, 18+1, 27+1)]
    if (cls) {
        h <- heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "row", Colv = F, key.title = "",
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", margins = m, keysize = 1.2,
              colsep = c(8, 16, 24, 32), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    } else {
        h <- heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "none", Colv = F, Rowv = F, key.title = "",
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", margins = m, keysize = 1.2,
              colsep = c(8, 16, 24, 32), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    }
    invisible(h)
}




CalcExpressedGenes <- function(ca) {
    exp.AA <- rownames(ca$fpkm$AA)[apply(ca$fpkm$AA > FPKM_CUTOFF, 1, any)]
    exp.RR <- rownames(ca$fpkm$RR)[apply(ca$fpkm$RR > FPKM_CUTOFF, 1, any)]
    exp.IA <- rownames(ca$fpkm$IA)[apply(ca$fpkm$IA > FPKM_CUTOFF, 1, any)]
    exp.IR <- rownames(ca$fpkm$IR)[apply(ca$fpkm$IR > FPKM_CUTOFF, 1, any)]
    exp.all <- unique(c(exp.AA, exp.RR, exp.IA, exp.IR))
    print(paste0("Overlap(AA, RR): ",  round(length(intersect(exp.AA, exp.RR))/length(union(exp.AA, exp.RR)), 3)))
    print(paste0("Overlap(IA, IR): ",  round(length(intersect(exp.IA, exp.IR))/length(union(exp.IA, exp.IR)), 3)))
    print(paste0("Overlap(AA, IA): ",  round(length(intersect(exp.AA, exp.IA))/length(union(exp.AA, exp.IA)), 3)))
    print(paste0("Overlap(RR, IR): ",  round(length(intersect(exp.RR, exp.IR))/length(union(exp.RR, exp.IR)), 3)))
    
    png(paste0(path.lib, "/Exp-genes-venn.png"), 500, 500)
    v <- venn(list(AA = exp.AA, IA = exp.IA, IR = exp.IR, RR = exp.RR))
    dev.off()
    
    genes <- gofull <- gosimple <- vector("list", length = length(attr(v, "intersections")))
    names(genes) <- names(gofull) <- names(gosimple) <- c(names(attr(v, "intersections")))
    for (vn in names(attr(v, "intersections"))) {
        .g <- attr(v, "intersections")[[vn]]
        .go <- DoGO2(sig.genes = .g, all.genes = exp.all)
        gofull[[vn]] <- .go$GOTABLE
        gosimple[[vn]] <- .go$GOSMPLTABLE
        genes[[vn]] <- BindTairName(.g)
    }
    
    names(gofull)    <- gsub(":", "-", names(gofull))
    names(gosimple) <- gsub(":", "-", names(gosimple))
    names(genes)     <- gsub(":", "-", names(genes))
    SaveExcel(gofull, paste0(path.lib, "/Exp-genes-GO-full.xlsx"))
    gc()
    SaveExcel(gosimple, paste0(path.lib, "/Exp-genes-GO-simple.xlsx"))
    gc()
    SaveExcel(genes, paste0(path.lib, "/Exp-genes.xlsx"))
    gc() 
    
    # expressed genes
    ca$EXP$genes <- list(AA = exp.AA, RR = exp.RR,
                         IA = exp.IA, IR = exp.IR,
                         ALL = exp.all)
    
    neo.exp.IA <- setdiff(exp.IA, exp.AA)
    neo.exp.IR <- setdiff(exp.IR, exp.RR)
    sup.exp.IA <- setdiff(exp.AA, exp.IA)
    sup.exp.IR <- setdiff(exp.RR, exp.IR)
    
    print(paste0("Newly expressed in A: ", length(neo.exp.IA)))
    print(paste0("Suppressed in A:      ", length(sup.exp.IA)))
    print(paste0("Newly expressed in R: ", length(neo.exp.IR)))
    print(paste0("Suppressed in R:      ", length(sup.exp.IR)))
    
    print(length(intersect(neo.exp.IA, neo.exp.IR))/length(union(neo.exp.IA, neo.exp.IR)))
    print(length(intersect(sup.exp.IA, sup.exp.IR))/length(union(sup.exp.IA, sup.exp.IR)))
    
    neo.exp.IA.go <- DoGO2(sig.genes = neo.exp.IA, all.genes = exp.all)
    neo.exp.IR.go <- DoGO2(sig.genes = neo.exp.IR, all.genes = exp.all)
    sup.exp.IA.go <- DoGO2(sig.genes = sup.exp.IA, all.genes = exp.all)
    sup.exp.IR.go <- DoGO2(sig.genes = sup.exp.IR, all.genes = exp.all)
    
    neo.IA.df <- list(gene = BindTairName(neo.exp.IA),
                      GOFULL = neo.exp.IA.go$GOTABLE,
                      GOSMPL = neo.exp.IA.go$GOSMPLTABLE)
    neo.IR.df <- list(gene = BindTairName(neo.exp.IR),
                      GOFULL = neo.exp.IR.go$GOTABLE,
                      GOSMPL = neo.exp.IR.go$GOSMPLTABLE)
    sup.IA.df <- list(gene = BindTairName(sup.exp.IA),
                      GOFULL= sup.exp.IA.go$GOTABLE,
                      GOSMPL = sup.exp.IA.go$GOSMPLTABLE)
    sup.IR.df <- list(gene = BindTairName(sup.exp.IR),
                      GOFULL= sup.exp.IR.go$GOTABLE,
                      GOSMPL = sup.exp.IR.go$GOSMPLTABLE)
    
    SaveExcel(neo.IA.df, paste0(path.lib, "/Newly-expressed-genes-A-genome.xlsx"))
    gc()
    SaveExcel(sup.IA.df, paste0(path.lib, "/Suppresed-genes-A-genome.xlsx"))
    gc()
    SaveExcel(neo.IR.df, paste0(path.lib, "/Newly-expressed-genes-R-genome.xlsx"))
    gc()
    SaveExcel(sup.IR.df, paste0(path.lib, "/Suppresed-genes-R-genome.xlsx"))
    gc()
    
    invisible(ca)
}

CalcNeoExpressedGenes <- function(ca) {
    exp.AA <- ca$EXP$genes$AA
    exp.RR <- ca$EXP$genes$RR
    exp.IA <- rownames(ca$fpkm$IA[apply(ca$fpkm$IA, 1, mean) > FPKM_CUTOFF, ])
    exp.IR <- rownames(ca$fpkm$IR[apply(ca$fpkm$IR, 1, mean) > FPKM_CUTOFF, ])
    exp.all <- ca$EXP$genes$ALL
    
    neo.exp.IA.from.AA <- setdiff(exp.IA, exp.AA)
    neo.exp.IA.from.RR <- intersect(neo.exp.IA.from.AA, exp.RR)
    neo.exp.IR.from.RR <- setdiff(exp.IR, exp.RR)
    neo.exp.IR.from.AA <- intersect(neo.exp.IR.from.RR, exp.AA)
    exp.IA.and.AA <- intersect(exp.AA, exp.IA)
    exp.IR.and.RR <- intersect(exp.RR, exp.IR)
    
    rep.exp.IA.from.AA <- setdiff(exp.AA, exp.IA)
    rep.exp.IR.from.RR <- setdiff(exp.RR, exp.IR)
    
    png(paste0(path.lib, "/ExpNeo-genes-venn.png"), 500, 500)
    v <- venn(list(AA = exp.AA, IA = exp.IA, IR = exp.IR, RR = exp.RR))
    dev.off()
    
    genes <- gofull <- gosimple <- vector("list", length = length(attr(v, "intersections")) + 1)
    names(genes) <- names(gofull) <- names(gosimple) <- c(names(attr(v, "intersections")), "II")
    for (vn in names(attr(v, "intersections"))) {
        .g <- attr(v, "intersections")[[vn]]
        e <- try(.go <- DoGO2(sig.genes = .g, all.genes = exp.all))
        if (class(e) == "try-error") {
            gofull[[vn]] <- data.frame(Desc = "no results.")
            gosimple[[vn]] <- data.frame(Desc = "no results.")
        } else {
            gofull[[vn]] <- .go$GOTABLE
            gosimple[[vn]] <- .go$GOSMPLTABLE
        }
        genes[[vn]] <- BindTairName(.g)
    }
    
    .go <- DoGO2(sig.genes = exp.ii.union, all.genes = exp.all)
    genes[["II"]] <- BindTairName(exp.ii.union)
    gofull[["II"]] <- .go$GOTABLE
    gosimple[["II"]] <- .go$GOSMPLTABLE
    
    names(gofull)    <- gsub(":", "-", names(gofull))
    names(gosimple) <- gsub(":", "-", names(gosimple))
    names(genes)     <- gsub(":", "-", names(genes))
    SaveExcel(gofull, paste0(path.lib, "/ExpNeo-genes-GO-full.xlsx"))
    gc()
    SaveExcel(gosimple, paste0(path.lib, "/ExpNeo-genes-GO-simple.xlsx"))
    gc()
    SaveExcel(genes, paste0(path.lib, "/ExpNeo-genes.xlsx"))
    gc() 
    
    invisible(ca)
}


IdentifyVEH <- function(ca) {
    .calchighcvgenes <- function(f, tag) {
        log10.f <- log10(f + 1)
        sd.log10.f <- apply(log10.f, 1, sd, na.rm = T)
        mu.log10.f <- apply(log10.f, 1, mean, na.rm = T)
        cv.log10.f <- sd.log10.f / mu.log10.f
        
        g <- ggplot(data.frame(cv = cv.log10.f, mean = mu.log10.f, gene = ifelse( (cv.log10.f > 0.20) & (mu.log10.f > 1.0), 'VEG', 'nonVEG')),
                    aes(x = mean, y = cv, color = gene))
        g <- g + geom_point(alpha = 0.2, size = 2) + scale_color_manual(values =  c('#999999', '#E41A1C'))
        g <- g + xlab('mean(log10FPKM)') + ylab('cv(log10FPKM)')
        png(paste0(path.expbias, '/CV.disp.', tag, '.png'), 500, 400)
        plot(g)
        dev.off()
        
        highcv <- (cv.log10.f > 0.20) & (mu.log10.f > 1.0)
        highcv.genes <- names(highcv)[highcv]
        goobj <- DoGO2(sig.genes = highcv.genes, all.genes = ca$EXP$genes$ALL)
        SaveExcel(list(full=goobj$GOTABLE, simlify=goobj$GOSMPLTABLE),
                  file.name = paste0(path.expbias, "/GO_", tag, "_CV0.20_and_mu1.0.xls"))
        invisible(list(go = goobj$GOSMPLTABLE, gene = highcv.genes,
                       df = data.frame(gene = highcv.genes, sd = sd.log10.f[highcv.genes],
                                  mean = mu.log10.f[highcv.genes], cv = cv.log10.f[highcv.genes])))
    }
    
    A <- ca$fpkm$AA[ca$EXP$genes$ALL, ]
    R <- ca$fpkm$RR[ca$EXP$genes$ALL, ]
    I <- ca$fpkm$II[ca$EXP$genes$ALL, ]
    a <- ca$fpkm$IA[ca$EXP$genes$ALL, ]
    r <- ca$fpkm$IR[ca$EXP$genes$ALL, ]
    A.go <- .calchighcvgenes(A, "AA")
    R.go <- .calchighcvgenes(R, "RR")
    a.go <- .calchighcvgenes(a, "IA")
    r.go <- .calchighcvgenes(r, "IR")
    I.go <- .calchighcvgenes(I, "II")
    print(paste0("# highCV genes in A: ", length(A.go$gene)))
    print(paste0("# highCV genes in R: ", length(R.go$gene)))
    print(paste0("# highCV genes in a: ", length(a.go$gene)))
    print(paste0("# highCV genes in r: ", length(r.go$gene)))
    print(paste0("# highCV genes in I: ", length(I.go$gene)))
    
    vlist <- list(IA = a.go$gene, IR = r.go$gene)#), A = A.go$gene, R = R.go$gene)
    
    highcv.genes.df <- list(Sheet1 = BindTairName(a.go$df),
                            Sheet2 = BindTairName(r.go$df),
                            Sheet3 = BindTairName(A.go$df),
                            Sheet4 = BindTairName(R.go$df))
    highcv.genes.df <- list(IA = BindTairName(a.go$df),
                            IR = BindTairName(r.go$df),
                            AA = BindTairName(A.go$df),
                            RR = BindTairName(R.go$df))
    SaveExcel(highcv.genes.df, file.name = paste0(path.expbias, "/VEHs.IA_IR_AA_RR.xlsx"))
 
    goterms <- unique(c(a.go$go$ID[a.go$go$qvalue < 0.1],
                        r.go$go$ID[r.go$go$qvalue < 0.1],
                        A.go$go$ID[A.go$go$qvalue < 0.1],
                        R.go$go$ID[R.go$go$qvalue < 0.1]))
    go2desc <- rep("", length = length(goterms))
    gomat <- matrix(NA, ncol = 4, nrow = length(goterms))

    rownames(gomat) <- names(go2desc) <- goterms
    colnames(gomat) <- c("IA", "IR", "A", "R")
    go2desc[a.go$go$ID[a.go$go$qvalue < 0.1]] <- a.go$go$Description[a.go$go$qvalue < 0.1]
    go2desc[r.go$go$ID[r.go$go$qvalue < 0.1]] <- r.go$go$Description[r.go$go$qvalue < 0.1]
    go2desc[A.go$go$ID[A.go$go$qvalue < 0.1]] <- A.go$go$Description[A.go$go$qvalue < 0.1]
    go2desc[R.go$go$ID[R.go$go$qvalue < 0.1]] <- R.go$go$Description[R.go$go$qvalue < 0.1]
    
    tpm <- a.go$go$qvalue; names(tpm) <- a.go$go$ID
    gomat[, "IA"] <- -log10(tpm[rownames(gomat)])
    tpm <- r.go$go$qvalue; names(tpm) <- r.go$go$ID
    gomat[, "IR"] <- -log10(tpm[rownames(gomat)])
    tpm <- A.go$go$qvalue; names(tpm) <- A.go$go$ID
    gomat[, "A"] <- -log10(tpm[rownames(gomat)])
    tpm <- R.go$go$qvalue; names(tpm) <- R.go$go$ID
    gomat[, "R"] <- -log10(tpm[rownames(gomat)])
    
    #gomat[a.go$go$ID[a.go$go$qvalue < 0.1], "IA"] <- -log10(a.go$go$qvalue[a.go$go$qvalue < 0.1])
    #gomat[r.go$go$ID[r.go$go$qvalue < 0.1], "IR"] <- -log10(r.go$go$qvalue[r.go$go$qvalue < 0.1])
    #gomat[A.go$go$ID[A.go$go$qvalue < 0.1], "A"] <- -log10(A.go$go$qvalue[A.go$go$qvalue < 0.1])
    #gomat[R.go$go$ID[R.go$go$qvalue < 0.1], "R"] <- -log10(R.go$go$qvalue[R.go$go$qvalue < 0.1])

    df <- data.frame(gomat, desc = go2desc)
    SaveExcel(list(Sheet1 = df), file.name = paste0(path.expbias, "/VEH_GOresults.xlsx"))
    
    dfrawgo <- list(VEH_AA = A.go$go, VEH_IA = a.go$go, VEH_IR = r.go$go, VEH_RR = R.go$go)
    SaveExcel(dfrawgo, file.name = paste0(path.expbias, "/VEH_GOresults_AA_IA_IR_RR.xlsx"))
        

    ca$VEHGO <- gomat
    ca$VEH <- highcv.genes.df
    ca
}


plotFCRsquared <- function(ca) {
    A <- log10(ca$fpkm$AA + 1)[ca$EXP$genes$ALL, ]
    R <- log10(ca$fpkm$RR + 1)[ca$EXP$genes$ALL, ]
    a <- log10(ca$fpkm$IA + 1)[ca$EXP$genes$ALL, ]
    r <- log10(ca$fpkm$IR + 1)[ca$EXP$genes$ALL, ]
    A.fc <- (A - A[, 1])[, -1]
    R.fc <- (R - R[, 1])[, -1]
    a.fc <- (a - a[, 1])[, -1]
    r.fc <- (r - r[, 1])[, -1]
    
    f <- log10(ca$fpkm$II + 1)[ca$EXP$genes$ALL, ]
    sd.log10.f <- apply(f, 1, sd, na.rm = T)
    mu.log10.f <- apply(f, 1, mean, na.rm = T)
    cv.log10.f <- sd.log10.f / mu.log10.f
    
    fcmat <- cbind(A.fc, a.fc, r.fc, R.fc)
    fcmat <- fcmat[cv.log10.f > 0.6, ]
    colnames(fcmat) <- paste0(rep(c("A", "IA", "IR", "R"), each = 8),
                              " (", TIMELABELS[-1], ")")
    
    rsqmatrix <- matrix(NA, ncol = ncol(fcmat), nrow = ncol(fcmat))
    colnames(rsqmatrix) <- rownames(rsqmatrix) <- colnames(fcmat)
    for (ci in 1:ncol(rsqmatrix)) {
        for (ri in 1:nrow(rsqmatrix)) {
            rsqmatrix[ri, ci] <- summary(lm(fcmat[, ri] ~ fcmat[, ci]))$r.squared
        }
    }
    colnames(rsqmatrix) <- rownames(rsqmatrix) <- rep(TIMELABELS[-1], times = 4)

    png(paste0(path.de, "/log10FC.fpkm.Rsquared.png"), 800, 800)
    heatmap.2(rsqmatrix, trace = "none", scale="none", dendrogram="none", Colv = F, Rowv = F,
              density.info = "none", col = COLSFUNC(), margins = c(4, 4),
              key.xlab = "R-squared value", key.title = "", sepcolor = "white",
              colsep = c(1:4 * 8), rowsep = c(1:4 * 8), sepwidth = c(0.1, 0.1))
    dev.off()
}



plotScatter <- function(ca, g, path.prefix) {
    for (.g in g) {
        .t <- CARHR2TAIR[[.g]]
        .s <- CARHR2NAME[[.g]]
        if (is.null(.t)) {
            gname <- .g
        } else {
            if (is.null(.s)) {.s <- ""
            } else { .s <- paste0(" / ", .s)}
            gname <- paste0(.g, " (", .t, .s, ")")
            
        }
        exp.A <- log10(ca$fpkm$AA[.g, ] + 1)
        exp.a <- log10(ca$fpkm$IA[.g, ] + 1)
        exp.r <- log10(ca$fpkm$IR[.g, ] + 1)
        exp.R <- log10(ca$fpkm$RR[.g, ] + 1)
        exp.I <- log10(ca$fpkm$IA[.g, ] * (1/2) + ca$fpkm$IA[.g, ] * (2/2) + 1)
        ARatio <- ca$aorigin.ratio[.g, ]
        lb <- gsub(" hr", "", TIMELABELS)
        png(paste0(path.prefix, "/", .g, ".png"), 730, 320)
        plot(0, 0, type = "n",  xlim = c(-0.2, 1.2), axes = F,
             ylim = c(-0.2, max(3, max(c(exp.a, exp.r, exp.A, exp.R)))) + 0.2,
             xlab = "A-origin ratio", ylab = "log10FPKM",
             main = gname)
        grid()
        axis(1, c(0, 0.2, 1/3, 0.4, 0.6, 0.8, 1.0), c("0", "0.2", "1/3", "0.4", "0.6", "0.8", "1.0"))
        axis(2)
        #abline(v = c(0, 1/3, 2/3, 1), col = "#6f6f6f", lty = 2)
        abline(v = c(0, 1), col = "#000000", lty = 1)
        text(rev(-0.02 * 1:9), exp.R, labels = lb, col = brewer.pal(8, "Dark2"))
        lines(rev(-0.02 * 1:9), exp.R, col = "#999999")
        text(1 + 0.02*1:9, exp.A, labels = lb, col = brewer.pal(8, "Dark2"))
        lines(1 + 0.02 * 1:9, exp.A, col = "#999999")
        
        lines(ARatio, exp.I, col = "#999999", lwd = 1.4)
        text(ARatio, exp.I, labels = lb, col = brewer.pal(8, "Dark2"), cex = 2)
        dev.off()
    }
    

}






DoGO <- function(sig.genes = NULL, all.genes = NULL, p.values = NULL, ontology = "BP",
                 p.cutoff = 1, q.cutoff = 0.1, method = "over-representation") {
    all.genes <- as.character(all.genes)
    sig.genes <- as.character(sig.genes)
    all.genes.at <- as.character(unlist(CARHR2TAIR[all.genes]))
    sig.genes.at <- as.character(unlist(CARHR2TAIR[sig.genes]))
    all.genes.at <- all.genes.at[all.genes.at != ""]
    sig.genes.at <- sig.genes.at[sig.genes.at != ""]
    ego <- enrichGO(gene = sig.genes.at, universe = all.genes.at, OrgDb = org.At.tair.db,
                    ont = ontology, pAdjustMethod = "BH", keyType = "TAIR",
                    pvalueCutoff = p.cutoff, qvalueCutoff = q.cutoff, readable = FALSE)
    rgo <- as.data.frame(ego)
    list(TABLE = rgo, OBJ = ego)
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








PlotRsqHeatmap <- function(f, type = "all") {
    A <- log10(ca$fpkm$AA + 1)
    R <- log10(ca$fpkm$RR + 1)
    a <- log10(ca$fpkm$IA + 1)
    r <- log10(ca$fpkm$IR + 1)
    y <- cbind(A, a, r, R)
    
    I <- log10(ca$fpkm$II + 1)
    cv.I <- apply(I, 1, sd) / apply(I, 1, mean)
    
    y <- y[cv.I > 0.8, ]
    rsqmatrix <- NULL
    rsqmatrix <- matrix(NA, ncol = 4 * 9, nrow = 4 * 9)
    colnames(rsqmatrix) <- rownames(rsqmatrix) <- paste0(rep(SPLABELS, each = 9), " (", rep(TIMELABELS, times = 4), ")")
    for (ci in 1:ncol(rsqmatrix)) {
        for (ri in 1:nrow(rsqmatrix)) {
            rsqmatrix[ri, ci] <- summary(lm(y[, ri] ~ y[, ci]))$r.squared
        }
    }
    heatmap.2(rsqmatrix, trace = "none", scale="none", dendrogram="none", Colv = F, Rowv = F,
              density.info = "none", col = COLSFUNCC2(), margins = c(10, 10),
              key.xlab = "R-squared value", key.title = "", sepcolor = "white",
              colsep = c(1:4 * 9), rowsep = c(1:4 * 9), sepwidth = c(0.1, 0.1))
}








plotGeneExpTimecourse <- function(ca, g) {
    gexp <- rbind(ca$fpkm$AA[g, ], ca$fpkm$RR[g, ], ca$fpkm$IA[g, ], ca$fpkm$IR[g, ])
    rownames(gexp) <- c('C. amara', 'C. rivularis', 'IA', 'IR')
    gexpdf <- melt(log10(gexp + 1))
    colnames(gexpdf) <- c('sample', 'timepoint', 'log10FPKM')
    
    gexpg <- ggplot(gexpdf, aes(x = timepoint, y = log10FPKM, group = sample, color = sample))
    gexpg <- gexpg + geom_line() + scale_color_nejm() + ylim(0, 3)
    gexpg <- gexpg + xlab('time point') + ylab('log10(FPKM)')
    gexpg <- gexpg + theme_bw() + guides(color = guide_legend(title = ''))
    gexpg
    
    rexpdf <- data.frame(Aorigin = ca$aorigin.ratio[g, ], timepoint = names(ca$aorigin.ratio[g, ]), gene = 'gene')
    rexpg <- ggplot(rexpdf, aes(x = timepoint, y = Aorigin, group = gene))
    rexpg <- rexpg + geom_line()# + scale_color_nejm()
    rexpg <- rexpg + xlab('time point') + ylab('A-origin ratio')
    rexpg <- rexpg + theme_bw() + ylim(0, 1)
    rexpg
 
    list(exp = gexpg, ratio = rexpg)   
}



## save GO genes into Excel
annotVEH <- function(glist) {
    VEH.A <- (glist %in% ca$VEH$AA$gene)
    VEH.a <- (glist %in% ca$VEH$IA$gene)
    VEH.r <- (glist %in% ca$VEH$IR$gene)
    VEH.R <- (glist %in% ca$VEH$RR$gene)
    df <- data.frame(VEH_AA = VEH.A, VEH_IA = VEH.a, VEH_IR = VEH.r, VEH_RR = VEH.R)
    df
}
 

Plot.VEH.GO.Profile <- function(ca, goid) {
    go.enriched.genes <- GO2CARHR[[goid]]
    veh.genes.A <- ca$VEH$AA$gene
    veh.genes.a <- ca$VEH$IA$gene
    veh.genes.r <- ca$VEH$IR$gene
    veh.genes.R <- ca$VEH$RR$gene
    
    gene.name.list <- list(
        AA = intersect(go.enriched.genes, veh.genes.A),
        IA = intersect(go.enriched.genes, veh.genes.a),
        IR = intersect(go.enriched.genes, veh.genes.r),
        RR = intersect(go.enriched.genes, veh.genes.R)
    )
    venn(gene.name.list)
    
    gene.fpkm.list <- NULL
    for (catename in names(gene.name.list)) {
        gene.fpkm.list[[catename]] <- ca$fpkm[[catename]][gene.name.list[[catename]], ]
        gene.name.x <- rownames(gene.fpkm.list[[catename]])
        gene.name.x.df <- BindTairName(gene.name.x)
        gene.name.x <- gene.name.x.df$name
        rownames(gene.fpkm.list[[catename]]) <- gene.name.x
    }
    
    gene.fc.list <- NULL
    all.gene.set <- NULL
    for (catename in names(gene.name.list)) {
        gene.fpkm.x <- log2(gene.fpkm.list[[catename]] + 1)
        gene.fc.list[[catename]] <- gene.fpkm.x[, -1] - gene.fpkm.x[, 1]
        all.gene.set <- c(all.gene.set, rownames(gene.fc.list[[catename]]))
    }
    all.gene.set <- unique(all.gene.set)
    all.gene.sheet <- matrix(NA, ncol = 4 * 8, nrow = length(all.gene.set))
    rownames(all.gene.sheet) <- all.gene.set
    all.gene.sheet[rownames(gene.fc.list[['AA']]), 1:8] <- gene.fc.list[['AA']]
    all.gene.sheet[rownames(gene.fc.list[['IA']]), 9:16] <- gene.fc.list[['IA']]
    all.gene.sheet[rownames(gene.fc.list[['IR']]), 17:24] <- gene.fc.list[['IR']]
    all.gene.sheet[rownames(gene.fc.list[['RR']]), 25:32] <- gene.fc.list[['RR']]

    rownames(all.gene.sheet) <- all.gene.set
    SaveExcel(list(Sheet1=all.gene.sheet),
              file = paste0(path.de, '/gene.prof.diff.go', gsub(':', '', goid) , '.xlsx'))
    
    
    all.gene.set <- unique(unlist(gene.name.list))
    fcsheet <- matrix(NA, ncol = 4 * 8, nrow = length(all.gene.set))
    f.A <- log2(ca$fpkm$AA[all.gene.set, ] + 1)
    f.a <- log2(ca$fpkm$IA[all.gene.set, ] + 1)
    f.r <- log2(ca$fpkm$IR[all.gene.set, ] + 1)
    f.R <- log2(ca$fpkm$RR[all.gene.set, ] + 1)
    fc.A <- f.A[, -1] - f.A[, 1]
    fc.a <- f.a[, -1] - f.a[, 1]
    fc.r <- f.r[, -1] - f.r[, 1]
    fc.R <- f.R[, -1] - f.R[, 1]



    fcsheet <- cbind(fc.A[all.gene.set, ], fc.a[all.gene.set, ], fc.r[all.gene.set, ], fc.R[all.gene.set, ])
    fcsheet <- BindTairName(fcsheet)
    fcsheet <- cbind(fcsheet, annotVEH(rownames(fcsheet)))
    SaveExcel(list(Sheet1=fcsheet),
              file = paste0(path.de, '/gene.prof.diff.go', gsub(':', '', goid) , '.alldata.xlsx'))
    
    matfc <- cbind(fc.A, fc.a, fc.r, fc.R)
    colnames(matfc) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 8), ' (',colnames(matfc), ')')
    rownames(matfc) <- fcsheet$name
    png(paste0(path.de, '/gene.prof.diff.go', gsub(':', '', goid) , '.alldata.png'), 2000, 2200, res=220)
    heatmap.2(as.matrix(matfc), trace = "none", scale = "none", Rowv = T, Colv = T,
              hclustfun = function(x) hclust(x,method = 'ward.D2'), lhei = c(8, 30),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", key.title = '',
              cexRow = 0.7, cexCol = 0.7)
    dev.off()

}



MeltToFC <- function(mat) {
    mat <- log2(mat + 1)
    mat.A <- mat[, 2:9] - mat[, 1]
    mat.a <- mat[, 11:18] - mat[, 10]
    mat.r <- mat[, 20:27] - mat[, 19]
    mat.R <- mat[, 29:36] - mat[, 28]
    
    matfc <- cbind(mat.A, mat.a, mat.r, mat.R)
    matfc
}



.common.RData  <- paste0(PATH$lib, "/src/common.RData")
.ca.RData <- paste0(PATH$lib, "/src/ca.RData")
load(.common.RData)
load(.ca.RData)
stop()


if (file.exists(.ca.RData)) {
    load(.ca.RData)
} else {
    design <- GetDesign()
    ca <- NULL
    ca <- GetCounts(ca, design)
    ca <- CalcFPKM(ca)
    ca <- CalcExpressedGenes(ca)
    #ca <- CalcNeoExpressedGenes(ca)
    ca <- IdentifyVEH(ca)
    ca <- PlotAoriginRatioDist(ca)
    save(ca, design, file = .ca.RData)
}



## A-origin ratio datasheet
write.table(BindTairName(ca$aorigin.ratio), sep = '\t', quote = F, row.names = T,
            file = paste0(path.expbias, '/A.origin.ratio.all_homeologs.xls'))



## scatter plots between homeologs at each time points
for (i in 1:9) {
    fa <- ca$fpkm$IA[, i]
    fr <- ca$fpkm$IR[, i]
    keep <- (fa > 1) | (fr > 1)
    df <- data.frame(gene = names(keep)[keep], A = log10(fa[keep] + 1), R = log10(fr[keep] + 1))
    g <- ggplot(df, aes(x = A, y = R))
    g <- g + geom_point(alpha = .3) + coord_fixed() + theme_bw()
    g <- g + xlab('log10(IA-FPKM)') + ylab('log10(IR-FPKM)')
    pdf(paste0(path.expbias, "/FPKM-cinuseta-scatter-", TIMELABELS[i], ".pdf"), 4, 4)
    print(g)
    dev.off()
}
cci <- rep(NA, 9)
for (i in 1:9) {
    xa <- ca$count$IA[, i]
    xr <- ca$count$IR[, i]
    fa <- ca$count$IA[, i]
    fr <- ca$count$IR[, i]
    keep <- (fa > 1) | (fr > 1)
    df <- data.frame(gene = names(keep)[keep], A = log10(xa[keep] + 1), R = log10(xr[keep] + 1))
    g <- ggplot(df, aes(x = A, y = R))
    g <- g + geom_point(alpha = .3) + coord_fixed() + theme_bw()
    g <- g + xlab('log10(count) of A-origin reads') + ylab('log10(count) of R-origin reads')
    g <- g + geom_abline(slope = 1, intercept = log10(2), color = '#DF8F44')
    pdf(paste0(path.expbias, "/count-cinuseta-scatter-", TIMELABELS[i], ".pdf"), 4, 4)
    print(g)
    dev.off()
    cci[i] <- cor(df$A, df$R)
}

## scatter plots between homeologs at each time points
ccp <- rep(NA, 9)
for (i in 1:9) {
    fa <- ca$fpkm$AA[, i]
    fr <- ca$fpkm$RR[, i]
    keep <- (fa > 1) | (fr > 1)
    df <- data.frame(gene = names(keep)[keep], A = log10(fa[keep] + 1), R = log10(fr[keep] + 1))
    g <- ggplot(df, aes(x = A, y = R))
    g <- g + geom_point(alpha = .3) + coord_fixed() + theme_bw()
    g <- g + xlab('log10(C. amara FPKM)') + ylab('log10(C. rivularis FPKM)')
    g <- g + xlim(0, 5) + ylim(0, 5) 
    pdf(paste0(path.expbias, "/FPKM-parent-scatter-", TIMELABELS[i], ".pdf"), 4, 4)
    print(g)
    dev.off()
    ccp[i] <- cor(df$A, df$R)
}




## scatterplot of A-origin ratio between any two samples
cc <- matrix(NA, ncol = 9, nrow = 9)
colnames(cc) <- rownames(cc) <- TIMELABELS
for (i in 1:8) {
    for (j in (i + 1):9) {
        fa <- ca$fpkm$IA[, c(i, j)]
        fr <- ca$fpkm$IR[, c(i, j)]
        keep <- (rowSums(cbind(fa, fr) > 1) > 0)
        ardf <- data.frame(samplei = ca$aorigin.ratio[keep, i],
                           samplej = ca$aorigin.ratio[keep, j])
        ardf <- ardf[rowSums(is.na(ardf)) == 0, ]
        g <- ggplot(ardf, aes(samplei, samplej))
        g <- g + geom_point(size = 1, alpha = 0.5)
        g <- g + theme_bw()
        g <- g + xlab(TIMELABELS[i]) + ylab(TIMELABELS[j])
        pdf(paste0(path.expbias, "/", TIMELABELS[i], "-", TIMELABELS[j],
                   "-aorigin-ratio.scatter.pdf"), 4.5, 4.5)
        print(ggMarginal(g, type = "histogram", bins = 40, col = "#363636", fill ="#363636"))
        dev.off()
        cc[i, j] <- cor(ardf$samplei, ardf$samplej, method = 'pearson')
    }
}





## number of newly-expressed genes and silenced genes
lapply(ca$EXP$genes, length)
length(intersect(ca$EXP$genes$AA, ca$EXP$genes$IA)) / length(union(ca$EXP$genes$AA, ca$EXP$genes$IA)) # A - Ia
length(intersect(ca$EXP$genes$RR, ca$EXP$genes$IR)) / length(union(ca$EXP$genes$RR, ca$EXP$genes$IR)) # R - Ir
length(intersect(ca$EXP$genes$AA, ca$EXP$genes$RR)) / length(union(ca$EXP$genes$AA, ca$EXP$genes$RR)) # A - R
length(intersect(ca$EXP$genes$IA, ca$EXP$genes$IR)) / length(union(ca$EXP$genes$IA, ca$EXP$genes$IR)) # Ia - Ir
length(setdiff(ca$EXP$genes$AA, ca$EXP$genes$IA))
length(setdiff(ca$EXP$genes$IA, ca$EXP$genes$AA))
length(setdiff(ca$EXP$genes$RR, ca$EXP$genes$IR))
length(setdiff(ca$EXP$genes$IR, ca$EXP$genes$RR))

exp.AA <- ca$EXP$genes$AA
exp.IA <- ca$EXP$genes$IA
exp.RR <- ca$EXP$genes$RR
exp.IR <- ca$EXP$genes$IR

neo.IA.go <- DoGO(sig.genes = setdiff(exp.IA, exp.AA), all.genes = union(exp.IA, exp.AA))
sil.IA.go <- DoGO(sig.genes = setdiff(exp.AA, exp.IA), all.genes = union(exp.IA, exp.AA))
neo.IR.go <- DoGO(sig.genes = setdiff(exp.IR, exp.RR), all.genes = union(exp.IR, exp.RR))
sil.IR.go <- DoGO(sig.genes = setdiff(exp.RR, exp.IR), all.genes = union(exp.IR, exp.RR))
write.table(neo.IA.go, file = 'result_files/neoA.expressed.genes.go.xls', sep = '\t', col.names = TRUE)
write.table(neo.IR.go, file = 'result_files/neoR.expressed.genes.go.xls', sep = '\t', col.names = TRUE)
write.table(sil.IA.go, file = 'result_files/silenceA.expressed.genes.go.xls', sep = '\t', col.names = TRUE)
write.table(sil.IR.go, file = 'result_files/silenceR.expressed.genes.go.xls', sep = '\t', col.names = TRUE)






## number of newly-expressed genes and silenced genes (single replicate)
for (i in 1:9) {
    exp.AA <- ca$fpkm$AA[, i] > FPKM_CUTOFF
    exp.RR <- ca$fpkm$RR[, i] > FPKM_CUTOFF
    exp.IA <- ca$fpkm$IA[, i] > FPKM_CUTOFF
    exp.IR <- ca$fpkm$IR[, i] > FPKM_CUTOFF
    exp.AA <- names(exp.AA)[exp.AA]
    exp.RR <- names(exp.RR)[exp.RR]
    exp.IA <- names(exp.IA)[exp.IA]
    exp.IR <- names(exp.IR)[exp.IR]
    
    print('========================================')
    print(i)
    print(paste0('A : ', length(exp.AA)))
    print(paste0('R : ', length(exp.RR)))
    print(paste0('IA: ', length(exp.IA)))
    print(paste0('IR: ', length(exp.IR)))
    print(paste0('A -> IA (-) : ', length(setdiff(exp.AA, exp.IA))))
    print(paste0('A -> IA (+) : ', length(setdiff(exp.IA, exp.AA))))
    print(paste0('R -> IR (-) : ', length(setdiff(exp.RR, exp.IR))))
    print(paste0('R -> IR (+) : ', length(setdiff(exp.IR, exp.RR))))
    print(paste0('Overlap(A, Ia):  ', length(intersect(exp.AA, exp.IA)) / length(union(exp.AA, exp.IA))))  # A - Ia
    print(paste0('Overlap(R, Ir):  ', length(intersect(exp.RR, exp.IR)) / length(union(exp.RR, exp.IR))))  # R - Ir
    print(paste0('Overlap(A, R) :  ', length(intersect(exp.AA, exp.RR)) / length(union(exp.AA, exp.RR))))  # A - R
    print(paste0('Overlap(Ia, Ir): ', length(intersect(exp.IA, exp.IR)) / length(union(exp.IA, exp.IR))))  # Ia - Ir
}


## venn diagrams of VEHs
veh.set <- list(A = ca$VEH$AA$gene, R = ca$VEH$RR$gene, IA = ca$VEH$IA$gene, IR = ca$VEH$IR$gene)
venn.diagram(veh.set, fill = pal_nejm()(4), alpha = rep(0.6, 4),
             height = 1600, width = 1600, main.cex = 0.8,
             filename = paste0(path.de, '/VEH.set.tiff'))







## gene expression profiles (z-scored)
g <- c('CARHR049560', # STM/BUM
       'CARHR279690', # BAM1
       'CARHR158740', # BAM2
       'CARHR227100', # BAM3
       'CARHR202780', # PLT3
       'CARHR274370', # REV
       'CARHR137920') # PDF1
g <- rev(g)
g.fpkm <- cbind(ca$fpkm$AA[g, ], ca$fpkm$IA[g, ], ca$fpkm$IR[g, ], ca$fpkm$RR[g, ])
rownames(g.fpkm) <- rev(c('STM/BUM', 'BAM1', 'BAM2', 'BAM3', 'PLT3', 'REV', 'PDF1'))
colnames(g.fpkm) <- paste0(rep(c('C. amara', 'IA', 'IR', 'C. rivularis'), each = 9),
                           '  ', colnames(g.fpkm))
.g.fpkm <- t(apply(log10(g.fpkm + 1), 1, scale))
colnames(.g.fpkm) <- colnames(g.fpkm)
g.fpkm <- .g.fpkm
g.fpkm.df <- melt(g.fpkm)
colnames(g.fpkm.df) <- c('gene', 'library', 'zscore')
ghm <- ggplot(g.fpkm.df, aes(x = library, y = gene, fill = zscore))
ghm <- ghm + geom_tile() + scale_fill_gradientn('value',
            colours = rev(brewer.pal(9, 'Spectral')))
ghm <- ghm + xlab('') + ylab('') + coord_fixed()
ghm <- ghm + theme_bw() + guides(fill = guide_legend(title = 'z-score'))
ghm <- ghm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   strip.background = element_rect(fill = "white", colour = "white"))
pdf(paste0(path.geneprof, '/meristem.selected.homeologs.pdf'), 12.3, 3.6)
plot(ghm)
dev.off()






## meristem associated genes (gene expression)
meristem.genes <- unique(c(GO2CARHR[['GO:0035266']], GO2CARHR[['GO:0010075']], GO2CARHR[['GO:0048509']]))
meristem.genes <- unique(GO2CARHR[['GO:0048507']])
meristem.genes <- intersect(meristem.genes, rownames(ca$fpkm$AA))
meristem.fpkm <- cbind(ca$fpkm$AA[meristem.genes, ], ca$fpkm$IA[meristem.genes, ], 
                      ca$fpkm$IR[meristem.genes, ], ca$fpkm$RR[meristem.genes, ])
colnames(meristem.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                           '  ', colnames(meristem.fpkm))
meristem.fpkm <- log10(meristem.fpkm + 1)
meristem.fpkm <- BindTairName(meristem.fpkm)

## plot water deprivation
g.fpkm <- meristem.fpkm[, 1:36]
g.fpkm.names <- meristem.fpkm[, 38]
names(g.fpkm.names) <- rownames(g.fpkm)
g.fpkm.names['CARHR027940'] <- 'HSL1(1)'
g.fpkm.names['CARHR239910'] <- 'HSL1(2)'
rownames(g.fpkm) <- g.fpkm.names
pdf(paste0(path.geneprof, '/meristem.homeologs.pdf'), 13, 45)
heatmap.2(as.matrix(g.fpkm), trace = "none", scale = "row", dendrogram = "row", Colv = F,
          hclustfun = function(x) hclust(x,method = 'ward.D2'), lhei = c(1, 20),
          density.info = "none", col = COLSFUNC2(), key.xlab = "z-score",
          colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
dev.off()








## water deprivation
watdep.genes <- GO2CARHR[['GO:0009414']]
watdep.genes <- intersect(watdep.genes, rownames(ca$fpkm$AA))
watdep.fpkm <- cbind(ca$fpkm$AA[watdep.genes, ], ca$fpkm$IA[watdep.genes, ], 
                     ca$fpkm$IR[watdep.genes, ], ca$fpkm$RR[watdep.genes, ])
colnames(watdep.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                                '  ', colnames(watdep.fpkm))
watdep.fpkm <- log10(watdep.fpkm + 1)
watdep.fpkm <- BindTairName(watdep.fpkm)


## plot water deprivation
g.fpkm <- watdep.fpkm[, 1:36]
rownames(g.fpkm) <- watdep.fpkm[, 38]
pdf(paste0(path.geneprof, '/water.deprivation.homeologs.pdf'), 13, 36)
heatmap.2(as.matrix(g.fpkm), trace = "none", scale = "row", dendrogram = "row", Colv = F,
          hclustfun = function(x) hclust(x,method = 'ward.D2'), lhei = c(1, 18),
          density.info = "none", col = COLSFUNC2(), key.xlab = "z-score",
          colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
dev.off()





## 	alcohol metabolic process
alcohol.genes <- GO2CARHR[['GO:0006066']]
alcohol.genes <- intersect(alcohol.genes, rownames(ca$fpkm$AA))
alcohol.fpkm <- cbind(ca$fpkm$AA[alcohol.genes, ], ca$fpkm$IA[alcohol.genes, ], 
                     ca$fpkm$IR[alcohol.genes, ], ca$fpkm$RR[alcohol.genes, ])
colnames(alcohol.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                                '  ', colnames(alcohol.fpkm))
alcohol.fpkm <- log10(alcohol.fpkm + 1)
alcohol.fpkm <- BindTairName(alcohol.fpkm)


## plot 	alcohol metabolic process
g.fpkm <- alcohol.fpkm[, 1:36]
rownames(g.fpkm) <- alcohol.fpkm[, 38]
pdf(paste0(path.geneprof, '/alcohol_metabolic_processs.homeologs.pdf'), 13, 36)
heatmap.2(as.matrix(g.fpkm), trace = "none", scale = "row", dendrogram = "row", Colv = F,
          hclustfun = function(x) hclust(x,method = 'ward.D2'), lhei = c(1, 18),
          density.info = "none", col = COLSFUNC2(), key.xlab = "z-score",
          colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
dev.off()








## response to ethylene
ethylene.genes <- GO2CARHR[['GO:0009723']]
ethylene.genes <- intersect(ethylene.genes, rownames(ca$fpkm$AA))
ethylene.fpkm  <- cbind(ca$fpkm$AA[ethylene.genes, ], ca$fpkm$IA[ethylene.genes, ], 
                        ca$fpkm$IR[ethylene.genes, ], ca$fpkm$RR[ethylene.genes, ])
colnames(ethylene.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                                '  ', colnames(ethylene.fpkm))
ethylene.fpkm <- log10(ethylene.fpkm + 1)
ethylene.fpkm <- BindTairName(ethylene.fpkm)





meristem.fpkm <- cbind(meristem.fpkm, annotVEH(rownames(meristem.fpkm)))
watdep.fpkm <- cbind(watdep.fpkm, annotVEH(rownames(watdep.fpkm)))
ethylene.fpkm <- cbind(ethylene.fpkm, annotVEH(rownames(ethylene.fpkm)))
alcohol.fpkm <- cbind(alcohol.fpkm, annotVEH(rownames(alcohol.fpkm)))

SaveExcel(list(meristem = meristem.fpkm, waterdep = watdep.fpkm, ethylene = ethylene.fpkm,
               alcoholmetabo = alcohol.fpkm),
          file = paste0(path.geneprof, '/gene.prof.xls'))




## save all genes into Excel
allexpgene.fpkm <- cbind(ca$fpkm$AA, ca$fpkm$IA, 
                         ca$fpkm$IR, ca$fpkm$RR)
colnames(allexpgene.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                                '  ', colnames(allexpgene.fpkm))
allexpgene.fpkm <- log10(allexpgene.fpkm + 1)
allexpgene.fpkm <- BindTairName(allexpgene.fpkm)

SaveExcel(list(all.exp.genes = allexpgene.fpkm),
          file = paste0(path.geneprof, '/gene.prof.all.xls'))







erf1.exp <- plotGeneExpTimecourse(ca, 'CARHR099190') # ERF1
cca1.exp <- plotGeneExpTimecourse(ca, 'CARHR142180') # CCA1
pdf1.exp <- plotGeneExpTimecourse(ca, 'CARHR137920') # PDF1
stp6.exp <- plotGeneExpTimecourse(ca, 'CARHR080460') # STP6

pdf(paste0(path.geneprof, '/ERF1.geneexp.pdf'), 6, 2.5)
plot(erf1.exp$exp)
dev.off()
pdf(paste0(path.geneprof, '/ERF1.ratio.pdf'), 5, 2.1)
plot(erf1.exp$ratio)
dev.off()

pdf(paste0(path.geneprof, '/CCA1.geneexp.pdf'), 6, 2.5)
plot(cca1.exp$exp)
dev.off()
pdf(paste0(path.geneprof, '/CCA1.ratio.pdf'), 5, 2.1)
plot(cca1.exp$ratio)
dev.off()

pdf(paste0(path.geneprof, '/PDF1.geneexp.pdf'), 6, 2.5)
plot(pdf1.exp$exp)
dev.off()
pdf(paste0(path.geneprof, '/PDF1.ratio.pdf'), 5, 2.1)
plot(pdf1.exp$ratio)
dev.off()


pdf(paste0(path.geneprof, '/STP6.geneexp.pdf'), 6, 2.5)
plot(stp6.exp$exp)
dev.off()
pdf(paste0(path.geneprof, '/STP61.ratio.pdf'), 5, 2.1)
plot(stp6.exp$ratio)
dev.off()





golist <- list(
    A = rownames(ca$VEHGO)[!is.na(ca$VEHGO[, "A"])],
    R = rownames(ca$VEHGO)[!is.na(ca$VEHGO[, "R"])],
    IA = rownames(ca$VEHGO)[!is.na(ca$VEHGO[, "IA"])],
    IR = rownames(ca$VEHGO)[!is.na(ca$VEHGO[, "IR"])]
)

venn.diagram(golist, fill = pal_nejm()(4), alpha = rep(0.6, 4),
             height = 1600, width = 1600, main.cex = 1,
             filename = paste0(path.geneprof, '/govenn.tiff'))




Plot.VEH.GO.Profile(ca, 'GO:0009414') # water deprivation
Plot.VEH.GO.Profile(ca, 'GO:0006066') # alchol
Plot.VEH.GO.Profile(ca, 'GO:0009723') # ethylene
Plot.VEH.GO.Profile(ca, 'GO:0009738') # abscisic acid-activated signaling pathway
Plot.VEH.GO.Profile(ca, 'GO:0006979') # response to ABS

Plot.VEH.GO.Profile(ca, 'GO:0009737') # response to ABS
Plot.VEH.GO.Profile(ca, 'GO:0071215') #
Plot.VEH.GO.Profile(ca, 'GO:0009873') # ehtylene activate pathway
Plot.VEH.GO.Profile(ca, 'GO:0009688') # ABS acid biosynthetic process



# LEC1, STM, LAS, BBM/PLT4, CUC1, CUC2, CUC3
fl <- c('CARHR022320', 'CARHR049560', 'CARHR045280', 'CARHR195780', 'CARHR090150',
        'CARHR267090', 'CARHR096710', 'CARHR041230', 'CARHR202780', 'CARHR094250',
        'CARHR245020', 'CARHR279490', 'CARHR259260', 'CARHR063720', 'CARHR188260',
        'CARHR202230', 'CARHR239940', 'CARHR023820', 'CARHR048650', 'CARHR015300',
        'CARHR066610', 'CARHR274370', 'CARHR004640', 'CARHR131950', 'CARHR069980',
        'CARHR066610', 'CARHR064200', 'CARHR143940', 'CARHR023460',
        'CARHR246030', 'CARHR164140', 'CARHR085820',
        'CARHR099190', 'CARHR142180', 'CARHR137920', # ERF1, CCA1, PDF1 
        'CARHR128870', 'CARHR029900', 'CARHR069240', 'CARHR057050', 'CARHR279690',
        'CARHR158740', 'CARHR200150', 'CARHR093320', 'CARHR285700', 'CARHR036610',
        'CARHR085820', 'CARHR128040', 'CARHR021160', 'CARHR077710', 'CARHR195360',
        'CARHR243890', 'CARHR021150')
        

dir.create('result_files/timecourse/meristem', recursive = TRUE)
.getlabel <- function(g) {
    caa <- ifelse(g %in% ca$VEH$AA$gene, 'is', 'isnot')
    cia <- ifelse(g %in% ca$VEH$IA$gene, 'is', 'isnot')
    cir <- ifelse(g %in% ca$VEH$IR$gene, 'is', 'isnot')
    crr <- ifelse(g %in% ca$VEH$RR$gene, 'is', 'isnot')
    list(labels = c('C.~amara', 'I[A]', 'I[R]', 'C.~rivularis'),
         colors = c(caa, cia, cir, crr))
}

for (.fl in fl) {
    .flb <- .getlabel(.fl)
    if (.fl %in% rownames(ca$fpkm$AA)) {
    fAA <- log10(ca$fpkm$AA[.fl, ] + 1)
    fIA <- log10(ca$fpkm$IA[.fl, ] + 1)
    fRR <- log10(ca$fpkm$RR[.fl, ] + 1)
    fIR <- log10(ca$fpkm$IR[.fl, ] + 1)
    xdf <- data.frame(
        fpkm = c(fAA, fIA, fIR, fRR),
        #time = c(names(fAA), names(fIA), names(fIR), names(fRR)),
        time = rep(c(0, 2, 4, 8, 12, 24, 48, 72, 96), times = 4),
        species = rep(.flb$labels, each = length(fAA)),
        VEH = rep(.flb$colors, each = length(fAA))
    )
    xdf$species <- factor(xdf$species, levels = .flb$labels)
    xdf$VEH <- factor(xdf$VEH, levels = c('isnot', 'is'))
    g <- ggplot(xdf, aes(x = time, y = fpkm, group = species, colour = VEH))
    g <- g + geom_line(size=2)
    g <- g + scale_x_continuous(breaks = c(0,  4,  12, 24, 48, 72, 96))
    g <- g + facet_grid(. ~ species, labeller = label_parsed)
    g <- g + ggtitle(paste0(.fl, ' / ', CARHR2TAIR[[.fl]], ' (', CARHR2NAME[[.fl]], ')'))
    g <- g + xlab('hours after submergence') + ylab('log10(FPKM + 1)')
	g <- g + theme(legend.text = element_text(size = 38), plot.title = element_text(size = 40),
                   axis.text.x = element_text(size = 32),
                   axis.text.y = element_text(size = 32),
                   axis.title.x = element_text(size = 38),
                   axis.title.y = element_text(size = 38),
                   strip.text.x = element_text(size = 38),
                   legend.position = 'none')
    cols <- c('#1b1919', '#cc0c00')
    if (length(unique(xdf$VEH)) == 1) {
        cols <- ifelse(unique(xdf$VEH) == 'is', '#cc0c00', '#1b1919')
    }
    g <- g + scale_color_manual(values = cols)
    png(paste0('result_files/timecourse/meristem/', .fl, '.png'), 2600, 480)
    print(g)
    dev.off()
    }
}






fl <- c('CARHR163320', 'CARHR171770', 'CARHR187720', 'CARHR285680',
        'CARHR271420', 'CARHR064080', 'CARHR066510', 'CARHR131280')
dir.create('result_files/timecourse/waterdep', recursive = TRUE)
for (.fl in fl) {
    .flb <- .getlabel(.fl)
    fAA <- log10(ca$fpkm$AA[.fl, ] + 1)
    fIA <- log10(ca$fpkm$IA[.fl, ] + 1)
    fRR <- log10(ca$fpkm$RR[.fl, ] + 1)
    fIR <- log10(ca$fpkm$IR[.fl, ] + 1)
    xdf <- data.frame(
        fpkm = c(fAA, fIA, fIR, fRR),
        #time = c(names(fAA), names(fIA), names(fIR), names(fRR)),
        time = rep(c(0, 2, 4, 8, 12, 24, 48, 72, 96), times = 4),
        species = rep(.flb$labels, each = length(fAA)),
        VEH = rep(.flb$colors, each = length(fAA))
    )
    xdf$species <- factor(xdf$species, levels = .flb$labels)
    xdf$VEH <- factor(xdf$VEH, levels = c('isnot', 'is'))
    g <- ggplot(xdf, aes(x = time, y = fpkm, group = species, colour = VEH))
    g <- g + geom_line(size=2)
    g <- g + scale_x_continuous(breaks = c(0,  4,  12, 24, 48, 72, 96))
    g <- g + facet_grid(. ~ species, labeller = label_parsed)
    g <- g + ggtitle(paste0(.fl, ' / ', CARHR2TAIR[[.fl]], ' (', CARHR2NAME[[.fl]], ')'))
    g <- g + xlab('hours after submergence') + ylab('log10(FPKM + 1)')
	g <- g + theme(legend.text = element_text(size = 38), plot.title = element_text(size = 40),
                   axis.text.x = element_text(size = 32),
                   axis.text.y = element_text(size = 32),
                   axis.title.x = element_text(size = 38),
                   axis.title.y = element_text(size = 38),
                   strip.text.x = element_text(size = 38),
                   legend.position = 'none')
    g <- g + scale_color_manual(values = c('#1b1919', '#cc0c00'))
    png(paste0('result_files/timecourse/waterdep/', .fl, '.png'), 2600, 480)
    print(g)
    dev.off()
}












