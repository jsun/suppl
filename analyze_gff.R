library(tidyverse)
library(ggsci)
options(stringsAsFactors = FALSE)





DWCODE <- c('iwgsc', '0.5k', '1.0k', '2.0k', '3.0k', '4.0k')
DWCODE_FMT <- c('IWGSC', '+0.5k', '+1k', '+2k', '+3k', '+4k')



# Obtained this numbers from calculating of cleaned FASTQ directly.
# Since featureCounts double counts the multiple mapped reads,
# we cannot just calculate the total reads from the log.
#                 RS1_1    RS1_16   RS1_17   RS1_18   RS1_2    RS1_3
INPUT_FQSIZE <- c(2078810, 1393116, 3385387, 3217619, 1241292, 2874064)




W <- 680
H <- 460

plot_bar <- function(df) {
    df$annotation <- factor(df$annotation, levels = DWCODE_FMT)
    dfsum <- df %>% group_by(annotation) %>% dplyr::summarise(mean = mean(value), sd = sd(value))

    gp <- ggplot() +
            geom_line(aes(x = annotation, y = mean, group = 1), data = dfsum) +
            geom_errorbar(aes(x = annotation, ymin = mean - sd, ymax = mean + sd, width = 0.3), data = dfsum) +
            geom_point(aes(x = annotation, y = value), data = df) +
            ylab('') + xlab('')
    gp
}

plot_log <- function() {
    
    tcode <- c('Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity')
    
    assigned_df <-  nofeature_df <- ambig_df <- data.frame(annotation = NULL, sample = NULL, value = NULL)

    for (wi in 1:length(DWCODE)) {
        x <- read.table(paste0('counts/counts_', DWCODE[wi], '/all.counts.gene.tsv.summary.gz'),
                        sep = '\t', header = TRUE)
        y <- x
        y[, -1] <- sweep(y[, -1], 2, 100 / INPUT_FQSIZE, '*')

        x <- x[x[, 1] %in% tcode, ]
        y <- y[y[, 1] %in% tcode, ]
        # print(DWCODE_FMT[di])      # print for Excel data
        # print(t(x[, -1]))  # print for Excel data
        
        assigned_df <- rbind(assigned_df,
              data.frame(annotation = DWCODE_FMT[wi], sample = colnames(y[, -1]), value = as.numeric(y[1, -1])))
        nofeature_df <- rbind(nofeature_df,
              data.frame(annotation = DWCODE_FMT[wi], sample = colnames(y[, -1]), value = as.numeric(y[2, -1])))
        ambig_df <- rbind(ambig_df,
              data.frame(annotation = DWCODE_FMT[wi], sample = colnames(y[, -1]), value = as.numeric(y[3, -1])))
    }
    
    assigned_df %>% select(annotation, value) %>%
            group_by(annotation) %>% summarise(mean = mean(value), sd = sd(value))
    nofeature_df %>% select(annotation, value) %>%
            group_by(annotation) %>% summarise(mean = mean(value), sd = sd(value))
    ambig_df %>% select(annotation, value) %>%
            group_by(annotation) %>% summarise(mean = mean(value), sd = sd(value))
    
    png(paste0('results/plots/assigned_reads.png'), W, H, res = 220)
    print(plot_bar(assigned_df))
    dev.off()
    
    png(paste0('results/plots/nofeature_reads.png'), W, H, res = 220)
    print(plot_bar(nofeature_df))
    dev.off()
    
    png(paste0('results/plots/ambig_reads.png'), W, H, res = 220)
    print(plot_bar(ambig_df))
    dev.off()
}


plot_gene <- function() {
 
    tcode <- c('Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity')
    
    gene_df <- data.frame(annotation = NULL, sample = NULL, value = NULL)
    for (wi in 1:length(DWCODE)) {
        x <- read.table(paste0('counts/counts_', DWCODE[wi], '/all.counts.gene.tsv.gz'),
                        sep = '\t', header = TRUE)
        x <- x[, -1]
        ngenes <- as.numeric(colSums(x > 0))
        gene_df <- rbind(gene_df,
              data.frame(annotation = DWCODE_FMT[wi], sample = colnames(x), value = ngenes))
    }
    
    gene_df %>% select(annotation, value) %>%
            group_by(annotation) %>% summarise(mean = mean(value), sd = sd(value))
    
    png(paste0('results/plots/expressed_genes.png'), W, H, res = 220)
    print(plot_bar(gene_df))
    dev.off()
 
    
    # use CPM
    gene_df <- data.frame(annotation = NULL, sample = NULL, value = NULL)
    for (wi in 1:length(DWCODE)) {
        x <- read.table(paste0('counts/counts_', DWCODE[wi], '/all.counts.gene.tsv.gz'),
                        sep = '\t', header = TRUE)
        x <- x[, -1]
        x <- sweep(x, 2, 1e6 / colSums(x), '*')
        ngenes <- as.numeric(colSums(x > 1))
        gene_df <- rbind(gene_df,
              data.frame(annotation = DWCODE_FMT[wi], sample = colnames(x), value = ngenes))
    }
    
 
   
}




find_gene <- function() {
    
    
    for (wi in 1:length(DWCODE)) {

        x0k <- read.table(paste0('counts/counts_iwgsc/all.counts.gene.tsv.gz'),
                        sep = '\t', header = TRUE)
        x1k <- read.table(paste0('counts/counts_', DWCODE[wi], '/all.counts.gene.tsv.gz'),
                        sep = '\t', header = TRUE)
    
        rownames(x0k) <- rownames(x1k) <- x1k[, 1]
        x0k <- as.matrix(x0k[, -1])[, c(1, 5, 6, 2, 3, 4)]
        x1k <- as.matrix(x1k[, -1])[, c(1, 5, 6, 2, 3, 4)]
    
        colnames(x0k) <- colnames(x1k) <- c('control #1', 'control #2', 'control #3',
                                            'cold #1', 'cold #2', 'cold #3')
        dfcnt <- data.frame(lib_name = NULL, x = NULL, y = NULL)
        for (i in 1:ncol(x0k)) {
            dfcnt <- rbind(dfcnt, data.frame(lib_name = colnames(x0k)[i],
                           x = log10(x0k[, i] + 1), y = log10(x1k[, i] + 1)))
        }
        dfcnt$lib_name <- factor(dfcnt$lib_name, levels = colnames(x0k))
        gp <- ggplot(dfcnt, aes(x = x, y = y)) +
                geom_point() + facet_wrap(~ lib_name, ncol = 3) +
                coord_fixed() + 
                xlab(paste0(DWCODE_FMT[1], ' log10(count)')) +
                ylab(paste0(DWCODE_FMT[1], ' log10(count)'))
        
        png(paste0('results/plots/count_relations_iwgsc_vs_', DWCODE[wi], '.png'), 1800, 1200, res = 220)
        print(gp)
        dev.off()
        
        if (wi == 3) {
            gp <- dfcnt %>% filter(lib_name == 'control #1') %>%
                        ggplot(aes(x = x, y = y)) +
                        geom_point() +
                        coord_fixed() + 
                        xlab(paste0(DWCODE_FMT[1], ' log10(count)')) +
                        ylab(paste0('IWGSC', DWCODE_FMT[wi], ' log10(count)'))
            png(paste0('results/plots/count_relations_iwgsc_vs_', DWCODE[wi], '_control1.png'), 600, 650, res = 220)
            print(gp)
            dev.off()
            
            df <- data.frame(x = log10(1 + x0k[, 1]), y = log10(1 + x1k[, 1]),
                             x2 = x0k[, 1], y2 = x1k[, 2])
            df[df$x > 1 & df$x < 2 & df$y > 3, ]
            # TraesCS4D02G145400 chr4D:134546507-134549866
            df[df$x == 0 & df$y > 2.5, ]
            # TraesCS2B02G330500 chr2B:474899493-474901687
            df[df$y == 0 & df$x > 2, ]
            # TraesCS4D02G263300 chr4D:434272089-434277206
            df[df$y < df$x, ]
            
            df <- df[!(df$x == 0 & df$y == 0), ]
            nrow(df[df$y < df$x, ]) / nrow(df)
            nrow(df[df$y > df$x, ]) / nrow(df)
            nrow(df[df$y == df$x, ]) / nrow(df)
        }
    }
    


}




plot_log()
plot_gene()
find_gene()




