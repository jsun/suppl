library(tidyverse)
library(ggsci)



plot_dist <- function(len_list) {
    df <- NULL
    for (seqname in names(len_list)) {
        df <- rbind(df, data.frame(replicate = seqname, seqlen = len_list[[seqname]][, 1]))
    }
    
    g <- ggplot(df, aes(x = seqlen)) +
            geom_histogram(binwidth = 1) + scale_y_log10() +
            xlab('read length') + ylab('frequency') +
            facet_wrap(~ replicate, ncol = 3)
    
    g
}




#
# seq length distribution of 3'-end RNA-Seq
#
tagseq_len <- vector('list', 6)
tagseq_len[[1]] <- read.table('aaic_data/tagseq_seqlen/20181109.A-TaeRS_1_Tae_RS1_1.cleaned.fastq.gz.len.tsv.gz', 
                           sep = '\t', header = FALSE)
tagseq_len[[2]] <- read.table('aaic_data/tagseq_seqlen/20181109.A-TaeRS_1_Tae_RS1_2.cleaned.fastq.gz.len.tsv.gz', 
                           sep = '\t', header = FALSE)
tagseq_len[[3]] <- read.table('aaic_data/tagseq_seqlen/20181109.A-TaeRS_1_Tae_RS1_3.cleaned.fastq.gz.len.tsv.gz', 
                           sep = '\t', header = FALSE)
tagseq_len[[4]] <- read.table('aaic_data/tagseq_seqlen/20181109.A-TaeRS_1_Tae_RS1_16.cleaned.fastq.gz.len.tsv.gz', 
                           sep = '\t', header = FALSE)
tagseq_len[[5]] <- read.table('aaic_data/tagseq_seqlen/20181109.A-TaeRS_1_Tae_RS1_17.cleaned.fastq.gz.len.tsv.gz', 
                           sep = '\t', header = FALSE)
tagseq_len[[6]] <- read.table('aaic_data/tagseq_seqlen/20181109.A-TaeRS_1_Tae_RS1_18.cleaned.fastq.gz.len.tsv.gz', 
                           sep = '\t', header = FALSE)


prob <- matrix(NA, ncol = 6, nrow = 110 - 39)
rownames(prob) <- as.character(40:110)

for (i in seq(tagseq_len)) {
    x <- tagseq_len[[i]][, 1]
    freq <- table(x)
    prob[names(freq), i] <- freq / sum(freq)
}

probdf <- data.frame(length = rownames(prob), prob = rowMeans(prob))

# the weights used for shortening paired-end RNA-Seq
print(probdf)

names(tagseq_len) <- paste0('replicate ', 1:6)


#
# seq length distribution of (original) paired-end RNA-Seq
#
fullseq_len <- vector('list', 6)
fullseq_len[[1]] <- read.table('aaic_data/fullseq_seqlen/20181221.A-ZH_W2017_1_CS_2_R1.cleaned.fastq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseq_len[[2]] <- read.table('aaic_data/fullseq_seqlen/20181221.A-ZH_W2017_1_CS_3_R1.cleaned.fastq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseq_len[[3]] <- read.table('aaic_data/fullseq_seqlen/20181221.A-ZH_W2017_1_CS_4_R1.cleaned.fastq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseq_len[[4]] <- read.table('aaic_data/fullseq_seqlen/20181221.A-ZH_W2017_1_CS_cold_2_R1.cleaned.fastq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseq_len[[5]] <- read.table('aaic_data/fullseq_seqlen/20181221.A-ZH_W2017_1_CS_cold_3_R1.cleaned.fastq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseq_len[[6]] <- read.table('aaic_data/fullseq_seqlen/20181221.A-ZH_W2017_1_CS_cold_4_R1.cleaned.fastq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)


names(fullseq_len) <- paste0('replicate ', 1:6)


#
# seq length distribution of (short) paired-end RNA-Seq
#
fullseqs_len <- vector('list', 6)
fullseqs_len[[1]] <- read.table('aaic_data/fullseqshortR1_seqlen/20181221.A-ZH_W2017_1_CS_2_R1.cleaned.short.fq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseqs_len[[2]] <- read.table('aaic_data/fullseqshortR1_seqlen/20181221.A-ZH_W2017_1_CS_3_R1.cleaned.short.fq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseqs_len[[3]] <- read.table('aaic_data/fullseqshortR1_seqlen/20181221.A-ZH_W2017_1_CS_4_R1.cleaned.short.fq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseqs_len[[4]] <- read.table('aaic_data/fullseqshortR1_seqlen/20181221.A-ZH_W2017_1_CS_cold_2_R1.cleaned.short.fq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseqs_len[[5]] <- read.table('aaic_data/fullseqshortR1_seqlen/20181221.A-ZH_W2017_1_CS_cold_3_R1.cleaned.short.fq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)
fullseqs_len[[6]] <- read.table('aaic_data/fullseqshortR1_seqlen/20181221.A-ZH_W2017_1_CS_cold_4_R1.cleaned.short.fq.gz.len.tsv.gz',
                                sep = '\t', header = FALSE)

names(fullseqs_len) <- paste0('replicate ', 1:6)








tagdist <- plot_dist(tagseq_len)
fulldist <- plot_dist(fullseq_len)
fullsdist <- plot_dist(fullseqs_len)


png('results/plots/seqLen_dist_tagseq.png', 1400, 800, res = 220)
print(tagdist)
dev.off()

png('results/plots/seqLen_dist_fullseq.png', 1400, 800, res = 220)
print(fulldist)
dev.off()

png('results/plots/seqLen_dist_fullseqshortR1.png', 1400, 800, res = 220)
print(fullsdist)
dev.off()



tag_multimap <- c(0.737244873, 0.734866085, 0.707786078, 0.757859627, 0.745625696, 0.72090912)
fulls_multimap <- c(0.693053258, 0.705167695, 0.696864113, 0.739947523, 0.73099401, 0.71668429)
fulls_multimap <- c(0.693422829, 0.701584505, 0.694187239, 0.731994172, 0.721804586, 0.711665759)
t.test(tag_multimap, fulls_multimap, alternative = 'two.sided', paired = FALSE, var.equal = TRUE)






