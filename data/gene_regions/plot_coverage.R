library(tidyverse)
library(ggsci)


plot_coverage <- function(gene) {
    cov_plus  <- read.table(paste0(gene, '.coverage.plus.bed'), sep = '\t', header = FALSE)
    cov_minus <- read.table(paste0(gene, '.coverage.minus.bed'), sep = '\t', header = FALSE)
    
    df <- rbind(data.frame(pos = cov_plus$V7, coverage = cov_plus$V8, strand = '+'),
                data.frame(pos = cov_minus$V7, coverage = - cov_minus$V8, strand = '-'))
    df$strand <- factor(df$strand, levels = c('+', '-'))
    
    fig <- ggplot(df, aes(x = pos, y = coverage, group = strand, fill = strand)) +
            geom_area() + ylim(min(c(df$coverage, -10)), max(df$coverage)) +
            scale_fill_nejm() + xlab('') + ylab('')
    fig
}

target_genes <- c('TraesCS2B02G330500', 'TraesCS4D02G145400', 'TraesCS4D02G263300')
for (target_gene in target_genes) {
    fig <- plot_coverage(target_gene)
    png(paste0(target_gene, '.png'), 1200, 800, res = 320)
    print(fig)
    dev.off()
}





