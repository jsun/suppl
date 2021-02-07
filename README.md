# HuoberBrezel

## RNA-Seq data analysis

To analysis paired-end RNA-Seq data, see `quant_shells/fullseq`.
To anlaysis 3'-end RNA-Seq data, see `quant_shells/tagseq`.


## Data summarization

To summarize RNA-Seq data anlysis results, use R scripts.

- analyze_gff.R
- analyze_seqlendist.R
- analyze_correlation.R
- analyze_downsampling.R
- bam_overlap.R


## Mapping coverage

```
for gene in TraesCS2B02G330500 TraesCS4D02G145400 TraesCS4D02G263300
do
    bedtools coverage \
             -a ${gene}.bed \
             -b /Users/jsun/Desktop/HB/data/tagseq/bam/20181109.A-TaeRS_1_Tae_RS1_1.bam \
             -bed -d -s  > ${gene}.coverage.plus.bed
    
    bedtools coverage \
             -a ${gene}.bed \
             -b /Users/jsun/Desktop/HB/data/tagseq/bam/20181109.A-TaeRS_1_Tae_RS1_1.bam \
             -bed -d -S  > ${gene}.coverage.minus.bed
done

Rscript plot_coverage.R
```



