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
cd data/gene_regions

for gene in TraesCS2B02G330500 TraesCS4D02G145400 TraesCS4D02G263300 TraesCS7A02G115400 TraesCS7B02G013100
do
    bedtools coverage \
             -a ${gene}.bed \
             -b ~/projects/HuoberBrezel/data/tagseq/bam/TaeRS2728_g172.bam \
             -bed -d -s  > ${gene}.tcs.coverage.plus.bed
    bedtools coverage \
             -a ${gene}.bed \
             -b ~/projects/HuoberBrezel/data/tagseq/bam/TaeRS2728_g172.bam \
             -bed -d -S  > ${gene}.tcs.coverage.minus.bed
    bedtools coverage \
             -a ${gene}.bed \
             -b ~/projects/HuoberBrezel/data/tagseq/bam/20181109.A-TaeRS_1_Tae_RS1_1.bam \
             -bed -d -s  > ${gene}.cs.coverage.plus.bed
    bedtools coverage \
             -a ${gene}.bed \
             -b ~/projects/HuoberBrezel/data/tagseq/bam/20181109.A-TaeRS_1_Tae_RS1_1.bam \
             -bed -d -S  > ${gene}.cs.coverage.minus.bed
done

Rscript plot_coverage.R
```


## Gene annotations

```
cd data
# download iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB.gz from IWGSC
python parse_funcanno_v1.py
```

