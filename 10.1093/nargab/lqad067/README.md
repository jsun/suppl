# A low-coverage 3' RNA-seq to detect homeolog expression in polyploid wheat


## Expression Quantification with RNA-Seq Data

1. See `quant_shells/genome_index/README.md` to check how to prepare reference sequence of Chinese Scpring (IWGSC v1.1) and how to identify homeolog triads.
2. See `quant_shells/fullseq/README.md` to check the analysis processes (e.g., QC, read mapping, read counting) of the conventional paired-end RNA-seq data.
3. See `quant_shells/tagseq/README.md` to check the analysis processes of 3' RNA-seq data.


## Data summarization

To summarize the gene and homeolog expression and perform the downstream analysis, run the following R scripts in this directory.

- `analyze_gff.R`
- `analyze_seqlendist.R`
- `analyze_correlation.R`
- `analyze_downsampling.R`
- `bam_overlap.R`
- `analyze_deg.R`


