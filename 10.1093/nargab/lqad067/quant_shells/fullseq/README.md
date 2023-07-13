# paired-end RNA-Seq

## Quantification

```bash
qsub 01_qc.sh
qsub 02_quant_hisat.sh
qsub 03_quant_star.sh
qsub 04_quant_eaglerc.sh
```

```bash
qsub shortage_fastq.sh
cd ~/projects/HuoberBrezel/data
mkdir -p shortfullseq/clean_fastq
mv fullseq/clean_fastq/*.short.fastq.gz  shortfullseq/clean_fastq/

qsub 12_quant_hisat.sh
```



## Summarizing of mapping statistics


```bash
qsub calc_nmappedreads.sh
```


## Downsampling

```bash
qsub downsampling.sh
```




