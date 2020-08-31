# paired-end RNA-Seq

## quantification

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



## summarizing of mapping statistics

This step takes a lot of time since the CS dataset are very large.
Use `calc_nmappedreads.sh` instead of the following scripts.


```bash
qsub calc_nmappedreads.sh
```







