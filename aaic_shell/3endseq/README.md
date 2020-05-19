# Quantification of homeolog expression

FASTQ quality control and homeolog expression quantification.

```
qsub 01_qc.sh
qsub 02_quant_hisat.sh
qsub 03_count_hisat.sh 
qsub 04_eagle.sh
```



# Commands for summarizing statistics

Count the number of reads in FASTQ file with `wc` command.
The output should be divided by 4.

```
gzip -dc input.fastq.gz | wc -l
```


Count the number of reads that uniquely mapped on A, B, and D-subgenomes
with HISAT2 by using `samtools`.

```
log_fpath=n_reads_ABDgenome.tsv

for bam_fpath in `ls bam | grep -v csi`
do
     echo ${bam_fpath} >> ${log_fpath}
     samtools view -h bam/${bam_fpath} | grep 'chr.A' | grep 'NH:i:1' | wc -l >> ${log_fpath}
     samtools view -h bam/${bam_fpath} | grep 'chr.B' | grep 'NH:i:1' | wc -l >> ${log_fpath}
     samtools view -h bam/${bam_fpath} | grep 'chr.D' | grep 'NH:i:1' | wc -l >> ${log_fpath}
done
```

Count the number of read pairs that uniquely mapped on A, B, and D-subgenomes
with HISAT2 by using `samtools`.
Note that, the output should be divided by 2 for counting the pairs.



