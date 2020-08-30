# tagseq


## quantification of tagseq data

Quantification of tagseq data sequenced from CS and TCS samples.

```
qsub 01_qc.sh
qsub 02_quant_hisat.sh
qsub 03_quant_star.sh
qsub 04_quant_eaglerc.sh
```


## downsampling experiment

Downsample CS data for investigating of the number of genes that can be detected by tagseq.


```
qsub 20_dwsample_fastq.sh
qsub 21_dwsample_quant.sh
```




## summarizing of mapping statistics

To count the number of mapped reads mapped on A, B, and D-subgenome
from HISAT2 mapping results, use the `samtools` with the following options.

```bash
cd ~/projects/HuoberBrezel/data
for dpath in `ls -d tagseq*`; do
    cd ${dpath}/bam
    
    log_fpath=nreads_mapped_stats.tsv
    if [ -f ${log_fpath} ]; then
        rm ${log_fpath}
    fi
    
    for bam_fpath in `ls *.bam`; do
        echo ${bam_fpath} >> ${log_fpath}
        samtools view -h ${bam_fpath} | grep 'chr.A' | grep 'NH:i:1' | wc -l >> ${log_fpath}
        samtools view -h ${bam_fpath} | grep 'chr.B' | grep 'NH:i:1' | wc -l >> ${log_fpath}
        samtools view -h ${bam_fpath} | grep 'chr.D' | grep 'NH:i:1' | wc -l >> ${log_fpath}
    done
    
    cd -
done
```

To summarize STAR results, use the following commands.

```bash
cd ~/projects/HuoberBrezel/data
cd tagseq/bamstar

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi
    
for bam_fpath in `ls *.bam`; do
    echo ${bam_fpath} >> ${log_fpath}
    samtools view -h -q 255 ${bam_fpath} | grep 'chr.A' | grep -v "@SQ" | wc -l >> ${log_fpath}
    samtools view -h -q 255 ${bam_fpath} | grep 'chr.B' | grep -v "@SQ" | wc -l >> ${log_fpath}
    samtools view -h -q 255 ${bam_fpath} | grep 'chr.D' | grep -v "@SQ" | wc -l >> ${log_fpath}
done
```


To summarize EAGLERC results, use the following commands.

```bash
cd ~/projects/HuoberBrezel/data
cd tagseq/eaglerc

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for list_dpath in `ls -d *`; do
    echo ${list_dpath} >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrA.list >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrB.list >> ${log_fpath}
    grep -c REF ${list_dpath}/${list_dpath}.ref.chrD.list >> ${log_fpath}
done
```




