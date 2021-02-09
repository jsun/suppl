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
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrA.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrB.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrD.list | grep -c REF >> ${log_fpath}
done


cd tagseq/eaglercngi

log_fpath=nreads_mapped_stats.tsv
if [ -f ${log_fpath} ]; then
    rm ${log_fpath}
fi

for list_dpath in `ls -d *`; do
    echo ${list_dpath} >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrA.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrB.list | grep -c REF >> ${log_fpath}
    grep -v REVERSE ${list_dpath}/${list_dpath}.ref.chrD.list | grep -c REF >> ${log_fpath}
done

```





To investigate the overlaps between HISAT2 and EGALE-RC results, use the following commands.


```bash
cd ~/projects/HuoberBrezel/data/tagseq/bam

for bam_fpath in `ls *.bam`
do
    samtools view -S ${bam_fpath} | grep "NH:i:1" | grep "chr.A" | cut -f1 > ${bam_fpath%.bam}.chrA.read_name.txt
    samtools view -S ${bam_fpath} | grep "NH:i:1" | grep "chr.B" | cut -f1 > ${bam_fpath%.bam}.chrB.read_name.txt
    samtools view -S ${bam_fpath} | grep "NH:i:1" | grep "chr.D" | cut -f1 > ${bam_fpath%.bam}.chrD.read_name.txt
done



cd ~/projects/HuoberBrezel/data/tagseq/eaglerc

for bam_dpath in `ls | grep -v nreads`
do
    cd ${bam_dpath}
    grep REF ${bam_dpath}.ref.chrA.list | grep -v REVERSE | cut -f1 > ${bam_dpath}.chrA.read_name.txt
    grep REF ${bam_dpath}.ref.chrB.list | grep -v REVERSE | cut -f1 > ${bam_dpath}.chrB.read_name.txt
    grep REF ${bam_dpath}.ref.chrD.list | grep -v REVERSE | cut -f1 > ${bam_dpath}.chrD.read_name.txt
    cd ..
done
```


