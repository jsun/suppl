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

Use `calc_nmappedreads.sh` instead of the following scripts.


```bash
qsub calc_nmappedreads.sh
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


