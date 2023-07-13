# 3' RNA-Seq


## Expression quantification

```
qsub 01_qc.sh
qsub 02_quant_hisat.sh
qsub 03_quant_star.sh
qsub 04_quant_eaglerc.sh
```


## Downsampling experiment


```
qsub 20_dwsample_fastq.sh
qsub 21_dwsample_quant.sh
```




## Summarizing of mapping statistics


```bash
qsub calc_nmappedreads.sh
```




