#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=4
#PBS -l elapstim_req=72:00:00
#PBS -N MKIDX_LAST


####   PBS -l memsz_job=128gb


BIN=/home/jqsun/local/bin
cd /home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS

CPU=4

${BIN}/lastdb -uNEAR -R01 index/cds.A_last chrA.cds.fa
${BIN}/lastdb -uNEAR -R01 index/cds.B_last chrB.cds.fa
${BIN}/lastdb -uNEAR -R01 index/cds.D_last chrD.cds.fa

${BIN}/lastal index/cds.A_last -P$CPU -D10000000000 chrB.cds.fa | last-map-probs -m 0.49 > A.B.maf
${BIN}/lastal index/cds.B_last -P$CPU -D10000000000 chrA.cds.fa | last-map-probs -m 0.49 > B.A.maf

${BIN}/lastal index/cds.B_last -P$CPU -D10000000000 chrD.cds.fa | last-map-probs -m 0.49 > B.D.maf
${BIN}/lastal index/cds.D_last -P$CPU -D10000000000 chrB.cds.fa | last-map-probs -m 0.49 > D.B.maf

${BIN}/lastal index/cds.A_last -P$CPU -D10000000000 chrD.cds.fa | last-map-probs -m 0.49 > A.D.maf
${BIN}/lastal index/cds.D_last -P$CPU -D10000000000 chrA.cds.fa | last-map-probs -m 0.49 > D.A.maf


