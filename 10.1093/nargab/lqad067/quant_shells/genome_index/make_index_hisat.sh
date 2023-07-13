#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq-long
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=196gb
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -N MKIDX_HISAT

nCPU=16
cd ~/projects/data/genome/IWGSC_RefSeq_v1.1_CS
mkdir -p index

dna_fasta=iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta
hisat2-build  -p ${nCPU} ${dna_fasta} index/dna_hisat2


