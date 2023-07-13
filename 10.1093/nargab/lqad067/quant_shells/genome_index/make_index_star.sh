#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq-long
#PBS -l gpunum_job=0
#PBS -l cpunum_job=8
#PBS -l memsz_job=500gb
#PBS -l elapstim_req=72:00:00
#PBS -N MKIDX_STAR



cd /home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS


BIN=/home/jqsun/local/bin
dna_fasta=iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
transcript_fasta=iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_transcripts.fasta
gtf=iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gtf


mkdir -p index


# 
# subgenome index (for EAGLERC)
# 
for chr in chrA chrB chrD
do
    mkdir -p index/dna_${chr}_star
    ${BIN}/STAR --runMode genomeGenerate --genomeDir index/dna_${chr}_star \
                                         --genomeFastaFiles ${chr}.fa      \
                                         --sjdbGTFfile ${gtf} --runThreadN 8
done



dna_fasta=refseqv1.0_part/CSv1.0.part.fasta
gff=refseqv1.0_part/IWGSC_v1.1_HC_20170706.ext_1.0k.part.gff3


#
# whole genome index
#
mkdir -p index/dnapart_star
${BIN}/STAR --runMode genomeGenerate --genomeDir index/dnapart_star \
                                     --genomeFastaFiles ${dna_fasta} \
                                     --limitGenomeGenerateRAM 50000000000 \
                                     --runThreadN 8





