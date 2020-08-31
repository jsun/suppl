#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -N HB_SFUL_QUANT_HISAT

pHISAT=1
pQuant=1
nCPU=16


PROJECT_DIR=~/projects/HuoberBrezel
BIN=~/local/bin
UTILS=~/local/utilfunc
DATA_DIR=${PROJECT_DIR}/data/shortfullseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bam
GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2


COUNTS_DIR=${DATA_DIR}/counts
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706


cd ${PROJECT_DIR}



if [ ${pHISAT} -eq 1 ]; then
    
mkdir -p ${BAM_DIR}
cd ${CLEAN_FASTQ_DIR}

for fastq_filepath in `ls *.fastq.gz `; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} _R1.cleaned.short.fastq.gz`
    
    time ${BIN}/hisat2 -p ${nCPU} -x ${GENOME_INDEX} \
                       -U ${fastq_filepath} \
                       -S ${BAM_DIR}/${fastq_filename}.sam
        
    cd ${BAM_DIR}
    ${BIN}/samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
    ${BIN}/samtools index -c ${fastq_filename}.bam
    rm ${fastq_filename}.sam
    cd -
done
fi






if [ ${pQuant} -eq 1 ]; then
#
# count mapped reads
#
mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
gff_versions=("iwgsc" "ext_1.0k")
for gff in ${gff_versions[@]}; do
    ${BIN}/featureCounts -p -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH}.${gff}.gff3 \
                         -o ${COUNTS_DIR}/counts.gene.${gff}.tsv *.bam
done
fi





