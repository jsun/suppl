#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=256gb
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -N HB_FUL_QUANT_STARCOUNT

pSTAR=0
pQuant=1
nCPU=16


PROJECT_DIR=~/projects/HuoberBrezel
BIN=~/local/bin
UTILS=~/local/utilfunc
DATA_DIR=${PROJECT_DIR}/data/fullseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bamstar
GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dnapart_star


COUNTS_DIR=${DATA_DIR}/countsstar
GTF_FILEPATH=${GENOME_DIR}/refseqv1.0_part/IWGSC_v1.1_HC_20170706.ext_1.0k.part.gff3





cd ${PROJECT_DIR}



if [ ${pSTAR} -eq 1 ]; then
    

mkdir -p ${BAM_DIR}
cd ${BAM_DIR}

for fastq_filepath in `ls ${CLEAN_FASTQ_DIR}/*_R1.cleaned.fastq.gz `; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} _R1.cleaned.fastq.gz`
        
    echo 'start mapping ...'
    echo `date`
    time ${BIN}/STAR --runThreadN ${nCPU} \
                     --genomeDir ${GENOME_INDEX} \
                     --readFilesCommand zcat \
                     --readFilesIn ${CLEAN_FASTQ_DIR}/${fastq_filename}_R1.cleaned.fastq.gz ${CLEAN_FASTQ_DIR}/${fastq_filename}_R2.cleaned.fastq.gz \
                     --genomeLoad NoSharedMemory \
                     --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                     --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
                     --outFileNamePrefix ${BAM_DIR}/${fastq_filename}_
    
    echo `date`
    echo 'finished.'

    ${BIN}/samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}_Aligned.out.sam
    ${BIN}/samtools index -c ${fastq_filename}.bam
    rm ${fastq_filename}_Aligned.out.sam
    
done
fi






if [ ${pQuant} -eq 1 ]; then
#
# count mapped reads
#
mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
${BIN}/featureCounts -p -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} \
                     -o ${COUNTS_DIR}/counts.gene.ext_1.0k.tsv *.bam
fi





