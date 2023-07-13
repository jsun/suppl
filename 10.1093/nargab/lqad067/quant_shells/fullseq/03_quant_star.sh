#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=256gb
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -N HB_FUL_QUANT_STAR


nCPU=16


PROJECT_DIR=~/projects/HuoberBrezel
DATA_DIR=${PROJECT_DIR}/data/fullseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bamstar
COUNTS_DIR=${DATA_DIR}/countsstar
GENOME_DIR=~/projects/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_star
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.gff3



cd ${PROJECT_DIR}




# mapping
mkdir -p ${BAM_DIR}
cd ${BAM_DIR}

for fastq_filepath in `ls ${CLEAN_FASTQ_DIR}/*_R1.cleaned.fastq.gz `; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} _R1.cleaned.fastq.gz`
        
    echo 'start mapping ...'
    echo `date`
    time STAR --runThreadN ${nCPU} \
              --genomeDir ${GENOME_INDEX} \
              --readFilesCommand zcat \
              --readFilesIn ${CLEAN_FASTQ_DIR}/${fastq_filename}_R1.cleaned.fastq.gz ${CLEAN_FASTQ_DIR}/${fastq_filename}_R2.cleaned.fastq.gz \
              --genomeLoad NoSharedMemory \
              --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
              --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
              --outFileNamePrefix ${BAM_DIR}/${fastq_filename}_
    echo `date`
    echo 'finished.'

    samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}_Aligned.out.sam
    samtools index -c ${fastq_filename}.bam
    #rm ${fastq_filename}_Aligned.out.sam
    
done






# count mapped reads

mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
featureCounts -p --countReadPairs -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH} \
              -o ${COUNTS_DIR}/counts.gene.ext_1.0k.tsv *.bam
cd ${COUNTS_DIR}
gzip *.tsv
gzip *.tsv.summary





