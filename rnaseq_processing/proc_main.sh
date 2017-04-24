#!/bin/bash

## Run this script as
##
## $bash proc_main.sh
##
## The variances of RUN_INDEXGENME, RUN_STAR, RUN_READCLASSIFY, RUN_COUNT should
## be 'true' or false.
##
## RUN_INDEXGENME:   Create genome index with STAR.
## RUN_STAR:         Mapping with STAR.
## RUN_READCLASSIFY: Classify read origins.
## RUN_COUNT:        Calculate read counts.
##

RUN_INDEXGENOME=false
RUN_STAR=false
RUN_READCLASSIFY=false
RUN_COUNT=false





## Path settings
##------------------------------------------------------------------------------
source ../local/global.sh
JOB_ROOT=${PROJECT_ROOT}/rnaseq_processing
BAM_DIR=${JOB_ROOT}/results/bam
RC_DIR=${JOB_ROOT}/results/readclassify
COUNT_DIR=${JOB_ROOT}/results/counts
TMP_DIR=${JOB_ROOT}/tmp
mkdir -p ${BAM_DIR}
mkdir -p ${RC_DIR}
mkdir -p ${COUNT_DIR}
mkdir -p ${TMP_DIR}




## Assembly parameter settings
##------------------------------------------------------------------------------
TAGS=(orig common other)
GENOMES=(camara crivularis)




## Create genome index
##------------------------------------------------------------------------------
if [ ${RUN_INDEXGENOME} ]; then
    cd ${GENOMELEAFDIR}/camara
    mkdir -p star_index
    ${BIN}/STAR --runMode genomeGenerate \
                --genomeDir ./star_index \
                --genomeFastaFiles genome.fa \
                --runThreadN ${THREADS} --limitGenomeGenerateRAM 10000000000

    cd ${GENOMELEAFDIR}/crivularis
    mkdir -p star_index
    ${BIN}/STAR --runMode genomeGenerate \
                --genomeDir ./star_index \
                --genomeFastaFiles genome.fa \
                --runThreadN ${THREADS} --limitGenomeGenerateRAM 10000000000
fi




## Mapping
##------------------------------------------------------------------------------
cd ${JOB_ROOT}
if [ ${RUN_STAR} ]; then
    for seq_lib in ${RNASEQ_LIBS[@]}; do
        fq1_file=${DATADIR}/fastq/${seq_lib}_filtered_R1.fq.gz
        fq2_file=${DATADIR}/fastq/${seq_lib}_filtered_R2.fq.gz
        
        for ref_genome in ${GENOMES[@]}; do
            tmp_dir=${TMP_DIR}/star/${seq_lib}__on__${ref_genome}
            mkdir -p ${tmp_dir}
            cd ${tmp_dir}
            
            star_index_path=${GENOMELEAFDIR}/${ref_genome}/star_index
            out_bam_name=${seq_lib}__on__${ref_genome}
            
            ${BIN}/STAR \
                --runThreadN ${THREADS} \
                --genomeDir ${star_index_path} \
                --readFilesIn ${fq1_file} ${fq2_file} \
                --readFilesCommand zcat \
                --genomeLoad NoSharedMemory
            
            ${BIN}/samtools-0.1.19/samtools view -bS Aligned.out.sam > Aligned.out.bam
            ${BIN}/samtools-0.1.19/samtools sort Aligned.out.bam ${out_bam_name}
            ${BIN}/samtools-0.1.19/samtools index ${out_bam_name}.bam
            
            mv ${out_bam_name}.bam ${BAM_DIR}/
            mv ${out_bam_name}.bam.bai ${BAM_DIR}/
            cd -
        done
    done
fi






## Read classification
##------------------------------------------------------------------------------
cd ${JOB_ROOT}
if [ ${RUN_READCLASSIFY} ]; then
    for seq_lib in ${RNASEQ_LIBS[@]}; do
        map_on_amara=${seq_lib}__on__camara
        map_on_rivularis=${seq_lib}__on__crivularis
        
        tmp_dir=${TMP_DIR}/readclassify/${seq_id}
        mkdir -p ${tmp_dir}
        cd ${tmp_dir}
        
        python ${SCRIPTDIR}/rc_filter.py ${BAM_DIR}/${map_on_amara}.bam \
                                         ${tmp_dir}/${map_on_amara}.filtered.bam \
                                         ${tmp_dir}/${map_on_amara}.error.bam


        python ${SCRIPTDIR}/rc_filter.py ${BAM_DIR}/${map_on_rivularis}.bam \
                                         ${tmp_dir}/${map_on_rivularis}.filtered.bam \
                                         ${tmp_dir}/${map_on_rivularis}.error.bam
                    
        python ${SCRIPTDIR}/read_classify.py ${tmp_dir}/${map_on_amara}.filtered.bam \
                                             ${tmp_dir}/${map_on_rivularis}.filtered.bam \
                                             ${tmp_dir}/${map_on_amara}.filtered.rc \
                                             ${tmp_dir}/${map_on_rivularis}.filtered.rc
        
    
        for tag in ${TAGS[@]}; do
            fname=${tmp_dir}/${map_on_amara}.filtered.rc_${tag}
            ${BIN}/samtools-0.1.19/samtools sort ${fname}.bam ${fname}.sorted
            ${BIN}/samtools-0.1.19/samtools index ${fname}.sorted.bam
            mv ${fname}.sorted.bam ${RC_DIR}/
            mv ${fname}.sorted.bam.bai ${RC_DIR}/
            
            fname=${tmp_dir}/${map_on_rivularis}.filtered.rc_${tag}
            ${BIN}/samtools-0.1.19/samtools sort ${fname}.bam ${fname}.sorted
            ${BIN}/samtools-0.1.19/samtools index ${fname}.sorted.bam
            mv ${fname}.sorted.bam ${RC_DIR}/
            mv ${fname}.sorted.bam.bai ${RC_DIR}/
        done
        
    done
fi






## Read counting
##------------------------------------------------------------------------------
cd ${JOB_ROOT}
if [ ${RUN_HTSEQ} ]; then
    tmp_dir=${TMP_DIR}/counts/
    mkdir -p ${tmp_dir}
    cd ${tmp_dir}
    
    for seq_lib in ${RNASEQ_LIBS[@]}; do
        for ref_genome in ${GENOMES[@]}; do
            for tag in ${TAGS[@]}; do
                bam_file=${RC_DIR}/${seq_lib}__on__${ref_genome}.filtered.rc_${tag}.sorted.bam
                gff_file=${GENOMELEAFDIR}/${ref_genome}/genome.modified.gff
                cnt_file=${COUNT_DIR}/${seq_lib}__on__${ref_genome}.filtered.rc_${tag}.sorted.txt

                python ${SCRIPTDIR}/htseq-count2.py -m divided -a 0 -r pos -s reverse \
                                                   -f bam -t gene -i ID \
                                                   ${bam_file} ${gff_file} > ${cnt_file}
            done
        done
    done
fi























