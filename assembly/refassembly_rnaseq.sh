#!/bin/bash

## Run this script as
## 
## $bash refassembly_rnaseq.sh camara
## $bash refassembly_rnaseq.sh crivularis
##



## Bash arguments
## Set up 'camara' or 'crivularis' to assembly C. amara or C. rivularis genome.
##------------------------------------------------------------------------------
if [ $# -ne 1 ]; then
    echo 'Set the species name as camara or crivularis'
    exit 1
fi
species=$1




## Assembly parameter settings
##------------------------------------------------------------------------------

MINCOVER=5
MAXCOVER=100000000




## Path settings
##------------------------------------------------------------------------------

source ../local/global.sh
JOB_ROOT=${PROJECT_ROOT}/assembly
SPDIR=${JOB_ROOT}/${species}
chirsuta_genome=${GENOMEDIR}/chirsuta



## Path settings
##------------------------------------------------------------------------------

if [ ${species} = 'camara' ]; then
    RNASEQ_R1=${FASTQDIR}/859_BG_filtered_R1.fq.gz,${FASTQDIR}/859_BN_filtered_R1.fq.gz,${FASTQDIR}/859_BH_filtered_R1.fq.gz,${FASTQDIR}/859_BI_filtered_R1.fq.gz,${FASTQDIR}/859_BJ_filtered_R1.fq.gz,${FASTQDIR}/859_BK_filtered_R1.fq.gz,${FASTQDIR}/859_BL_filtered_R1.fq.gz,${FASTQDIR}/859_BM_filtered_R1.fq.gz,${FASTQDIR}/859_BW_filtered_R1.fq.gz
    RNASEQ_R2=${FASTQDIR}/859_BG_filtered_R2.fq.gz,${FASTQDIR}/859_BN_filtered_R2.fq.gz,${FASTQDIR}/859_BH_filtered_R2.fq.gz,${FASTQDIR}/859_BI_filtered_R2.fq.gz,${FASTQDIR}/859_BJ_filtered_R2.fq.gz,${FASTQDIR}/859_BK_filtered_R2.fq.gz,${FASTQDIR}/859_BL_filtered_R2.fq.gz,${FASTQDIR}/859_BM_filtered_R2.fq.gz,${FASTQDIR}/859_BW_filtered_R2.fq.gz
    #RNASEQ_R1=${FASTQDIR}/test_R1.fq.gz
    #RNASEQ_R2=${FASTQDIR}/test_R2.fq.gz
elif [ ${species} = 'crivularis' ]; then
    RNASEQ_R1=${FASTQDIR}/859_AZ_filtered_R1.fq.gz,${FASTQDIR}/859_AF_filtered_R1.fq.gz,${FASTQDIR}/859_BA_filtered_R1.fq.gz,${FASTQDIR}/859_BB_filtered_R1.fq.gz,${FASTQDIR}/859_BC_filtered_R1.fq.gz,${FASTQDIR}/859_BD_filtered_R1.fq.gz,${FASTQDIR}/859_BE_filtered_R1.fq.gz,${FASTQDIR}/859_BF_filtered_R1.fq.gz,${FASTQDIR}/859_BV_filtered_R1.fq.gz
    RNASEQ_R2=${FASTQDIR}/859_AZ_filtered_R2.fq.gz,${FASTQDIR}/859_AF_filtered_R2.fq.gz,${FASTQDIR}/859_BA_filtered_R2.fq.gz,${FASTQDIR}/859_BB_filtered_R2.fq.gz,${FASTQDIR}/859_BC_filtered_R2.fq.gz,${FASTQDIR}/859_BD_filtered_R2.fq.gz,${FASTQDIR}/859_BE_filtered_R2.fq.gz,${FASTQDIR}/859_BF_filtered_R2.fq.gz,${FASTQDIR}/859_BV_filtered_R2.fq.gz
else
    echo 'Set the species name as camara or crivularis'
    exit 1
fi







mkdir -p ${SPDIR}
cd ${SPDIR}


TMPDIR=${SPDIR}/tmp1
mkdir -p ${TMPDIR}
cd ${TMPDIR}

cp ${chirsuta_genome}/genome.fa genome.fa
cp ${chirsuta_genome}/genome.gff genome.gff




## Main processes
##------------------------------------------------------------------------------


for i in 1 2 3 4 5 6 7 8 9 10; do
    
    j=$(( ${i} + 1 ))
    
    TMPDIR=${SPDIR}/tmp${i}
    NEXT_TMPDIR=${SPDIR}/tmp${j}
    mkdir -p ${TMPDIR}
    mkdir -p ${NEXT_TMPDIR}
    
    
    cd ${TMPDIR}
    
    date '+%Y-%m-%d %H:%M'
    mkdir -p ${TMPDIR}/star_index
    ${BIN}/STAR --runMode genomeGenerate \
                --genomeDir ${TMPDIR}/star_index \
                --genomeFastaFiles genome.fa \
                --runThreadN ${THREADS} --limitGenomeGenerateRAM 10000000000

    date '+%Y-%m-%d %H:%M'
    ${BIN}/STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${TMPDIR}/star_index \
        --readFilesIn ${RNASEQ_R1} ${RNASEQ_R2} \
        --readFilesCommand zcat \
        --genomeLoad NoSharedMemory
    
    ${BIN}/samtools-0.1.19/samtools view -@ ${THREADS} -bS Aligned.out.sam > Aligned.out.bam
    ${BIN}/samtools-0.1.19/samtools sort Aligned.out.bam Aligned.out.sorted
    ${BIN}/samtools-0.1.19/samtools index Aligned.out.sorted.bam
    ${BIN}/samtools-0.1.19/samtools depth Aligned.out.sorted.bam > Aligned.out.sorted.bam.depth

    ${BIN}/samtools-0.1.19/samtools mpileup -uf genome.fa Aligned.out.sorted.bam | ${BIN}/samtools-0.1.19/bcftools/bcftools view -bvcg - > snps_samt.bcf
    ${BIN}/samtools-0.1.19/bcftools/bcftools view snps_samt.bcf > snps_samt.vcf

    python ${SCRIPTDIR}/sortvcf.py snps_samt.vcf > snps_samt.sorted.vcf 
    python ${SCRIPTDIR}/vcffilter.py snps_samt.sorted.vcf ${MINCOVER} ${MAXCOVER} > snps_samt.filtered.vcf
    python ${SCRIPTDIR}/replace.py genome.fa snps_samt.filtered.vcf 1> assembled_genome.fa 2> assembled_genome.hist
    python ${SCRIPTDIR}/add_species_and_depth.py genome.fa Aligned.out.sorted.bam.depth ${MINCOVER} > assembled_genome.depth.fa
    python ${SCRIPTDIR}/replace_fa_gff.py assembled_genome.depth.fa snps_samt.filtered.vcf genome.gff assembled_genome.final.fa assembled_genome.final.gff 2> replacedpy.log

    
    mv assembled_genome.final.fa ${NEXT_TMPDIR}/genome.fa
    mv assembled_genome.final.gff ${NEXT_TMPDIR}/genome.gff
done



cp ${NEXT_TMPDIR}/genome.fa ${GENOMELEAFDIR}/${species}/
cp ${NEXT_TMPDIR}/genome.gff ${GENOMELEAFDIR}/${species}/
    
cd ${GENOMELEAFDIR}/${species}
python ${SCRIPTDIR}/chfmt_gff.py genome.gff genome.modified.gff



cd ${JOB_ROOT}

