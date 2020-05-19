#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l memsz_job=64gb
#PBS -l elapstim_req=72:00:00
#PBS -N WGENOME_EAGEL_TEST3

nCPU=16
PROJECT_DIR=/home/jqsun/research/wheat_genome

BIN=/home/jqsun/local/bin
UTILS=/home/jqsun/local/utilfunc
EAGLE_SCRIPT_DIR=/home/jqsun/local/src/eagle/scripts

BATCH_DIR=${PROJECT_DIR}/data/eaglerc
BAM_DIR=${BATCH_DIR}/bam
COUNTS_DIR=${BATCH_DIR}/homeolog_counts
UNIQ_COUNTS_DIR=${BATCH_DIR}/uniq_counts

FASTQ_DIR=${PROJECT_DIR}/data/clean_fastq
EAGLE_DIR=${BATCH_DIR}/eagle

GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS/index/dna_hisat2
GENOME_INDEX_PREFIX=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS/index/dna_
GENOME_INDEX_SUFFIX=_star
GTF_A=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrA.ext.gff3
GTF_B=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrB.ext.gff3
GTF_D=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrD.ext.gff3


pSTAR=0
pEAGLE=0
pEAGLERC=0
pCountHomeolog=1
pCountUnique=1


echo "shell: 02_tagseq.eaglerc.sh"
echo "pSTAR: ${pSTAR}"
echo "pEAGLE: ${pEAGLE}"
echo "pEAGLERC: ${pEAGLERC}"
echo "pCountHomeolog: ${pCountHomeolog}"
echo "pCountUnique: ${pCountUnique}"
echo "nCPU: ${nCPU}"
echo "----- start -----"
echo `date`


cd ${PROJECT_DIR}


## mapping with STAR
if [ ${pSTAR} -eq 1 ]; then
    mkdir -p ${BAM_DIR}
    
    cd ${BAM_DIR}
    for fastq_filepath in `ls ${FASTQ_DIR}/*.clean.fastq.gz`; do
        echo ${fastq_filepath}

        fastq_filename=`basename ${fastq_filepath} _R1.clean.fastq.gz`
        
        # mapping
        chr_names=(chrA chrB chrD)
        for chr_name in ${chr_names[@]}; do
            starsam=${BAM_DIR}/${fastq_filename}__on__${chr_name}
            mkdir -p ${starsam}
            cd ${starsam}
            
            echo "map ${fastq_filename} on ${chr_name} chromosome."
            echo `date`
            time ${BIN}/STAR \
                --genomeDir ${GENOME_INDEX_PREFIX}${chr_name}${GENOME_INDEX_SUFFIX} \
                --readFilesCommand zcat \
                --readFilesIn ${fastq_filepath} \
                --runThreadN ${nCPU} --genomeLoad NoSharedMemory \
                --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout
            echo `date`
            echo "finished."
            
            echo "start samtools for sorting BAM."
            echo `date`
            ${BIN}/samtools view -@ ${nCPU} -Shb Aligned.out.sam > Aligned.out.bam
            ${BIN}/samtools sort -@ ${nCPU} -o Aligned.out.refsort.bam Aligned.out.bam
            ${BIN}/samtools index -c Aligned.out.refsort.bam
            echo `date`
            echo "fnished."
            
            rm Aligned.out.sam
            cd -
        done
    done
fi





## EAGLE processes (multi-CPUs)
if [ ${pEAGLE} -eq 1 ]; then
    
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq` ; do
        echo ${seqlib}
        chr_names=(A B D)
        for chr_name in ${chr_names[@]}; do
            echo ${chr_name}
            cd ${seqlib}__on__chr${chr_name}
            for vcf_file in `ls ${GENOME_DIR}/${chr_name}.vs.*.gtf.vcf`; do
                vcf_basename=`basename ${vcf_file} .gtf.vcf`
                echo "start eagle"
                echo `date`
                time ${BIN}/eagle -t ${nCPU} -a Aligned.out.refsort.bam \
                     -r ${GENOME_DIR}/chr${chr_name}.fa \
                     -v ${vcf_file} \
                     --splice --rc 1> eagle.out.${vcf_basename}.txt 2> eagle.out.${vcf_basename}.readinfo.txt
                echo `date`
                echo "finished eagle"
            done
            cd -
        done
    done
    
fi





## EAGLE-RC processes (single-CPUs)
if [ ${pEAGLERC} -eq 1 ]; then
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq ` ; do
        echo ${seqlib}
        chr_names=(A B D)
        for chr_name in ${chr_names[@]}; do
            echo ${chr_name}
            cd ${seqlib}__on__chr${chr_name}
            for vcf_file in `ls ${GENOME_DIR}/${chr_name}.vs.*.gtf.vcf`; do
                vcf_basename=`basename ${vcf_file} .gtf.vcf`
                echo "start eaglerc"
                echo `date`
                time ${BIN}/eagle-rc --listonly -a Aligned.out.refsort.bam \
                                -o eaglerc.${vcf_basename} \
                                -v eagle.out.${vcf_basename}.txt eagle.out.${vcf_basename}.readinfo.txt > eaglerc.out.${vcf_basename}.list
                echo `date`
                echo "finished eaglerc"
            done
            cd ../
        done
    done
fi






# EAGLE-RC & featureCounts
cd ${PROJECT_DIR}
pyenv local eagle

if [ ${pCountHomeolog} -eq 1 ]; then
    mkdir -p ${EAGLE_DIR}
    cd ${EAGLE_DIR}
    
    for seqlib in `ls  ${BAM_DIR} | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq`; do
        bam_on__A=${BAM_DIR}/${seqlib}__on__chrA
        bam_on__B=${BAM_DIR}/${seqlib}__on__chrB
        bam_on__D=${BAM_DIR}/${seqlib}__on__chrD

        mkdir -p ${seqlib}
        cd ${seqlib}
        python ${EAGLE_SCRIPT_DIR}/ref3_consensus.py -d -u -o ${seqlib}.ref        \
           -A ${bam_on__A}/eaglerc.out.A.vs.B.list ${bam_on__A}/eaglerc.out.A.vs.D.list \
           -B ${bam_on__B}/eaglerc.out.B.vs.A.list ${bam_on__B}/eaglerc.out.B.vs.D.list \
           -D ${bam_on__D}/eaglerc.out.D.vs.A.list ${bam_on__D}/eaglerc.out.D.vs.B.list

        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam -o ${seqlib}.chrA ${seqlib}.ref.chrA.list
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam -o ${seqlib}.chrB ${seqlib}.ref.chrB.list
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam -o ${seqlib}.chrD ${seqlib}.ref.chrD.list

        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.txt ${seqlib}.chrA.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.txt ${seqlib}.chrB.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.txt ${seqlib}.chrD.ref.bam
        cd -
    done


    cd ${EAGLE_DIR}
    mkdir ${COUNTS_DIR}
    for d in `ls`; do
        cp ${d}/*.chrA.counts.txt ${COUNTS_DIR}/
        cp ${d}/*.chrB.counts.txt ${COUNTS_DIR}/
        cp ${d}/*.chrD.counts.txt ${COUNTS_DIR}/
    done

    # cd ${COUNTS_DIR}
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -skip 1 -a -i 0 -c 6 *.chrA.counts.txt > eagle.chrA.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -skip 1 -a -i 0 -c 6 *.chrB.counts.txt > eagle.chrB.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -skip 1 -a -i 0 -c 6 *.chrD.counts.txt > eagle.chrD.tsv

    # python ${EAGLE_SCRIPT_DIR}/tablize.py -a ${GENOME_DIR}/homeolog.A.list ./eagle.chrA.tsv | sort -k1 > eagle.chrA.homeolog.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -a ${GENOME_DIR}/homeolog.B.list ./eagle.chrB.tsv | cut -f 2,3- | sort -k1 > eagle.chrB.homeolog.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -a ${GENOME_DIR}/homeolog.D.list ./eagle.chrD.tsv | cut -f 2,3- | sort -k1 > eagle.chrD.homeolog.tsv

    # perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.chrA.homeolog.tsv > temp.txt
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -add temp.txt > eagle.chrA.homeolog.genelevel.tsv
    # perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.chrB.homeolog.tsv > temp.txt
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -add temp.txt > eagle.chrB.homeolog.genelevel.tsv
    # perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.chrD.homeolog.tsv > temp.txt
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -add temp.txt > eagle.chrD.homeolog.genelevel.tsv
fi






# Subgenome unique mapped reads
if [ ${pCountUnique} -eq 1 ]; then
    cd ${EAGLE_DIR}

    for seqlib in `ls  ${BAM_DIR} | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq`; do
        bam_on__A=${BAM_DIR}/${seqlib}__on__chrA
        bam_on__B=${BAM_DIR}/${seqlib}__on__chrB
        bam_on__D=${BAM_DIR}/${seqlib}__on__chrD

        cd ${seqlib}
        echo "" > dummy.txt
        echo ${seqlib}
        echo 'start eaglerc for qunily reads.'
        echo `date`
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam \
                                  -u ${bam_on__B}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrA.uniq dummy.txt
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrB.uniq dummy.txt
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__B}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrD.uniq dummy.txt
        echo `date`
        echo 'finished.'

        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.uniq.counts.txt ${seqlib}.chrA.uniq.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.uniq.counts.txt ${seqlib}.chrB.uniq.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.uniq.counts.txt ${seqlib}.chrD.uniq.ref.bam
        cd -
    done


    cd ${EAGLE_DIR}
    mkdir ${UNIQ_COUNTS_DIR}
    for d in `ls`; do
        cp ${d}/*.uniq.counts.txt ${UNIQ_COUNTS_DIR}/
    done
    cd ${UNIQ_COUNTS_DIR}

    # python ${EAGLE_SCRIPT_DIR}/tablize.py -skip 1 -a -i 0 -c 6 *.chrA.uniq.counts.txt > eagle.chrA.uniq.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -skip 1 -a -i 0 -c 6 *.chrB.uniq.counts.txt > eagle.chrB.uniq.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -skip 1 -a -i 0 -c 6 *.chrD.uniq.counts.txt > eagle.chrD.uniq.tsv

    # Subgenome unique mapped reads in subgenome unique genes (i.e. non-homeologs)
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -a ${GENOME_DIR}/chrA.uniq.list ./eagle.chrA.uniq.tsv > eagle.chrA.uniq.genelevel.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -a ${GENOME_DIR}/chrB.uniq.list ./eagle.chrB.uniq.tsv > eagle.chrB.uniq.genelevel.tsv
    # python ${EAGLE_SCRIPT_DIR}/tablize.py -a ${GENOME_DIR}/chrD.uniq.list ./eagle.chrD.uniq.tsv > eagle.chrD.uniq.genelevel.tsv
fi



echo `date`
echo "----- finish -----"
