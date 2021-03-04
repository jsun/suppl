# Indexing IWGSC v1.1 genome sequence


## Data preparation


Download sequences and annotations from IWGSC website.

```
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/iwgsc_refseqv1.0_all_chromosomes.zip
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/iwgsc_refseqv1.0_all_chromosomes.zip.md5.txt
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/iwgsc_refseqv1.1_genes_2017July06.zip
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/iwgsc_refseqv1.1_genes_2017July06.zip.md5.txt
```




## HISAT index

```
cd /home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
mkdir -p index

cd shell
qsub make_index_hisat.sh

cd -


# then extend 3'-dicrection 1kbp
cd iwgsc_refseqv1.1_genes_2017July06
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 1000 > IWGSC_v1.1_HC_20170706.ext.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3  500 > IWGSC_v1.1_HC_20170706.ext_0.5k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 1000 > IWGSC_v1.1_HC_20170706.ext_1.0k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 1500 > IWGSC_v1.1_HC_20170706.ext_1.5k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 2000 > IWGSC_v1.1_HC_20170706.ext_2.0k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 2500 > IWGSC_v1.1_HC_20170706.ext_2.5k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 3000 > IWGSC_v1.1_HC_20170706.ext_3.0k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 3500 > IWGSC_v1.1_HC_20170706.ext_3.5k.gff3
python ../scripts/extend_gff.py IWGSC_v1.1_HC_20170706.gff3 4000 > IWGSC_v1.1_HC_20170706.ext_4.0k.gff3


grep '^chr.A' IWGSC_v1.1_HC_20170706.ext_1.0k.gff3 > IWGSC_v1.1_HC_20170706.ext_1.0k.chrA.gff3
grep '^chr.B' IWGSC_v1.1_HC_20170706.ext_1.0k.gff3 > IWGSC_v1.1_HC_20170706.ext_1.0k.chrB.gff3
grep '^chr.D' IWGSC_v1.1_HC_20170706.ext_1.0k.gff3 > IWGSC_v1.1_HC_20170706.ext_1.0k.chrD.gff3

```


## EAGLE index (STAR index for each chromosome)

```
cd shell
qsub make_index_star.sh
```




## Homeolog identification


The following steps is used for identificating homeologs
by using reciprocal best hit approach with LAST alignemnt results.



```sh
# split whole genome into A, B, and D subgenomes

cd iwgsc_refseqv1.1_genes_2017July06

gffread -T -o IWGSC_v1.1_HC_20170706.gtf IWGSC_v1.1_HC_20170706.gff

grep '^chr.A' IWGSC_v1.1_HC_20170706.gtf > IWGSC_v1.1_HC_20170706.chrA.gtf
grep '^chr.B' IWGSC_v1.1_HC_20170706.gtf > IWGSC_v1.1_HC_20170706.chrB.gtf
grep '^chr.D' IWGSC_v1.1_HC_20170706.gtf > IWGSC_v1.1_HC_20170706.chrD.gtf

grep $'mRNA\t' IWGSC_v1.1_HC_20170706.gff3 | grep $'chr.A\t' | perl -ne 'chomp; m/ID=(.*?);/; print "$1\n";' > IWGSC_v1.1_HC_20170706.chrA.cds.list
grep $'mRNA\t' IWGSC_v1.1_HC_20170706.gff3 | grep $'chr.B\t' | perl -ne 'chomp; m/ID=(.*?);/; print "$1\n";' > IWGSC_v1.1_HC_20170706.chrB.cds.list
grep $'mRNA\t' IWGSC_v1.1_HC_20170706.gff3 | grep $'chr.D\t' | perl -ne 'chomp; m/ID=(.*?);/; print "$1\n";' > IWGSC_v1.1_HC_20170706.chrD.cds.list


cd ..

gffread -g iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -w chrA.cds.fa ../iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrA.gtf
gffread -g iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -w chrB.cds.fa ../iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrB.gtf
gffread -g iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -w chrD.cds.fa ../iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrD.gtf

python scripts/calc_cds_length.py > cds.length.tsv


cd ..

# make sub-genome sequences (NOT CDS sequences !!)
python scripts/make_subgenome_fasta.py iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta \
                                       chrA.fa chrB.fa chrD.fa chrUnknown.fa


# make LAST index and perform reciprocal best hit

mkdir index
qsub shell/make_index_last.sh


# use EAGEL script, set the path
script=/home/jqsun/local/src/eagle/scripts

# SNP detection and homeolog identification
gtf=./iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gtf
python ${script}/homeolog_genotypes.py -o A.vs.B -f exon -g ${gtf} A.B.maf B.A.maf   # coordinates based on A
python ${script}/homeolog_genotypes.py -o B.vs.A -f exon -g ${gtf} B.A.maf A.B.maf   # coordinates based on B
python ${script}/homeolog_genotypes.py -o B.vs.D -f exon -g ${gtf} B.D.maf D.B.maf   # coordinates based on B
python ${script}/homeolog_genotypes.py -o D.vs.B -f exon -g ${gtf} D.B.maf B.D.maf   # coordinates based on D
python ${script}/homeolog_genotypes.py -o A.vs.D -f exon -g ${gtf} A.D.maf D.A.maf   # coordinates based on A
python ${script}/homeolog_genotypes.py -o D.vs.A -f exon -g ${gtf} D.A.maf A.D.maf   # coordinates based on D

perl ${script}/triple_homeolog.pl A.vs.B.reciprocal_best B.vs.D.reciprocal_best A.vs.D.reciprocal_best > homeolog.ABD.list

cat A.vs.B.reciprocal_best A.vs.D.reciprocal_best | cut -f1 | sort | uniq > A.vs.all.list
cat B.vs.A.reciprocal_best B.vs.D.reciprocal_best | cut -f1 | sort | uniq > B.vs.all.list
cat D.vs.A.reciprocal_best D.vs.B.reciprocal_best | cut -f1 | sort | uniq > D.vs.all.list

python ${script}/tablize.py -v0 A.vs.all.list iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrA.cds.list > chrA.uniq.list
python ${script}/tablize.py -v0 B.vs.all.list iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrB.cds.list > chrB.uniq.list
python ${script}/tablize.py -v0 D.vs.all.list iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.chrD.cds.list > chrD.uniq.list

# hash table for changing B, D genome id to A genome id
cut -f 1 homeolog.ABD.list > homeolog.A.list
awk '{print $2"\t"$1;}' homeolog.ABD.list > homeolog.B.list
awk '{print $3"\t"$1;}' homeolog.ABD.list > homeolog.D.list
```





