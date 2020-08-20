# Timecourse RNA-Seq data analysis of C. insueta


## QC

Qaulity controlling is performed by Trimmomatic in `qc` directory.

```
cd qc/
bash qc.sh
```



## Assembly

A-genome is assembled with RNA-Seq reads of C. amara (leaf, 0-96 hr),
and R-genome is assembled with RNA-Seq reads of C. rivularis (leaf, 0-96 hr).
The genome of C. hirsuta (H-genome) is used as reference-guide.
Sequences and annotations of H-genome is stored in `data/genome/chirsuta` directory.

```
cd assembly/
bash refassemble_rnaseq.sh camara
bash refassemble_rnaseq.sh crivularis
```

The assembled A-genome and R-genome is saved in
`data/genome_leaf/camara` and `data/genome_leaf/crivularis` directories, respectively.


The mapping qualities are calulated for evaluating the assembly qualities.


```
bash calc_stats.sh camara
bash calc_stats.sh crivularis
```




## Mapping

Mapping of RNA-Seq reads onto A-genome and R-genome is performed
with STAR, and read classification is performed with HomeoRoq.


```
cd rnaseq_processing/
bash proc_main.sh
```





## Expression analysis

Expression analysis is mainly performed by R script.
However, it is required to create gene ontology annotations before
expression analysis.

To prepare the gene ontology annotations for Cardamine speceis gene ID,
use the following command to download gene ontology annotations for
A. thaliana, and adjust the annotations to Cardamine species..


```
cd data_analysis/
bash create_go_annotations.sh
```

Move mapping and counting results form `rnaseq_processing` to `data_analysis` directory.

```
mv rnaseq_processing/counts/results/counts data_analysis/data/
```

Then, calculate exons length for each gene to calculate FPKM.


```
cd data_analysis/lib/src
python ${SCRIPTDIR}/calc_exons_length.py --gffa ${GENOMELEAFDIR}/camara/genome.modified.gff \
                                         --gffr ${GENOMELEAFDIR}/crivularis/genome.modified.gff \
                                         --output gene_length.tsv

python ${SCRIPTDIR}/create_gene_chr_list.py --gff ${GENOMELEAFDIR}/camara/genome.modified.gff > gene_chr.txt

```

Then, use the `main_analysis.R` to perform gene expression analysis.




