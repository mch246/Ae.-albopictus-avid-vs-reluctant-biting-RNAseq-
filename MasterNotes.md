# Workflow/pipeline details for the albopictus biting RNA-seq analysis
Note that this pipeline was modified based on previous work by Sarah Marzec analyzing RNAseq data comparing gene expression differences in biting vs non biting samples of two Culex subspecies.  See here: https://github.com/srmarzec/Culex_Biting_RNAseq  

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found in the following google sheets: ([link](https://docs.google.com/spreadsheets/d/1_RixzDGNsUlvhOMVCTXViuR_Rim8dNrS5joLMOwbsXo/edit#gid=0))

### Data Accession
Data was generated at Georgetown University using Ae. albopictus specimens collected from Manassas, Virginia in 2018 and maintained in the lab for 9 subsequent generations preceding sample collection

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/Albo_biting_trim_mRNA.SBATCH))

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/FastQC.SBATCH))

Preliminary trimming and fastqc showed a poor "per sequence base content" for the first ~15 bases.  Most samples also had warnings or fail for "per sequence GC content" and "sequence duplication levels. We used HEADCROP flag to remove the first 15 bases. All other flags (TRAILING, SLIDINGWINDOW, and MINLEN) were general/default for basic quality of bases and did not result in much difference of trimming.

### Mapping
We will be mapping with STAR (v2.7.1a).

Mapping was done using the following Aedes albopictus reference genome, GenBank assembly accession: GCA_018104305.1 (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/104/305/GCA_018104305.1_AalbF3/GCA_018104305.1_AalbF3_genomic.fna.gz) and the corresponding annotation file (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/104/305/GCA_018104305.1_AalbF3/GCA_018104305.1_AalbF3_genomic.gbff.gz)

STAR (v2.7.1a) was used for indexing the genome ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/STAR_index.SBATCH))

Reads were mapped in a two pass method. The first pass followed typical method with splice junctions from annotations ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/STAR_mapping.SBATCH)). The second pass is similar except that it additionally uses the output splice junctions info from the first pass (these would be novel splice junctions) to facilitate mapping ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/STAR_two_pass.SBATCH)).

Output sam files were converted to bam and then indexed ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/sam2bam.SBATCH))

### HTSeq Gene Counts
HTSeq (v0.13.5) was used to counts reads mapped to genes for downstream analyses ([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/htseq_count.SBATCH))

# Downstream analysis of gene counts using R (v4.1.3)

### DESeq2 (v1.34.0) package in R was used to identify differentially expressed genes
([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/albo_biting_DEseq.R))

### KEGGREST (v1.34.0) package in R was used to identify significantly enriched pathways (no sig. pathways identified)
([script](https://github.com/mch246/Ae.-albopictus-avid-vs-reluctant-biting-RNAseq-/blob/main/KEGG_Enrichment_BVN.R))
