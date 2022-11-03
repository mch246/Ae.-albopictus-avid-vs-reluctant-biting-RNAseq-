# KEGG analysis
# The pathway enrichment analysis for this I modified from "https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html"

# BiocManager::install("KEGGREST")

# Load libraries
install.packages("KEGGREST")
library(KEGGREST)

# Set some variables
setwd("/Users/mch284/Documents/Armbruster Lab/Auto_biting")
genes_with_pvalues <- "./BvN_LFCshrink_padj.csv" # with at least columns "gene", "padj"

keggGeneID_output <- "./BVN_keggID.txt"

# Read in data
gene_list <- read.csv(genes_with_pvalues, header = T)
View(gene_list)

# Make a list of all the ncbi to albopictus gene ids
convs <- keggConv("ncbi-geneid", "aalb")
head(convs)
#convs2 <- keggConv("aalb", "uniprot")
#convs3 <- keggConv("aalb", "ncbi-proteinid")

# Convert gene ids from list to something relevant to our work. Note that all our gene ids begin with "LOC". For NCBI LOC1234 is equivalent to GeneID = 1234. The "LOC"+GeneID is when orthologs have not yet been determined and a published symbol is not available. 
gene_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", gene_list$gene)
head(gene_list)

# Find the matching ncbi id in the conversion list and take the name associated (the kegg id) and assign that in the kegg_id column of the gene_list dataframe
gene_list$kegg_id = names(convs)[match(gene_list$ncbi_geneid, as.character(convs))]
#gene_list$uniprot_id = names(convs2)[match(gene_list$kegg_id, as.character(convs2))]
#gene_list$protein_id = names(convs3)[match(gene_list$kegg_id, as.character(convs3))] # getting the protein id does nothing for us downstream as these are not searchable in the fasta file. I think this is likely due to protein "versions" (basically all the protein ids in the fasta end with a ".#" following so I think they are alternative proteins for the genes)
head(gene_list)
# If you want to write out the KEGG ID list and do this online
write.table(gene_list$kegg_id, 
            file=keggGeneID_output, col.names = F, row.names = F, quote = F)
###
# You can input this KEGG list at "https://www.kegg.jp/kegg/tool/map_pathway1.html" to find enriched pathways
###


########
######
####
#Trying to automate this process with KEGG pathway enrichment
#all_genes_list <- read.csv(file="../misc/DESeq_results_pharatelarvae.csv")
all_genes_list <- read.csv(file="resOrdered_BVN.csv")
head(all_genes_list)
#all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$geneID)
all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$X)#note: need to rename x to gene
head(all_genes_list)
all_genes_list$kegg_id = names(convs)[match(all_genes_list$ncbi_geneid, as.character(convs))]


# Get the pathways list from KEGG
pathways.list <- keggList("pathway", "aalb")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
head(pathway.codes)
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

#geneList <- all_genes_list$X11d_padj #CHANGE THIS LINE
geneList <- all_genes_list$padj #CHANGE THIS LINE
head(geneList)
#geneLFClist <- all_genes_list$X11d_Log2FoldChange #CHANGE THIS LINE
geneLFClist <- all_genes_list$log2FoldChange 

#names(geneList) <- sub("aalb:","", all_genes_list$kegg_id) # get rid of the beginning "aalb:" since the gene list we brought from kegg doesn't have this
names(geneList) <- sub("aalb:","", all_genes_list$kegg_id) # get rid of the beginning "aalb:" since the gene list we brought from kegg doesn't have this
head(geneList)
View(all_genes_list)
View(geneList)
names(geneLFClist) <- sub("aalb:","", all_genes_list$kegg_id)
head(geneLFClist)

#genes.by.pathway_40d <- genes.by.pathway[-c(99, 120)]

pathway_pval <- data.frame()

for (pathway in 1:length(genes.by.pathway)){
  pathway.genes <- genes.by.pathway[[pathway]]
  if (!is.na(pathway.genes)){
    list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
    list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
    scores.in.pathway <- geneList[list.genes.in.pathway]
    scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
    if (length(scores.in.pathway) > 0){
      p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
    } else{
      p.value <- NA
    }
    new_row <- c(names(genes.by.pathway[pathway]), p.value, length(list.genes.in.pathway), sum(abs(geneLFClist[list.genes.in.pathway])>0.58&geneList[list.genes.in.pathway]<0.05), sum(geneLFClist[list.genes.in.pathway]>0.58&geneList[list.genes.in.pathway]<0.05), sum(geneLFClist[list.genes.in.pathway]< -0.58&geneList[list.genes.in.pathway]<0.05))
    pathway_pval <- rbind(pathway_pval, new_row)
  }
}

colnames(pathway_pval) <- c("pathwayCode", "pval", "annotated", "DEG", "up", "down")
pathway_pval <- pathway_pval[complete.cases(pathway_pval),]

pathway_pval$pathwayName = pathways.list[match(pathway_pval$pathwayCode, sub("path:","", names(pathways.list)))]

head(pathway_pval)
pathway_pval$pval <- as.numeric(pathway_pval$pval)

pathway_pval <- pathway_pval[order(pathway_pval$pval),]
head(pathway_pval)
View(pathway_pval)
# Write out a csv with these data
write.csv(pathway_pval, 
          file="./BvN_keggPathwayEnrichment_allgenes.csv", row.names = F)
