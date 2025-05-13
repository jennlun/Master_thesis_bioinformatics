## Set Working directory
setwd("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/3_GO_analysis")

## Load packages
library(goseq)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggh4x)

##--Upload data files--
#"Full" GO file
swissprot_go <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/1_file_generation/output_file_generation/swissprot_go.txt", 
                           header=T, stringsAsFactors=FALSE, sep="\t")

#"Long" GO file
ccar2swissprotGO_long <- read.delim("./GO_source_files/ccar2swissprotGO_long.txt", header=T, stringsAsFactors=FALSE, sep="\t")

#Make a named list for use with goseq package
ccar_gene2cat <- split(ccar2swissprotGO_long$GO_id, ccar2swissprotGO_long$gene_id) 

#Upload reformatted GO information file 
go_data = read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/1_file_generation/output_file_generation/go_data.txt")

#Read in the table with gene ids and mean, median gene length, and length of the longest isoform, and merged (exons) 
lengths_raw <- read.delim("./GO_source_files/carcar_annotation_v5.genelength", header=TRUE, stringsAsFactors=FALSE, sep="\t") 

##Upload DEG files from DESeq
#Upload Liver dataframes
deg_Liver_AvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_AvsN.txt")
deg_Liver_3hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_3hRvsN.txt")
deg_Liver_24hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_24hRvsN.txt")
deg_Liver_3hRvsA <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_3hRvsA.txt")
deg_Liver_24hRvsA <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_24hRvsA.txt")
deg_Liver_24hRvs3hR <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_24hRvs3hR.txt")

#Upload Gill dataframes
deg_Gill_AvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_AvsN.txt")
deg_Gill_3hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_3hRvsN.txt")
deg_Gill_24hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_24hRvsN.txt")
deg_Gill_3hRvsA <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_3hRvsA.txt")
deg_Gill_24hRvsA <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_24hRvsA.txt")
deg_Gill_24hRvs3hR <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_24hRvs3hR.txt")



#-----Step 1: Make Subsets Relating to Anoxia and Reoxygenation-----

#--Make subsets for Liver (use logFC 0 because it is already taken into account in the modelling)--

#Upregulated genes
#Store genes that are significantly upregulated (padj <0.05) and have a postitive logfold change, as they are then upregulated
#which() command ensures that both criteria are met
up_Liver_deg_AvsN=deg_Liver_AvsN$gene_id[which(deg_Liver_AvsN$padj < 0.05 & deg_Liver_AvsN$log2FoldChange > 0)]
up_Liver_deg_3hRvsN=deg_Liver_3hRvsN$gene_id[which(deg_Liver_3hRvsN$padj < 0.05 & deg_Liver_3hRvsN$log2FoldChange > 0)]
up_Liver_deg_24hRvsN=deg_Liver_24hRvsN$gene_id[which(deg_Liver_24hRvsN$padj < 0.05 & deg_Liver_24hRvsN$log2FoldChange > 0)]
up_Liver_deg_3hRvsA=deg_Liver_3hRvsA$gene_id[which(deg_Liver_3hRvsA$padj < 0.05 & deg_Liver_3hRvsA$log2FoldChange > 0)]
up_Liver_deg_24hRvsA=deg_Liver_24hRvsA$gene_id[which(deg_Liver_24hRvsA$padj < 0.05 & deg_Liver_24hRvsA$log2FoldChange > 0)]
up_Liver_deg_24hRvs3hR=deg_Liver_24hRvs3hR$gene_id[which(deg_Liver_24hRvs3hR$padj < 0.05 & deg_Liver_24hRvs3hR$log2FoldChange > 0)]

#Do the same for upregulated genes
#Store genes with negative logfold change, as they are downregulated
down_Liver_deg_AvsN=deg_Liver_AvsN$gene_id[which(deg_Liver_AvsN$padj < 0.05 & deg_Liver_AvsN$log2FoldChange < 0)]
down_Liver_deg_3hRvsN=deg_Liver_3hRvsN$gene_id[which(deg_Liver_3hRvsN$padj < 0.05 & deg_Liver_3hRvsN$log2FoldChange < 0)]
down_Liver_deg_24hRvsN=deg_Liver_24hRvsN$gene_id[which(deg_Liver_24hRvsN$padj < 0.05 & deg_Liver_24hRvsN$log2FoldChange < 0)]
down_Liver_deg_3hRvsA=deg_Liver_3hRvsA$gene_id[which(deg_Liver_3hRvsA$padj < 0.05 & deg_Liver_3hRvsA$log2FoldChange > 0)]
down_Liver_deg_24hRvsA=deg_Liver_24hRvsA$gene_id[which(deg_Liver_24hRvsA$padj < 0.05 & deg_Liver_24hRvsA$log2FoldChange > 0)]
down_Liver_deg_24hRvs3hR=deg_Liver_24hRvs3hR$gene_id[which(deg_Liver_24hRvs3hR$padj < 0.05 & deg_Liver_24hRvs3hR$log2FoldChange > 0)]

#identifiy the gene names that are significantly upregulated only in one comparison 
#(e.g. up_Liver_deg_3hRvsN) and not in another (e.g. up_Liver_deg_AvsN)
up_Liver_deg_3hR_only = setdiff(up_Liver_deg_3hRvsN, up_Liver_deg_AvsN)
up_Liver_deg_24hR_only =setdiff(up_Liver_deg_24hRvsN, union(up_Liver_deg_AvsN, up_Liver_deg_3hRvsN))

#identifiy the gene names that are significantly downregulated only in one comparison and not in another
down_Liver_deg_3hR_only = setdiff(down_Liver_deg_3hRvsN, down_Liver_deg_AvsN)
down_Liver_deg_24hR_only =setdiff(down_Liver_deg_24hRvsN, union(down_Liver_deg_AvsN, down_Liver_deg_3hRvsN))

#Combine all unique genes that are up- or downregulated from all the liver subsets
all_Liver_degs =unique(c(up_Liver_deg_AvsN, up_Liver_deg_3hRvsN, up_Liver_deg_3hRvsA, 
                         up_Liver_deg_24hRvsN, up_Liver_deg_24hRvs3hR, up_Liver_deg_24hRvsA, 
                         down_Liver_deg_AvsN, down_Liver_deg_3hRvsN, down_Liver_deg_3hRvsA, 
                         down_Liver_deg_24hRvsN, down_Liver_deg_24hRvs3hR, down_Liver_deg_24hRvsA))

#Save full file for plotting
writeLines(all_Liver_degs, "./GO_output/all_Liver_degs.txt")

#Select the top 10 genes with the smallest adjusted p-values from the DEG dataframes
toptags_Liver_AvsN = deg_Liver_AvsN %>% slice_min(padj, n = 10)
toptags_Liver_3hRvsN = deg_Liver_3hRvsN %>% slice_min(padj, n = 10)
toptags_Liver_24hRvsN = deg_Liver_24hRvsN %>% slice_min(padj, n = 10)
toptags_Liver_3hRvsA = deg_Liver_3hRvsA %>% slice_min(padj, n = 10)
toptags_Liver_24hRvs3hR = deg_Liver_24hRvs3hR %>% slice_min(padj, n = 10)
toptags_Liver_24hRvsA = deg_Liver_24hRvsA %>% slice_min(padj, n = 10)

#Combine all the "toptags" into one object
toptags_Liver <- unique(c(toptags_Liver_AvsN$gene_id, toptags_Liver_3hRvsN$gene_id, 
                          toptags_Liver_24hRvsN$gene_id, toptags_Liver_3hRvsA$gene_id, 
                          toptags_Liver_24hRvs3hR$gene_id, toptags_Liver_24hRvsA$gene_id))

#Save toptags for plotting
writeLines(toptags_Liver, "./GO_output/toptags_Liver.txt")


#--Make subsets for Gill (use logFC 0 because it is already taken into account in the modelling)--

#Upregulated genes
#Store genes that are significantly upregulated (padj <0.05) and have a postitive logfold change, as they are then upregulated
#which() command ensures that both criteria are met
up_Gill_deg_AvsN=deg_Gill_AvsN$gene_id[which(deg_Gill_AvsN$padj < 0.05 & deg_Gill_AvsN$log2FoldChange > 0)]
up_Gill_deg_3hRvsN=deg_Gill_3hRvsN$gene_id[which(deg_Gill_3hRvsN$padj < 0.05 & deg_Gill_3hRvsN$log2FoldChange > 0)]
up_Gill_deg_24hRvsN=deg_Gill_24hRvsN$gene_id[which(deg_Gill_24hRvsN$padj < 0.05 & deg_Gill_24hRvsN$log2FoldChange > 0)]
up_Gill_deg_3hRvsA=deg_Gill_3hRvsA$gene_id[which(deg_Gill_3hRvsA$padj < 0.05 & deg_Gill_3hRvsA$log2FoldChange > 0)]
up_Gill_deg_24hRvsA=deg_Gill_24hRvsA$gene_id[which(deg_Gill_24hRvsA$padj < 0.05 & deg_Gill_24hRvsA$log2FoldChange > 0)]
up_Gill_deg_24hRvs3hR=deg_Gill_24hRvs3hR$gene_id[which(deg_Gill_24hRvs3hR$padj < 0.05 & deg_Gill_24hRvs3hR$log2FoldChange > 0)]

#Do the same for upregulated genes
#Store genes with negative logfold change, as they are downregulated
down_Gill_deg_AvsN=deg_Gill_AvsN$gene_id[which(deg_Gill_AvsN$padj < 0.05 & deg_Gill_AvsN$log2FoldChange < 0)]
down_Gill_deg_3hRvsN=deg_Gill_3hRvsN$gene_id[which(deg_Gill_3hRvsN$padj < 0.05 & deg_Gill_3hRvsN$log2FoldChange < 0)]
down_Gill_deg_24hRvsN=deg_Gill_24hRvsN$gene_id[which(deg_Gill_24hRvsN$padj < 0.05 & deg_Gill_24hRvsN$log2FoldChange < 0)]
down_Gill_deg_3hRvsA=deg_Gill_3hRvsA$gene_id[which(deg_Gill_3hRvsA$padj < 0.05 & deg_Gill_3hRvsA$log2FoldChange > 0)]
down_Gill_deg_24hRvsA=deg_Gill_24hRvsA$gene_id[which(deg_Gill_24hRvsA$padj < 0.05 & deg_Gill_24hRvsA$log2FoldChange > 0)]
down_Gill_deg_24hRvs3hR=deg_Gill_24hRvs3hR$gene_id[which(deg_Gill_24hRvs3hR$padj < 0.05 & deg_Gill_24hRvs3hR$log2FoldChange > 0)]

#identifiy the gene names that are significantly upregulated only in one comparison 
#(e.g. up_Gill_deg_3hRvsN) and not in another (e.g. up_Gill_deg_AvsN)
up_Gill_deg_3hR_only = setdiff(up_Gill_deg_3hRvsN, up_Gill_deg_AvsN)
up_Gill_deg_24hR_only =setdiff(up_Gill_deg_24hRvsN, union(up_Gill_deg_AvsN, up_Gill_deg_3hRvsN))

#identifiy the gene names that are significantly downregulated only in one comparison and not in another
down_Gill_deg_3hR_only = setdiff(down_Gill_deg_3hRvsN, down_Gill_deg_AvsN)
down_Gill_deg_24hR_only =setdiff(down_Gill_deg_24hRvsN, union(down_Gill_deg_AvsN, down_Gill_deg_3hRvsN))

#Combine all unique genes that are up- or downregulated from all the liver subsets
all_Gill_degs =unique(c(up_Gill_deg_AvsN, up_Gill_deg_3hRvsN, up_Gill_deg_3hRvsA, up_Gill_deg_24hRvsN, 
                        up_Gill_deg_24hRvs3hR, up_Gill_deg_24hRvsA, down_Gill_deg_AvsN, down_Gill_deg_3hRvsN, 
                        down_Gill_deg_3hRvsA, down_Gill_deg_24hRvsN, down_Gill_deg_24hRvs3hR, down_Gill_deg_24hRvsA))

#Save full file for plotting
writeLines(all_Gill_degs, "./GO_output/all_Gill_degs.txt" )


#Select the top 10 genes with the smallest adjusted p-values from the DEG dataframes
toptags_Gill_AvsN = deg_Gill_AvsN %>% slice_min(padj, n = 10)
toptags_Gill_3hRvsN = deg_Gill_3hRvsN %>% slice_min(padj, n = 10)
toptags_Gill_24hRvsN = deg_Gill_24hRvsN %>% slice_min(padj, n = 10)
toptags_Gill_3hRvsA = deg_Gill_3hRvsA %>% slice_min(padj, n = 10)
toptags_Gill_24hRvs3hR = deg_Gill_24hRvs3hR %>% slice_min(padj, n = 10)
toptags_Gill_24hRvsA = deg_Gill_24hRvsA %>% slice_min(padj, n = 10)

#Combine all the "toptags" into one dataframe
toptags_Gill <- unique(c(toptags_Gill_AvsN$gene_id, toptags_Gill_3hRvsN$gene_id, 
                         toptags_Gill_24hRvsN$gene_id, toptags_Gill_3hRvsA$gene_id, 
                         toptags_Gill_24hRvs3hR$gene_id, toptags_Gill_24hRvsA$gene_id))

#Save toptags for plotting
writeLines(toptags_Gill, "./GO_output/toptags_Gill.txt")


#-----Step 2: Prepare files for GO analysis-----

#Create a vector containing all the unique gene IDs in the "full" swissprot GO file
go_genes <- unique(swissprot_go$gene_id)

#Only keep lengths for the gene ids for which there was a GO term and that are included in filtered count data
lengths_raw_go <- lengths_raw[lengths_raw$gene %in% go_genes, ] 

#Filter lengths_raw_go based on previously determined genes of interest
lengths_raw_go_filt <- lengths_raw_go[lengths_raw_go$gene %in% deg_Liver_AvsN$gene_id, ] 

# use one of the results tables from DESeq analysis to get a list of all the gene ids, 
#but only use the gene ids that has GO terms
all.genes <- intersect(go_genes,deg_Liver_AvsN$gene_id) 

#Save a file with all genes ("background set")
write.table(all.genes, "./GO_output/all_genes.txt", row.names = FALSE, col.names = FALSE, quote = F)

#Make a numeric vector of the column containing the median length (the median across the isoforms)
lengths <- as.numeric(lengths_raw_go_filt[,2]) 

#Convert the numeric vector to a named numeric vector, i.e. each value in the vector is named by the gene id
names(lengths) <- as.character(lengths_raw_go_filt$gene) 



#-----Step 3: Run GO analyses on the genes differentially regulated in AvsN-----

#--3.1: upregulated AvsN genes Liver--

#Create a factor that specifies if the genes in all.genes are present in th AvsN liver upreulated genes
#1 = present, 0 = not present
up_AvsN_genes_Liver <- factor(as.integer(all.genes %in% up_Liver_deg_AvsN))

#Make it so each element of up_AvsN_genes_Liver will have a name corresponding to the corresponding element in all.genes
names(up_AvsN_genes_Liver) = all.genes

#show a table of how many genes were assigned 1 or 0
table(up_AvsN_genes_Liver)

#show the first lines/elements of up_AvsN_genes_Liver
head(up_AvsN_genes_Liver)

#Use the nullp function to generate a probability weight function (pwf) that accounts for potential biases in the gene data
#Generates a plot
up_AvsN_pwf_Liver=nullp(up_AvsN_genes_Liver, bias.data = lengths)

#Run GO analysis on the upregulated Liver samples
up_AvsN_GO_Liver=goseq(up_AvsN_pwf_Liver, gene2cat = ccar_gene2cat)

#Make a dataframe containing only the gene IDs and go IDs that are statistically significally upregulated in AvsN in the liver
up_AvsN_go_genes_Liver=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% up_Liver_deg_AvsN, c(1,2)] 

#Group the data by GO category, i.e. column 1 is the category and column 2 has a comma-separated list of the DEGs with that category. 
#I.e. these are not all the genes that has this term, but only the ones that were also DEGs.
#Group the df based on the GO IDs, then combine all gene IDs with the same GO ID in one "frame" separated by commas
up_AvsN_go_genes_Liver = up_AvsN_go_genes_Liver %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Add the comma-separated lists of gene ids to the results table
#Combine the GO analysis results with the list of genes associated with each GO category 
#pairs statistical results with specific genes related to those GO categories, helping in the interpretation and understanding of the biological significance of the results
up_AvsN_GO_Liver = merge(up_AvsN_GO_Liver, up_AvsN_go_genes_Liver, by.x = "category", by.y = "GO_id", all = T)

#Add more information from the GO data file
up_AvsN_GO_Liver <- merge(up_AvsN_GO_Liver, go_data, by.x = "category", by.y = "id")

#each entry in up_AvsN_GO_Liver$over_represented_pvalue is adjusted to account for the rate of false discoveries, 
#and these adjusted values are stored in a new column called over_represented_FDR. 
#This makes the dataset more reliable for identifying truly significant GO categories, 
#reducing the likelihood of reporting false positives due to the sheer number of tests performed.
up_AvsN_GO_Liver$over_represented_FDR = p.adjust(up_AvsN_GO_Liver$over_represented_pvalue, method="BH")

#Same as above but with under-represented values
up_AvsN_GO_Liver$under_represented_FDR = p.adjust(up_AvsN_GO_Liver$under_represented_pvalue, method="BH")

#calculate the expected number of DEGs within each GO category and assign the result to a new column (expectedDEInCat)
#expectedDEInCat gives a baseline expectation against which observed counts of differentially expressed genes can be compared. 
#This helps identify categories that are statistically enriched or depleted for differentially expressed genes
up_AvsN_GO_Liver$expectedDEInCat = (up_AvsN_GO_Liver$numInCat/length(up_AvsN_genes_Liver))*length(intersect(up_Liver_deg_AvsN, all.genes))

#Calculate the fold enrichment of DEGs in each GO category and assigns the result to a new column (fold_enrichment)
#fold enrichment is used to quantify how much more (or less) frequent DEGs are in a GO category compared to what would be expected by chance.
#calculating fold enrichment helps identify GO categories that may be biologically relevant in terms of DEG
up_AvsN_GO_Liver$fold_enrichment = up_AvsN_GO_Liver$numDEInCat/up_AvsN_GO_Liver$expectedDEInCat

#Calculates the negative logarithm (base 10) of the false discovery rate (FDR) for each GO category and assigns the result to a new column (-log10(FDR))
#This column reflects the degree of statistical significance of each GO category
#Larger values imply greater statistical significance (i.e., lower FDRs), indicating potentially meaningful enrichment beyond random chance
#Values near zero suggest less significant categories
up_AvsN_GO_Liver$'-log10(FDR)' = -log10(up_AvsN_GO_Liver$over_represented_FDR)

#create a new column called description, which links the term name and the category identifier for each GO category
#results in each entry in the description column being a string formatted as "TermName (CategoryID)", 
#making it easier to understand both the general meaning and specific identifier of each GO category.
#Useful for plotting etc
up_AvsN_GO_Liver$description = paste(up_AvsN_GO_Liver$term," (",up_AvsN_GO_Liver$category,")", sep="")

#Creates a new column "group" containing "anoxia" in each row
#categorizes the terms
up_AvsN_GO_Liver$group = paste("Anoxia")

#Same as above, but with the column "direction" containing "up-regulated"
up_AvsN_GO_Liver$direction = paste("up-regulated")

#Filter the df to only include rows where the false discovery rate (FDR) for over-representation is less than 0.05
#Basically, only keep GO categories with statistically significant over-representation
up_AvsN_GO_Liver <- up_AvsN_GO_Liver[up_AvsN_GO_Liver$over_represented_FDR < 0.05, ]

#Order results by fold enrichment to make the plot look nice (it will have fold enrichment on x axis)
up_AvsN_GO_Liver <- up_AvsN_GO_Liver[order(up_AvsN_GO_Liver$fold_enrichment),] 



#--3.2: downregulated AvsN genes Liver--

#For detailed explanations of steps, see 3.1

#Specify whether genes are downregulated, and name them
down_AvsN_genes_Liver <- factor(as.integer(all.genes %in% down_Liver_deg_AvsN))
names(down_AvsN_genes_Liver) = all.genes
table(down_AvsN_genes_Liver)
head(down_AvsN_genes_Liver)

#Prepare data for analysis and run GO analysis
down_AvsN_pwf_Liver=nullp(down_AvsN_genes_Liver, bias.data = lengths)
down_AvsN_GO_Liver=goseq(down_AvsN_pwf_Liver, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
down_AvsN_go_genes_Liver=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% down_Liver_deg_AvsN, c(1,2)] 
down_AvsN_go_genes_Liver = down_AvsN_go_genes_Liver %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the down_AvsN_GO_Liver data frame, focusing on down-regulated genes under an "Anoxia" condition
down_AvsN_GO_Liver = merge(down_AvsN_GO_Liver, down_AvsN_go_genes_Liver, by.x = "category", by.y = "GO_id", all = T)
down_AvsN_GO_Liver <- merge(down_AvsN_GO_Liver, go_data, by.x = "category", by.y = "id")
down_AvsN_GO_Liver$over_represented_FDR = p.adjust(down_AvsN_GO_Liver$over_represented_pvalue, method="BH")
down_AvsN_GO_Liver$under_represented_FDR = p.adjust(down_AvsN_GO_Liver$under_represented_pvalue, method="BH")
down_AvsN_GO_Liver$expectedDEInCat = (down_AvsN_GO_Liver$numInCat/length(down_AvsN_genes_Liver))*length(intersect(down_Liver_deg_AvsN, all.genes))
down_AvsN_GO_Liver$fold_enrichment = down_AvsN_GO_Liver$numDEInCat/down_AvsN_GO_Liver$expectedDEInCat
down_AvsN_GO_Liver$'-log10(FDR)' = -log10(down_AvsN_GO_Liver$over_represented_FDR)
down_AvsN_GO_Liver$description = paste(down_AvsN_GO_Liver$term," (",down_AvsN_GO_Liver$category,")", sep="")
down_AvsN_GO_Liver$group = paste("Anoxia")
down_AvsN_GO_Liver$direction = paste("down-regulated")
down_AvsN_GO_Liver <- down_AvsN_GO_Liver[down_AvsN_GO_Liver$over_represented_FDR < 0.05, ]
down_AvsN_GO_Liver <- down_AvsN_GO_Liver[order(down_AvsN_GO_Liver$fold_enrichment),]


# combine the results for up- and down-regulated genes:
GO_AvsN_Liver = rbind(up_AvsN_GO_Liver, down_AvsN_GO_Liver) #combine the rows, since the columns are the same



#--3.3: upregulated AvsN genes Gill--

#For detailed explanations of steps, see 3.1

#Specify whether genes are upregulated, and name them
up_AvsN_genes_Gill <- factor(as.integer(all.genes %in% up_Gill_deg_AvsN))
names(up_AvsN_genes_Gill) = all.genes
table(up_AvsN_genes_Gill)
head(up_AvsN_genes_Gill)

#Prepare data for analysis and run GO analysis
up_AvsN_pwf_Gill=nullp(up_AvsN_genes_Gill, bias.data = lengths)
up_AvsN_GO_Gill=goseq(up_AvsN_pwf_Gill, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
up_AvsN_go_genes_Gill=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% up_Gill_deg_AvsN, c(1,2)] 
up_AvsN_go_genes_Gill = up_AvsN_go_genes_Gill %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 


up_AvsN_GO_Gill = merge(up_AvsN_GO_Gill, up_AvsN_go_genes_Gill, by.x = "category", by.y = "GO_id", all = T)
up_AvsN_GO_Gill <- merge(up_AvsN_GO_Gill, go_data, by.x = "category", by.y = "id")
up_AvsN_GO_Gill$over_represented_FDR = p.adjust(up_AvsN_GO_Gill$over_represented_pvalue, method="BH")
up_AvsN_GO_Gill$under_represented_FDR = p.adjust(up_AvsN_GO_Gill$under_represented_pvalue, method="BH")
up_AvsN_GO_Gill$expectedDEInCat = (up_AvsN_GO_Gill$numInCat/length(up_AvsN_genes_Gill))*length(intersect(up_Gill_deg_AvsN, all.genes))
up_AvsN_GO_Gill$fold_enrichment = up_AvsN_GO_Gill$numDEInCat/up_AvsN_GO_Gill$expectedDEInCat
up_AvsN_GO_Gill$'-log10(FDR)' = -log10(up_AvsN_GO_Gill$over_represented_FDR)
up_AvsN_GO_Gill$description = paste(up_AvsN_GO_Gill$term," (",up_AvsN_GO_Gill$category,")", sep="")
up_AvsN_GO_Gill$group = paste("Anoxia")
up_AvsN_GO_Gill$direction = paste("up-regulated")
up_AvsN_GO_Gill <- up_AvsN_GO_Gill[up_AvsN_GO_Gill$over_represented_FDR < 0.05, ]
up_AvsN_GO_Gill <- up_AvsN_GO_Gill[order(up_AvsN_GO_Gill$fold_enrichment),] 



#--3.4: downregulated AvsN genes Gill--

#For detailed explanations of steps, see 3.1

#Specify whether genes are downregulated, and name them
down_AvsN_genes_Gill <- factor(as.integer(all.genes %in% down_Gill_deg_AvsN))
names(down_AvsN_genes_Gill) = all.genes
table(down_AvsN_genes_Gill)
head(down_AvsN_genes_Gill)

#Prepare data for analysis and run GO analysis
down_AvsN_pwf_Gill=nullp(down_AvsN_genes_Gill, bias.data = lengths)
down_AvsN_GO_Gill=goseq(down_AvsN_pwf_Gill, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
down_AvsN_go_genes_Gill=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% down_Gill_deg_AvsN, c(1,2)] 
down_AvsN_go_genes_Gill = down_AvsN_go_genes_Gill %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the down_AvsN_GO_Gill data frame, focusing on down-regulated genes under an "Anoxia" condition
down_AvsN_GO_Gill = merge(down_AvsN_GO_Gill, down_AvsN_go_genes_Gill, by.x = "category", by.y = "GO_id", all = T)
down_AvsN_GO_Gill <- merge(down_AvsN_GO_Gill, go_data, by.x = "category", by.y = "id")
down_AvsN_GO_Gill$over_represented_FDR = p.adjust(down_AvsN_GO_Gill$over_represented_pvalue, method="BH")
down_AvsN_GO_Gill$under_represented_FDR = p.adjust(down_AvsN_GO_Gill$under_represented_pvalue, method="BH")
down_AvsN_GO_Gill$expectedDEInCat = (down_AvsN_GO_Gill$numInCat/length(down_AvsN_genes_Gill))*length(intersect(down_Gill_deg_AvsN, all.genes))
down_AvsN_GO_Gill$fold_enrichment = down_AvsN_GO_Gill$numDEInCat/down_AvsN_GO_Gill$expectedDEInCat
down_AvsN_GO_Gill$'-log10(FDR)' = -log10(down_AvsN_GO_Gill$over_represented_FDR)
down_AvsN_GO_Gill$description = paste(down_AvsN_GO_Gill$term," (",down_AvsN_GO_Gill$category,")", sep="")
down_AvsN_GO_Gill$group = paste("Anoxia")
down_AvsN_GO_Gill$direction = paste("down-regulated")
down_AvsN_GO_Gill <- down_AvsN_GO_Gill[down_AvsN_GO_Gill$over_represented_FDR < 0.05, ]
down_AvsN_GO_Gill <- down_AvsN_GO_Gill[order(down_AvsN_GO_Gill$fold_enrichment),]


# combine the results for up- and down-regulated genes:
GO_AvsN_Gill = rbind(up_AvsN_GO_Gill, down_AvsN_GO_Gill) #combine the rows, since the columns are the same



#-----Step 4: Run GO analyses on the genes differentially regulated in 24hr Reox-----

#Runs analysis on genes that are only upregulated when compared to 24hr reox, and not any other comparison
#only in 24hr vs N, not in 3hr vs N or A vs N


#--4.1: upregulated 24hr Reox genes Liver--

#For detailed explanations of steps, see 3.1

#Specify whether genes are upregulated, and name them
up_R24hr_genes_Liver <- factor(as.integer(all.genes %in% up_Liver_deg_24hR_only))
names(up_R24hr_genes_Liver) = all.genes
table(up_R24hr_genes_Liver)
head(up_R24hr_genes_Liver)

#Prepare data for analysis and run GO analysis
up_R24hr_pwf_Liver=nullp(up_R24hr_genes_Liver, bias.data = lengths)
up_R24hr_GO_Liver=goseq(up_R24hr_pwf_Liver, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
up_R24hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% up_Liver_deg_24hR_only, c(1,2)] 
up_R24hr_go_genes = up_R24hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the up_R24hr_GO_Liver data frame, focusing on up-regulated genes under an "24hr Reox" condition
up_R24hr_GO_Liver = merge(up_R24hr_GO_Liver, up_R24hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
up_R24hr_GO_Liver <- merge(up_R24hr_GO_Liver, go_data, by.x = "category", by.y = "id")
up_R24hr_GO_Liver$over_represented_FDR = p.adjust(up_R24hr_GO_Liver$over_represented_pvalue, method="BH")
up_R24hr_GO_Liver$under_represented_FDR = p.adjust(up_R24hr_GO_Liver$under_represented_pvalue, method="BH")
up_R24hr_GO_Liver$expectedDEInCat = (up_R24hr_GO_Liver$numInCat/length(up_R24hr_genes_Liver))*length(intersect(up_Liver_deg_24hR_only, all.genes))
up_R24hr_GO_Liver$fold_enrichment = up_R24hr_GO_Liver$numDEInCat/up_R24hr_GO_Liver$expectedDEInCat
up_R24hr_GO_Liver$'-log10(FDR)' = -log10(up_R24hr_GO_Liver$over_represented_FDR)
up_R24hr_GO_Liver$description = paste(up_R24hr_GO_Liver$term," (",up_R24hr_GO_Liver$category,")", sep="")
up_R24hr_GO_Liver$group = paste("24hr Reox")
up_R24hr_GO_Liver$direction = paste("up-regulated")
up_R24hr_GO_Liver <- up_R24hr_GO_Liver[up_R24hr_GO_Liver$over_represented_FDR < 0.05, ]
up_R24hr_GO_Liver <- up_R24hr_GO_Liver[order(up_R24hr_GO_Liver$fold_enrichment),] 



#--4.2: downregulated 24hr Reox genes Liver--

#For detailed explanations of steps, see 3.1

#Specify whether genes are downregulated, and name them
down_R24hr_genes_Liver <- factor(as.integer(all.genes %in% down_Liver_deg_24hR_only))
names(down_R24hr_genes_Liver) = all.genes
table(down_R24hr_genes_Liver)
head(down_R24hr_genes_Liver)

#Prepare data for analysis and run GO analysis
down_R24hr_pwf_Liver=nullp(down_R24hr_genes_Liver, bias.data = lengths)
down_R24hr_GO_Liver=goseq(down_R24hr_pwf_Liver, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
down_R24hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% down_Liver_deg_24hR_only, c(1,2)] 
down_R24hr_go_genes = down_R24hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the down_R24hr_GO_Liver data frame, focusing on down-regulated genes under an "24hr Reox" condition
down_R24hr_GO_Liver = merge(down_R24hr_GO_Liver, down_R24hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
down_R24hr_GO_Liver <- merge(down_R24hr_GO_Liver, go_data, by.x = "category", by.y = "id")
down_R24hr_GO_Liver$over_represented_FDR = p.adjust(down_R24hr_GO_Liver$over_represented_pvalue, method="BH")
down_R24hr_GO_Liver$under_represented_FDR = p.adjust(down_R24hr_GO_Liver$under_represented_pvalue, method="BH")
down_R24hr_GO_Liver$expectedDEInCat = (down_R24hr_GO_Liver$numInCat/length(down_R24hr_genes_Liver))*length(intersect(down_Liver_deg_24hR_only, all.genes))
down_R24hr_GO_Liver$fold_enrichment = down_R24hr_GO_Liver$numDEInCat/down_R24hr_GO_Liver$expectedDEInCat
down_R24hr_GO_Liver$'-log10(FDR)' = -log10(down_R24hr_GO_Liver$over_represented_FDR)
down_R24hr_GO_Liver$description = paste(down_R24hr_GO_Liver$term," (",down_R24hr_GO_Liver$category,")", sep="")
down_R24hr_GO_Liver$group = paste("24hr Reox")
down_R24hr_GO_Liver$direction = paste("down-regulated")
down_R24hr_GO_Liver <- down_R24hr_GO_Liver[down_R24hr_GO_Liver$over_represented_FDR < 0.05, ]
down_R24hr_GO_Liver <- down_R24hr_GO_Liver[order(down_R24hr_GO_Liver$fold_enrichment),] 


# combine the results for up- and down-regulated genes:
GO_R24hr_Liver = rbind(up_R24hr_GO_Liver, down_R24hr_GO_Liver) #combine the rows, since the columns are the same



#--4.3: upregulated 24hr Reox genes Gill--

#For detailed explanations of steps, see 3.1

#Specify whether genes are upregulated, and name them
up_R24hr_genes_Gill <- factor(as.integer(all.genes %in% up_Gill_deg_24hR_only))
names(up_R24hr_genes_Gill) = all.genes
table(up_R24hr_genes_Gill)
head(up_R24hr_genes_Gill)

#Prepare data for analysis and run GO analysis
up_R24hr_pwf_Gill=nullp(up_R24hr_genes_Gill, bias.data = lengths)
up_R24hr_GO_Gill=goseq(up_R24hr_pwf_Gill, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
up_R24hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% up_Gill_deg_24hR_only, c(1,2)] 
up_R24hr_go_genes = up_R24hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the up_R24hr_GO_Gill data frame, focusing on up-regulated genes under an "24hr Reox" condition
up_R24hr_GO_Gill = merge(up_R24hr_GO_Gill, up_R24hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
up_R24hr_GO_Gill <- merge(up_R24hr_GO_Gill, go_data, by.x = "category", by.y = "id")
up_R24hr_GO_Gill$over_represented_FDR = p.adjust(up_R24hr_GO_Gill$over_represented_pvalue, method="BH")
up_R24hr_GO_Gill$under_represented_FDR = p.adjust(up_R24hr_GO_Gill$under_represented_pvalue, method="BH")
up_R24hr_GO_Gill$expectedDEInCat = (up_R24hr_GO_Gill$numInCat/length(up_R24hr_genes_Gill))*length(intersect(up_Gill_deg_24hR_only, all.genes))
up_R24hr_GO_Gill$fold_enrichment = up_R24hr_GO_Gill$numDEInCat/up_R24hr_GO_Gill$expectedDEInCat
up_R24hr_GO_Gill$'-log10(FDR)' = -log10(up_R24hr_GO_Gill$over_represented_FDR)
up_R24hr_GO_Gill$description = paste(up_R24hr_GO_Gill$term," (",up_R24hr_GO_Gill$category,")", sep="")
up_R24hr_GO_Gill$group = paste("24hr Reox")
up_R24hr_GO_Gill$direction = paste("up-regulated")
up_R24hr_GO_Gill <- up_R24hr_GO_Gill[up_R24hr_GO_Gill$over_represented_FDR < 0.05, ]
up_R24hr_GO_Gill <- up_R24hr_GO_Gill[order(up_R24hr_GO_Gill$fold_enrichment),] 



#--4.4: downregulated 24hr Reox genes Gill--

#For detailed explanations of steps, see 3.1

#Specify whether genes are downregulated, and name them
down_R24hr_genes_Gill <- factor(as.integer(all.genes %in% down_Gill_deg_24hR_only))
names(down_R24hr_genes_Gill) = all.genes
table(down_R24hr_genes_Gill)
head(down_R24hr_genes_Gill)

#Prepare data for analysis and run GO analysis
down_R24hr_pwf_Gill=nullp(down_R24hr_genes_Gill, bias.data = lengths)
down_R24hr_GO_Gill=goseq(down_R24hr_pwf_Gill, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
down_R24hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% down_Gill_deg_24hR_only, c(1,2)] 
down_R24hr_go_genes = down_R24hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the down_R24hr_GO_Gill data frame, focusing on down-regulated genes under an "24hr Reox" condition
down_R24hr_GO_Gill = merge(down_R24hr_GO_Gill, down_R24hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
down_R24hr_GO_Gill <- merge(down_R24hr_GO_Gill, go_data, by.x = "category", by.y = "id")
down_R24hr_GO_Gill$over_represented_FDR = p.adjust(down_R24hr_GO_Gill$over_represented_pvalue, method="BH")
down_R24hr_GO_Gill$under_represented_FDR = p.adjust(down_R24hr_GO_Gill$under_represented_pvalue, method="BH")
down_R24hr_GO_Gill$expectedDEInCat = (down_R24hr_GO_Gill$numInCat/length(down_R24hr_genes_Gill))*length(intersect(down_Gill_deg_24hR_only, all.genes))
down_R24hr_GO_Gill$fold_enrichment = down_R24hr_GO_Gill$numDEInCat/down_R24hr_GO_Gill$expectedDEInCat
down_R24hr_GO_Gill$'-log10(FDR)' = -log10(down_R24hr_GO_Gill$over_represented_FDR)
down_R24hr_GO_Gill$description = paste(down_R24hr_GO_Gill$term," (",down_R24hr_GO_Gill$category,")", sep="")
down_R24hr_GO_Gill$group = paste("24hr Reox")
down_R24hr_GO_Gill$direction = paste("down-regulated")
down_R24hr_GO_Gill <- down_R24hr_GO_Gill[down_R24hr_GO_Gill$over_represented_FDR < 0.05, ]
down_R24hr_GO_Gill <- down_R24hr_GO_Gill[order(down_R24hr_GO_Gill$fold_enrichment),] 


# combine the results for up- and down-regulated genes:
GO_R24hr_Gill = rbind(up_R24hr_GO_Gill, down_R24hr_GO_Gill) #combine the rows, since the columns are the same



#-----Step 5: Run GO analyses on the genes differentially regulated in 3hr Reox-----

#Runs analysis on genes that are only upregulated when compared to 3hr reox, and not any other comparison
#only in 3hr vs N, not in 24hr vs N or A vs N


#--5.1: upregulated 3hr Reox genes Liver--

#For detailed explanations of steps, see 3.1

#Specify whether genes are upregulated, and name them
up_R3hr_genes_Liver <- factor(as.integer(all.genes %in% up_Liver_deg_3hR_only))
names(up_R3hr_genes_Liver) = all.genes
table(up_R3hr_genes_Liver)
head(up_R3hr_genes_Liver)

#Prepare data for analysis and run GO analysis
up_R3hr_pwf_Liver=nullp(up_R3hr_genes_Liver, bias.data = lengths)
up_R3hr_GO_Liver=goseq(up_R3hr_pwf_Liver, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
up_R3hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% up_Liver_deg_3hR_only, c(1,2)] 
up_R3hr_go_genes = up_R3hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the up_R3hr_GO_Liver data frame, focusing on up-regulated genes under an "3hr Reox" condition
up_R3hr_GO_Liver = merge(up_R3hr_GO_Liver, up_R3hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
up_R3hr_GO_Liver <- merge(up_R3hr_GO_Liver, go_data, by.x = "category", by.y = "id")
up_R3hr_GO_Liver$over_represented_FDR = p.adjust(up_R3hr_GO_Liver$over_represented_pvalue, method="BH")
up_R3hr_GO_Liver$under_represented_FDR = p.adjust(up_R3hr_GO_Liver$under_represented_pvalue, method="BH")
up_R3hr_GO_Liver$expectedDEInCat = (up_R3hr_GO_Liver$numInCat/length(up_R3hr_genes_Liver))*length(intersect(up_Liver_deg_3hR_only, all.genes))
up_R3hr_GO_Liver$fold_enrichment = up_R3hr_GO_Liver$numDEInCat/up_R3hr_GO_Liver$expectedDEInCat
up_R3hr_GO_Liver$'-log10(FDR)' = -log10(up_R3hr_GO_Liver$over_represented_FDR)
up_R3hr_GO_Liver$description = paste(up_R3hr_GO_Liver$term," (",up_R3hr_GO_Liver$category,")", sep="")
up_R3hr_GO_Liver$group = paste("3hr Reox")
up_R3hr_GO_Liver$direction = paste("up-regulated")
up_R3hr_GO_Liver <- up_R3hr_GO_Liver[up_R3hr_GO_Liver$over_represented_FDR < 0.05, ]
up_R3hr_GO_Liver <- up_R3hr_GO_Liver[order(up_R3hr_GO_Liver$fold_enrichment),] 



#--5.2: downregulated 3hr Reox genes Liver--

#For detailed explanations of steps, see 3.1

#Specify whether genes are downregulated, and name them
down_R3hr_genes_Liver <- factor(as.integer(all.genes %in% down_Liver_deg_3hR_only))
names(down_R3hr_genes_Liver) = all.genes
table(down_R3hr_genes_Liver)
head(down_R3hr_genes_Liver)

#Prepare data for analysis and run GO analysis
down_R3hr_pwf_Liver=nullp(down_R3hr_genes_Liver, bias.data = lengths)
down_R3hr_GO_Liver=goseq(down_R3hr_pwf_Liver, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
down_R3hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% down_Liver_deg_3hR_only, c(1,2)] 
down_R3hr_go_genes = down_R3hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the down_R3hr_GO_Liver data frame, focusing on down-regulated genes under an "3hr Reox" condition
down_R3hr_GO_Liver = merge(down_R3hr_GO_Liver, down_R3hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
down_R3hr_GO_Liver <- merge(down_R3hr_GO_Liver, go_data, by.x = "category", by.y = "id")
down_R3hr_GO_Liver$over_represented_FDR = p.adjust(down_R3hr_GO_Liver$over_represented_pvalue, method="BH")
down_R3hr_GO_Liver$under_represented_FDR = p.adjust(down_R3hr_GO_Liver$under_represented_pvalue, method="BH")
down_R3hr_GO_Liver$expectedDEInCat = (down_R3hr_GO_Liver$numInCat/length(down_R3hr_genes_Liver))*length(intersect(down_Liver_deg_3hR_only, all.genes))
down_R3hr_GO_Liver$fold_enrichment = down_R3hr_GO_Liver$numDEInCat/down_R3hr_GO_Liver$expectedDEInCat
down_R3hr_GO_Liver$'-log10(FDR)' = -log10(down_R3hr_GO_Liver$over_represented_FDR)
down_R3hr_GO_Liver$description = paste(down_R3hr_GO_Liver$term," (",down_R3hr_GO_Liver$category,")", sep="")
down_R3hr_GO_Liver$group = paste("3hr Reox")
down_R3hr_GO_Liver$direction = paste("down-regulated")
down_R3hr_GO_Liver <- down_R3hr_GO_Liver[down_R3hr_GO_Liver$over_represented_FDR < 0.05, ]
down_R3hr_GO_Liver <- down_R3hr_GO_Liver[order(down_R3hr_GO_Liver$fold_enrichment),] 


# combine the results for up- and down-regulated genes:
GO_R3hr_Liver = rbind(up_R3hr_GO_Liver, down_R3hr_GO_Liver) #combine the rows, since the columns are the same



#--5.3: upregulated 3hr Reox genes Gill--

#For detailed explanations of steps, see 3.1

#Specify whether genes are upregulated, and name them
up_R3hr_genes_Gill <- factor(as.integer(all.genes %in% up_Gill_deg_3hR_only))
names(up_R3hr_genes_Gill) = all.genes
table(up_R3hr_genes_Gill)
head(up_R3hr_genes_Gill)

#Prepare data for analysis and run GO analysis
up_R3hr_pwf_Gill=nullp(up_R3hr_genes_Gill, bias.data = lengths)
up_R3hr_GO_Gill=goseq(up_R3hr_pwf_Gill, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
up_R3hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% up_Gill_deg_3hR_only, c(1,2)] 
up_R3hr_go_genes = up_R3hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the up_R3hr_GO_Gill data frame, focusing on up-regulated genes under an "3hr Reox" condition
up_R3hr_GO_Gill = merge(up_R3hr_GO_Gill, up_R3hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
up_R3hr_GO_Gill <- merge(up_R3hr_GO_Gill, go_data, by.x = "category", by.y = "id")
up_R3hr_GO_Gill$over_represented_FDR = p.adjust(up_R3hr_GO_Gill$over_represented_pvalue, method="BH")
up_R3hr_GO_Gill$under_represented_FDR = p.adjust(up_R3hr_GO_Gill$under_represented_pvalue, method="BH")
up_R3hr_GO_Gill$expectedDEInCat = (up_R3hr_GO_Gill$numInCat/length(up_R3hr_genes_Gill))*length(intersect(up_Gill_deg_3hR_only, all.genes))
up_R3hr_GO_Gill$fold_enrichment = up_R3hr_GO_Gill$numDEInCat/up_R3hr_GO_Gill$expectedDEInCat
up_R3hr_GO_Gill$'-log10(FDR)' = -log10(up_R3hr_GO_Gill$over_represented_FDR)
up_R3hr_GO_Gill$description = paste(up_R3hr_GO_Gill$term," (",up_R3hr_GO_Gill$category,")", sep="")
up_R3hr_GO_Gill$group = paste("3hr Reox")
up_R3hr_GO_Gill$direction = paste("up-regulated")
up_R3hr_GO_Gill <- up_R3hr_GO_Gill[up_R3hr_GO_Gill$over_represented_FDR < 0.05, ]
up_R3hr_GO_Gill <- up_R3hr_GO_Gill[order(up_R3hr_GO_Gill$fold_enrichment),] 



#--5.4: downregulated 3hr Reox genes Gill--

#For detailed explanations of steps, see 3.1

#Specify whether genes are downregulated, and name them
down_R3hr_genes_Gill <- factor(as.integer(all.genes %in% down_Gill_deg_3hR_only))
names(down_R3hr_genes_Gill) = all.genes
table(down_R3hr_genes_Gill)
head(down_R3hr_genes_Gill)

#Prepare data for analysis and run GO analysis
down_R3hr_pwf_Gill=nullp(down_R3hr_genes_Gill, bias.data = lengths)
down_R3hr_GO_Gill=goseq(down_R3hr_pwf_Gill, gene2cat = ccar_gene2cat)

#Make a dataframe containing only statistically significant results and Group the data by GO category
down_R3hr_go_genes=ccar2swissprotGO_long[ccar2swissprotGO_long$gene_id %in% down_Gill_deg_3hR_only, c(1,2)] 
down_R3hr_go_genes = down_R3hr_go_genes %>% group_by(GO_id) %>% summarise(ccar_gene_id = paste(unique(gene_id), collapse=',')) 

#Processing the down_R3hr_GO_Gill data frame, focusing on down-regulated genes under an "3hr Reox" condition
down_R3hr_GO_Gill = merge(down_R3hr_GO_Gill, down_R3hr_go_genes, by.x = "category", by.y = "GO_id", all = T)
down_R3hr_GO_Gill <- merge(down_R3hr_GO_Gill, go_data, by.x = "category", by.y = "id")
down_R3hr_GO_Gill$over_represented_FDR = p.adjust(down_R3hr_GO_Gill$over_represented_pvalue, method="BH")
down_R3hr_GO_Gill$under_represented_FDR = p.adjust(down_R3hr_GO_Gill$under_represented_pvalue, method="BH")
down_R3hr_GO_Gill$expectedDEInCat = (down_R3hr_GO_Gill$numInCat/length(down_R3hr_genes_Gill))*length(intersect(down_Gill_deg_3hR_only, all.genes))
down_R3hr_GO_Gill$fold_enrichment = down_R3hr_GO_Gill$numDEInCat/down_R3hr_GO_Gill$expectedDEInCat
down_R3hr_GO_Gill$'-log10(FDR)' = -log10(down_R3hr_GO_Gill$over_represented_FDR)
down_R3hr_GO_Gill$description = paste(down_R3hr_GO_Gill$term," (",down_R3hr_GO_Gill$category,")", sep="")
down_R3hr_GO_Gill$group = paste("3hr Reox")
down_R3hr_GO_Gill$direction = paste("down-regulated")
down_R3hr_GO_Gill <- down_R3hr_GO_Gill[down_R3hr_GO_Gill$over_represented_FDR < 0.05, ]
down_R3hr_GO_Gill <- down_R3hr_GO_Gill[order(down_R3hr_GO_Gill$fold_enrichment),] 


# combine the results for up- and down-regulated genes:
GO_R3hr_Gill = rbind(up_R3hr_GO_Gill, down_R3hr_GO_Gill) #combine the rows, since the columns are the same



#-----Step 6: Plotting!-----

#--6.1: Plot prep--

#Combine all liver datasets and all gill datasets separately for plotting!
GO_all_Liver = rbind(GO_AvsN_Liver, GO_R24hr_Liver, GO_R3hr_Liver) 
GO_all_Gill = rbind(GO_AvsN_Gill, GO_R24hr_Gill, GO_R3hr_Gill) 

#Export both datasets
write.table(GO_all_Liver[, -c(9,16)], "./GO_output/GO_all_Liver.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(GO_all_Gill[, -c(9,16)], "./GO_output/GO_all_Gill.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#Clean up the description column by removing instances of certain strings
#Think this was used to remove "RNA" things, so look into if we need to use them now?
#GO.all$description <- gsub("(SSU-rGill, 5.8S rGill, LSU-rGill)", "", GO.all$description, fixed = T) 
#GO.all$term <- gsub("(SSU-rGill, 5.8S rGill, LSU-rGill)", "", GO.all$term, fixed = T)

#Remove "-regulated" from up-regulated and down-regulated in the "direction" column
GO_all_Liver$direction <- gsub("-regulated", "", GO_all_Liver$direction, fixed = T) 
GO_all_Gill$direction <- gsub("-regulated", "", GO_all_Gill$direction, fixed = T)

##Filter out the top 12 results based on fold enrichment

# Sort based on fold enrichment and then by FDR (for statistical significance)
top_GO_Liver_sorted <- GO_all_Liver %>%
  arrange(desc(fold_enrichment), "-log10(FDR)")

top_GO_Gill_sorted <- GO_all_Gill %>%
  arrange(desc(fold_enrichment), "-log10(FDR)")

#Get the top 12 GO terms
top_GO_Liver = top_GO_Liver_sorted %>% slice_max(fold_enrichment, n = 12, with_ties = FALSE)
top_GO_Gill = top_GO_Gill_sorted %>% slice_max(fold_enrichment, n = 12, with_ties = FALSE)




#--6.2: Liver GO plot--

# make the plot with ggplot
Liver_GO_plot = ggplot(top_GO_Liver, aes(x = term, y = fold_enrichment)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = "azure4", linewidth = .5) + 
  geom_point(data = top_GO_Liver, aes(x = reorder(term, fold_enrichment),
                                y = log2(fold_enrichment), group = ontology, size = numDEInCat,
                                colour = `-log10(FDR)`), alpha=.8) + 
  scale_colour_viridis_c(breaks = c(0,4,8,12,16,20,24), limits = c(-1,max(top_GO_Liver$`-log10(FDR)`)+1)) + #scale_y_continuous(limits = c(0,35)) +
  coord_flip() + theme_bw() + scale_x_discrete(labels = wrap_format(40)) +
  theme(axis.ticks.y = element_blank(), axis.ticks.length=unit(-0.1, "cm"), 
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"), size = 10), axis.title.y = element_blank(),
        axis.text.y = element_text(margin=margin(0,5,5,5,"pt"), size = 10), strip.text = element_text(size = 10),
        axis.text = element_text(color = "black"), 
        panel.grid.minor = element_blank(), 
        legend.position.inside = c(0.7, 0.7), legend.text=element_text(size=10), legend.title=element_text(size=10),
        legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.5, 'cm')) +
  ylab("log2(FE)") + labs(title = "Liver") +
  facet_nested(direction + ontology ~ factor(group, levels=c("Anoxia","3hr Reox", "24hr Reox")), scales = "free_y", space = "free_y") + # this is the important function that tells ggplot to group (nest) by certain variables, so that up and down genes are separated within each column, and separated by CC, MF and BP horizontally
  labs(color="-log10(FDR)", size="Number of genes") + guides(
    color = guide_colorbar(order = 0, title.position = "top"),
    size = guide_legend(order = 1, title.position = "top")
  )

#Plot the plot
Liver_GO_plot 

#Make a pdf of the plot
cairo_pdf("./GO_output/Liver_GO_plot.pdf", width = 8, height = 7)
Liver_GO_plot
dev.off()

#Make a tiff of the plot
tiff("./GO_output/Liver_GO_plot.tif", units = "in", width = 8, height = 7, compression = "none", res = 300) 
Liver_GO_plot
dev.off()



#--6.3: Gill GO plot--

# make the plot with ggplot
Gill_GO_plot = ggplot(top_GO_Gill, aes(x = term, y = fold_enrichment)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = "azure4", linewidth = .5) + 
  geom_point(data = top_GO_Gill, aes(x = reorder(term, fold_enrichment),
                                      y = log2(fold_enrichment), group = ontology, size = numDEInCat,
                                      colour = `-log10(FDR)`), alpha=.8) + 
  scale_colour_viridis_c(breaks = c(0,2,4,6,8,10), limits = c(-1,max(top_GO_Gill$`-log10(FDR)`)+1)) + #scale_y_continuous(limits = c(0,35)) +
  coord_flip() + theme_bw() + scale_x_discrete(labels = wrap_format(40)) +
  theme(axis.ticks.y = element_blank(), axis.ticks.length=unit(-0.1, "cm"), 
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"), size = 10), axis.title.y = element_blank(),
        axis.text.y = element_text(margin=margin(0,5,5,5,"pt"), size = 10), strip.text = element_text(size = 10),
        axis.text = element_text(color = "black"), 
        panel.grid.minor = element_blank(), 
        legend.position.inside = c(0.7, 0.7), legend.text=element_text(size=10), legend.title=element_text(size=10),
        legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.5, 'cm')) +
  ylab("log2(FE)") + labs(title = "Gill") +
  facet_nested(direction + ontology ~ factor(group, levels=c("Anoxia","3hr Reox", "24hr Reox")), scales = "free_y", space = "free_y") + # this is the important function that tells ggplot to group (nest) by certain variables, so that up and down genes are separated within each column, and separated by CC, MF and BP horizontally
  labs(color="-log10(FDR)", size="Number of genes") + guides(
    color = guide_colorbar(order = 0, title.position = "top"),
    size = guide_legend(order = 1, title.position = "top")
  )

#Plot the plot
Gill_GO_plot 

#Make a pdf of the plot
cairo_pdf("./GO_output/Gill_GO_plot.pdf", width = 8, height = 7)
Gill_GO_plot
dev.off()

#Make a tiff of the plot
tiff("./GO_output/Gill_GO_plot.tif", units = "in", width = 8, height = 7, compression = "none", res = 300) 
Gill_GO_plot
dev.off()













