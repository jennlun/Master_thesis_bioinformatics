## Set Working directory
setwd("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/4_plot_venn_heat")

## Load packages
library(VennDiagram)
library(scales)
library(edgeR)
library(viridis)
library(ComplexHeatmap)
library(tidyverse)


##--Upload data files--
## Read count matrices
Liver <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/source_files_DESeq/Liver_Matrix_Count.txt")
Liver <- Liver[,c(1,20:25,14:19,8:13,2:7)] #reorder

Gill <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/source_files_DESeq/Gill_Matrix_Count.txt")
Gill <- Gill[,c(1,20:25,14:19,8:13,2:7)] #reorder

#Import the files containing all genes that are up- or downregulated
all_Liver_degs = read_lines("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_03_12_RNAseq_analysis/2025_03_13_GO_analysis/GO_output/all_Liver_degs.txt")
all_Gill_degs = read_lines("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_03_12_RNAseq_analysis/2025_03_13_GO_analysis/GO_output/all_Gill_degs.txt")

#Import the files containing toptags (top 10 genes)
toptags_Liver = read_lines("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_03_12_RNAseq_analysis/2025_03_13_GO_analysis/GO_output/toptags_Liver.txt")
toptags_Gill = read_lines("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_03_12_RNAseq_analysis/2025_03_13_GO_analysis/GO_output/toptags_Gill.txt")

## Upload Annotation file
annotation <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_03_12_RNAseq_analysis/2025_03_12_file_generation/output_file_generation/blastp_hits_swissprot_genes.txt")

##Upload DEG files from DESeq
#Upload Liver dataframes
deg_Liver_AvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_AvsN.txt")
deg_Liver_3hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_3hRvsN.txt")
deg_Liver_24hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Liver_24hRvsN.txt")

#Upload Gill dataframes
deg_Gill_AvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_AvsN.txt")
deg_Gill_3hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_3hRvsN.txt")
deg_Gill_24hRvsN <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis/output_DESeq/deg_Gill_24hRvsN.txt")


#-----Step 1: Plot venn diagrams between treatment groups-----

#Look at the overlap between A, 3hR and 24hR

#Create a vector containing the names of our two tissue types
setnames <- c("Liver", "Gill")

#Create a vector containing the names of our DEGs for both tissues
listAvsN <- c("deg_Liver_AvsN", "deg_Gill_AvsN")
list3hRvsN <- c("deg_Liver_3hRvsN", "deg_Gill_3hRvsN")
list24hRvsN <- c("deg_Liver_24hRvsN", "deg_Gill_24hRvsN")

#Create a function for making the vennplots
vennplot <- function(listAvsN, list3hRvsN, list24hRvsN, setname){
  
  # get data
  
  setAvsN=get(listAvsN)
  set3hRvsN=get(list3hRvsN)
  set24hRvsN=get(list24hRvsN)
  
  #generate venn plot for upregulated genes
  
  vennplotup=venn.diagram(
    x = list(rownames(setAvsN)[which(setAvsN$padj < 0.05 & setAvsN$log2FoldChange > 0)], rownames(set3hRvsN)[which(set3hRvsN$padj < 0.05 &  set3hRvsN$log2FoldChange > 0)], rownames(set24hRvsN)[which(set24hRvsN$padj < 0.05 & set24hRvsN$log2FoldChange > 0)]),
    disable.logging = T,
    category.names = c("AvsN" , "3hRvsN" , "24hRvsN"),
    filename = NULL,
    lwd = 1,
    col=c('black', 'black', 'black'),
    fill = c(alpha("black",0.3), alpha('black',0.3), alpha('black',0.3)),
    cex = 0.8,
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.default.pos = "outer",
    cat.pos = c(-15, 15, 180),
    cat.dist = c(0.08, 0.08, 0.08),
    cat.fontfamily = "sans",
    cat.col = c('black', 'black', 'black'),
    rotation = 1,
    margin=0.1,
    euler.d=F,
    scaled=F)
  
  #save plot
  tiff(filename = paste0("./vennplot_up_",setname,".tiff"), width = 2, height = 2, units = "in", res = 300)
  print(grid.draw(vennplotup))
  dev.off()
  
  #generate venn plot for downregulated genes
  
  vennplotdown=venn.diagram(
    x = list(rownames(setAvsN)[which(setAvsN$padj < 0.05 & setAvsN$log2FoldChange < 0)], rownames(set3hRvsN)[which(set3hRvsN$padj < 0.05 &  set3hRvsN$log2FoldChange < 0)], rownames(set24hRvsN)[which(set24hRvsN$padj < 0.05 & set24hRvsN$log2FoldChange < 0)]),
    disable.logging = T,
    category.names = c("AvsN" , "3hRvsN" , "24hRvsN"),
    filename = NULL,
    lwd = 1,
    col=c('black', 'black', 'black'),
    fill = c(alpha("black",0.3), alpha('black',0.3), alpha('black',0.3)),
    cex = 0.8,
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.default.pos = "outer",
    cat.pos = c(-15, 15, 180),
    cat.dist = c(0.08, 0.08, 0.08),
    cat.fontfamily = "sans",
    cat.col = c('black', 'black', 'black'),
    rotation = 1,
    margin=0.1,
    euler.d=F,
    scaled=F)
  
  #save plot
  tiff(filename = paste0("./vennplot_down_",setname,".tiff"), width = 2, height = 2, units = "in", res = 300)
  print(grid.draw(vennplotdown))
  dev.off()
  
}

#Run function with our data to generate vennplot
mapply(vennplot, listAvsN, list3hRvsN, list24hRvsN, setnames)



#-----Step 2: Plot heat maps-----

#--2.1: Normalize datasets to prep for plotting--

#--2.1.1: Process count matrices--
#use the first column as rownames
rownames(Liver) <- Liver[,1] 
#remove the first column
Liver <- Liver[,-1] 

#filter counts so that there are at least 6 samples with non-zero counts for each gene
Liver_filt <- Liver[rowSums(cpm(Liver) > 0) >= 6, ] 


Gill_filt <- Gill

#filter the gill df to only include the rows where the gene ID matches the rownames in the liver_filt df
Gill_filt <- Gill[Gill$Geneid %in% row.names(Liver_filt),] 

#use the first column as rownames
rownames(Gill_filt) <- Gill_filt[,1] 
#remove the first column
Gill_filt <- Gill_filt[,-1] 

#use the first column as rownames
rownames(Gill) <- Gill[,1] 
#remove the first column
Gill <- Gill[,-1] 



#--2.1.2: Use TMM normalization on filtered Liver counts--
#Create a DGEList object (a data structure that stores information related to differential gene expression analysis. 
#Typically includes count data, sample information, and other metadata required for the analysis
Liver_tmm_filt <- DGEList(counts=Liver_filt, group = c("A_Normoxia","A_Normoxia","A_Normoxia","A_Normoxia",
                                                       "A_Normoxia","A_Normoxia","B_Anoxia","B_Anoxia","B_Anoxia",
                                                       "B_Anoxia","B_Anoxia","B_Anoxia","C_3h_Reox","C_3h_Reox",
                                                       "C_3h_Reox","C_3h_Reox","C_3h_Reox","C_3h_Reox",
                                                       "D_24h_Reox","D_24h_Reox","D_24h_Reox","D_24h_Reox",
                                                       "D_24h_Reox","D_24h_Reox"))

#calculate normalization factors using the TMM method to adjust the raw counts for differences in library composition across samples
Liver_tmm_filt <- calcNormFactors(Liver_tmm_filt, method = "TMM" )

#convert raw count data to counts per million (cpm)
#helps standardize gene expression measurements across samples by controlling for sequencing depth differences, 
#turning raw counts into something more comparable across differing library sizes
Liver_tmm_filt <- cpm(Liver_tmm_filt)

#Turn matrix into df
Liver_tmm_filt <- as.data.frame(Liver_tmm_filt)

#Use the rownames to create a new column calles gene_id
Liver_tmm_filt$gene_id <- row.names(Liver_tmm_filt)

#Remove rownames
row.names(Liver_tmm_filt)=c()

#Add information from the annotation file using the gene IDs as identifiers
Liver_tmm_filt = merge(annotation, Liver_tmm_filt, by.x="gene_id", by.y="gene_id", all.y = T)

#Export table of filtered TMM count matrix:
write.table(Liver_tmm_filt, "./output_files/Liver_genes_filtered_TMM.matrix", row.names=FALSE, 
            col.names=TRUE, sep="\t", quote = F)  ### note that this will be the rounded counts



#--2.1.3: Use TMM normalization on filtered Gill counts--
Gill_tmm_filt <- DGEList(counts=Gill_filt, group = c("A_Normoxia","A_Normoxia","A_Normoxia","A_Normoxia",
                                                     "A_Normoxia","A_Normoxia","B_Anoxia","B_Anoxia","B_Anoxia",
                                                     "B_Anoxia","B_Anoxia","B_Anoxia","C_3h_Reox","C_3h_Reox",
                                                     "C_3h_Reox","C_3h_Reox","C_3h_Reox","C_3h_Reox",
                                                     "D_24h_Reox","D_24h_Reox","D_24h_Reox","D_24h_Reox",
                                                     "D_24h_Reox","D_24h_Reox"))
Gill_tmm_filt <- calcNormFactors(Gill_tmm_filt, method = "TMM" )
Gill_tmm_filt <- cpm(Gill_tmm_filt)
Gill_tmm_filt <- as.data.frame(Gill_tmm_filt)
Gill_tmm_filt$gene_id <- row.names(Gill_tmm_filt)
row.names(Gill_tmm_filt)=c()
Gill_tmm_filt = merge(annotation, Gill_tmm_filt, by.x="gene_id", by.y="gene_id", all.y = T)

#Export table of filtered TMM count matrix:
write.table(Gill_tmm_filt, "./output_files/Gill_genes_filtered_TMM.matrix", 
            row.names=FALSE, col.names=TRUE, sep="\t", quote = F)  ### note that this will be the rounded counts



#--2.1.4: Use TMM normalization on full Liver and gill counts--

# TMM normalisation on Liver for all genes:
Liver_tmm <- DGEList(counts=Liver, group = c("A_Normoxia","A_Normoxia","A_Normoxia","A_Normoxia",
                                             "A_Normoxia","A_Normoxia","B_Anoxia","B_Anoxia","B_Anoxia",
                                             "B_Anoxia","B_Anoxia","B_Anoxia","C_3h_Reox","C_3h_Reox",
                                             "C_3h_Reox","C_3h_Reox","C_3h_Reox","C_3h_Reox",
                                             "D_24h_Reox","D_24h_Reox","D_24h_Reox","D_24h_Reox",
                                             "D_24h_Reox","D_24h_Reox"))
Liver_tmm <- calcNormFactors(Liver_tmm, method = "TMM" )
Liver_tmm <- cpm(Liver_tmm)
Liver_tmm <- as.data.frame(Liver_tmm)
Liver_tmm$gene_id <- row.names(Liver_tmm)
row.names(Liver_tmm)=c()
Liver_tmm = merge(annotation, Liver_tmm, by.x="gene_id", by.y="gene_id", all.y = T)

#Export table of non-filtered TMM count matrix:
write.table(Liver_tmm, "./output_files/Liver_genes_nonfiltered_TMM.matrix", 
            row.names=FALSE, col.names=TRUE, sep="\t", quote = F)  ### note that this will be the rounded counts


# TMM normalization of Gill counts for all genes:
Gill_tmm <- DGEList(counts=Gill, group = c("A_Normoxia","A_Normoxia","A_Normoxia","A_Normoxia",
                                           "A_Normoxia","A_Normoxia","B_Anoxia","B_Anoxia","B_Anoxia",
                                           "B_Anoxia","B_Anoxia","B_Anoxia","C_3h_Reox","C_3h_Reox",
                                           "C_3h_Reox","C_3h_Reox","C_3h_Reox","C_3h_Reox",
                                           "D_24h_Reox","D_24h_Reox","D_24h_Reox","D_24h_Reox",
                                           "D_24h_Reox","D_24h_Reox"))
Gill_tmm <- calcNormFactors(Gill_tmm, method = "TMM" )
Gill_tmm <- cpm(Gill_tmm)
Gill_tmm <- as.data.frame(Gill_tmm)
Gill_tmm$gene_id <- row.names(Gill_tmm)
row.names(Gill_tmm)=c()
Gill_tmm = merge(annotation, Gill_tmm, by.x="gene_id", by.y="gene_id", all.y = T)


#Export table of non-filtered TMM count matrix:
write.table(Gill_tmm, "./output_files/Gill_genes_nonfiltered_TMM.matrix", 
            row.names=FALSE, col.names=TRUE, sep="\t", quote = F)  ### note that this will be the rounded counts



#--2.1.4: make TMM dataset with all Liver and Gill toptags--

#Get the unique toptags from each toptag dataset
all_toptags <- unique(as.character(c(toptags_Liver, toptags_Gill)))

#Get the toptag gene names 
all_toptags_names <- Liver_tmm_filt[Liver_tmm_filt$gene_id %in% all_toptags, c(1:2)]

#Save full toptag names
write.table(all_toptags_names, "./output_files/toptags_names.txt", row.names = F, 
            col.names = T, quote = F, sep = "\t") 

# Remove the section behind the underscore
all_toptags_names$symbol <- sub("_.*", "", all_toptags_names$symbol)

#generate table that can be edited in excel (add defining symbol to duplicate names)
#remove # if changing parameters
#write.table(all_toptags_names, "./output_files/toptags_names_modified.txt", row.names = F, col.names = T, quote = F, sep = "\t") 


#import table with edited gene names
all_toptags_names <- read.table("./output_files/toptags_names_modified.txt", 
                                header = T, na.strings="", sep="\t", quote = "")

#keep only the rows from Tmm liver file where the gene ID is present in the toptags file
all_toptags_Liver_tmm <- Liver_tmm[Liver_tmm$gene_id %in% toptags_Liver, ]

#Add gene names from toptag name file
all_toptags_Liver_tmm <- merge(all_toptags_names, all_toptags_Liver_tmm, by.x = "gene_id", by.y = "gene_id")

#repeat for Gill
all_toptags_Gill_tmm <- Gill_tmm[Gill_tmm$gene_id %in% toptags_Gill, ]
all_toptags_Gill_tmm <- merge(all_toptags_names, all_toptags_Gill_tmm, by.x = "gene_id", by.y = "gene_id")

#Save tables with toptag names and info
write.table(all_toptags_Liver_tmm, "./output_files/all_toptags_Liver_tmm.txt", 
            row.names = F, col.names = T, quote = F, sep = "\t")
write.table(all_toptags_Gill_tmm, "./output_files/all_toptags_Gill_tmm.txt", 
            row.names = F, col.names = T, quote = F, sep = "\t")



#--2.2: Plot heatmap for all differentially expressed Liver genes--

#Create a color palette consisting of 20 colors
#PURPLE IS THE BEST COLOR
mypalette <- viridis_pal(option="magma")(20)

#Filter the full Liver tmm file to only keep rows where the gene ID is present in all DEGs from the liver samples
all_Liver_degs_tmm=Liver_tmm[Liver_tmm$gene_id %in% all_Liver_degs,]

#Remove duplicate rows so that each gene ID is only present once
all_Liver_degs_tmm=all_Liver_degs_tmm[!duplicated(all_Liver_degs_tmm$gene_id),]

#Name rows based on the gene_id column
rownames(all_Liver_degs_tmm) = all_Liver_degs_tmm$gene_id

#create a matrix containing only columns 4-27 of the full dataframe
all_Liver_degs_matrix=as.matrix(all_Liver_degs_tmm[,c(4:27)])

#Convert all data to numeric and maintain structure
#Ensures that we only have numeric entries in the matrix
all_Liver_degs_matrix=matrix(as.numeric(all_Liver_degs_matrix), ncol = ncol(all_Liver_degs_matrix))

#Name rows based on the gene_id column in the all degs df
row.names(all_Liver_degs_matrix)=all_Liver_degs_tmm$gene_id

#Scale and transpose data to prepare for plotting
scaled_all_Liver_degs = scale(t(all_Liver_degs_matrix))

#reversal of transposition to properly align the matrix
#genes in rows and samples in columns etc
scaled_all_Liver_degs = t(scaled_all_Liver_degs)


#Generate heatmap
Liver_heat_all=Heatmap(scaled_all_Liver_degs, col = mypalette, border = TRUE, cluster_columns = F, 
                       show_row_dend = F, row_km = 4, show_row_names = F, show_column_names = F, 
                       row_title = NULL, column_names_rot = 0, column_title_gp = gpar(fontsize = 10, fontface = "bold"),  
                       heatmap_legend_param = list(title_gp = gpar(fontface = "plain", fontsize = 10), 
                                                   legend_label_gp=gpar(fontsize=10), title = "z-score", 
                                                   labels_gp=gpar(fontsize=10), direction = "horizontal"), 
                       column_split = factor(c(rep("N", 6), rep("A", 6), rep("3hR", 6), rep("24hR", 6)), 
                                             levels = c("N", "A", "3hR", "24hR")), column_names_gp = gpar(fontsize = 10), 
                       cluster_column_slices = F, column_order = sort(colnames(scaled_all_Liver_degs)))

#render heatmap
Liver_heat_all <- draw(Liver_heat_all, heatmap_legend_side = "top")

#retrieve the order in which rows (genes) were arranged in the heatmap after any clustering operations had been applied
#Can be used in further analysis, ex clustering
clusters=row_order(Liver_heat_all)

#Create objects containing rows from different clusters
cluster1=all_Liver_degs_tmm[clusters$`1`,]
cluster2=all_Liver_degs_tmm[clusters$`2`,]
cluster3=all_Liver_degs_tmm[clusters$`3`,]
cluster4=all_Liver_degs_tmm[clusters$`4`,]

#Export heatmap
pdf("./Liver_heatmap_all.pdf", width = (17/2.56), height = (10/2.56))  # Adjust width and height as needed (unit in inches)
draw(Liver_heat_all)
dev.off()


tiff("./Liver_heatmap_all.tif", width = 4.5, height = 6,units = "in", res =300)
Liver_heat_all
dev.off()



#--2.3: Plot heatmap for the top 10 differentially expressed Liver genes--

#Create a matrix containing only genes that are specified in the toptags liver file
toptags_Liver_matrix <- all_toptags_Liver_tmm[all_toptags_Liver_tmm$gene_id %in% toptags_Liver, ]

#Assign the gene ids as row names
rownames(toptags_Liver_matrix) = as.vector(toptags_Liver_matrix$symbol.x) 

#Make a matrix with the top tags
toptags_Liver_matrix=as.matrix(toptags_Liver_matrix[,c(5:28)]) 

#Scale and transpose data to prepare for plotting
scaled_top_Liver_degs = scale(t(toptags_Liver_matrix))
scaled_top_Liver_degs = t(scaled_top_Liver_degs)

#Generate heatmap
Liver_heat_top = Heatmap(scaled_top_Liver_degs, col = mypalette, cluster_columns = F, 
                         cluster_rows = T, show_column_dend = F, show_row_dend = T, row_km = 4, 
                         show_row_names = T, row_title = NULL, 
                         column_title = "Liver",
                         column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                         row_names_gp = gpar(fontface = "italic", fontsize = 8), show_column_names = F, 
                         heatmap_legend_param = list(title_gp = gpar(fontface = "plain", fontsize = 10), 
                                                     legend_label_gp=gpar(fontsize=10), title = "z-score", 
                                                     labels_gp = gpar(fontsize = 10), 
                                                     direction = "horizontal"), 
                         column_split = factor(c(rep("N", 6), rep("A", 6),rep("3hR", 6),rep("24hR", 6)), 
                                               levels = c("N","A","3hR","24hR")), 
                         rect_gp = gpar(col = "grey", lwd = 0.5), 
                         border_gp = gpar(col = "white", lty = 0.5))



#render heatmap
Liver_heat_top <- draw(Liver_heat_top, heatmap_legend_side = "top")
Liver_heat_top

#Export heatmap
pdf("./Liver_heat_top.pdf", width = (17/2.56), height = (14/2.56))  # Adjust width and height as needed (unit in inches)
draw(Liver_heat_top)
dev.off()

tiff("./Liver_heat_top.tif", width = 3.7, height = 10.5,units = "in", res =300)
Liver_heat_top
dev.off()



#--2.4: Plot heatmap for all differentially expressed Gill genes--

#Create a color palette consisting of 20 colors
mypalette <- viridis_pal(option="magma")(20)

#Filter the full Gill tmm file to only keep rows where the gene ID is present in all DEGs from the liver samples
all_Gill_degs_tmm=Gill_tmm[Gill_tmm$gene_id %in% all_Gill_degs,]

#Remove duplicate rows so that each gene ID is only present once
all_Gill_degs_tmm=all_Gill_degs_tmm[!duplicated(all_Gill_degs_tmm$gene_id),]

#Name rows based on the gene_id column
rownames(all_Gill_degs_tmm) = all_Gill_degs_tmm$gene_id

#create a matrix containing only columns 4-27 of the full dataframe
all_Gill_degs_matrix=as.matrix(all_Gill_degs_tmm[,c(4:27)])

#Convert all data to numeric and maintain structure
#Ensures that we only have numeric entries in the matrix
all_Gill_degs_matrix=matrix(as.numeric(all_Gill_degs_matrix), ncol = ncol(all_Gill_degs_matrix))

#Name rows based on the gene_id column in the all degs df
row.names(all_Gill_degs_matrix)=all_Gill_degs_tmm$gene_id

#Scale and transpose data to prepare for plotting
scaled_all_Gill_degs = scale(t(all_Gill_degs_matrix))

#reversal of transposition to properly align the matrix
#genes in rows and samples in columns etc
scaled_all_Gill_degs = t(scaled_all_Gill_degs)


#Generate heatmap
Gill_heat_all=Heatmap(scaled_all_Gill_degs, col = mypalette, border = TRUE, cluster_columns = F, 
                       show_row_dend = F, row_km = 4, show_row_names = F, show_column_names = F, 
                       row_title = NULL, column_names_rot = 0, column_title_gp = gpar(fontsize = 10, fontface = "bold"),  
                       heatmap_legend_param = list(title_gp = gpar(fontface = "plain", fontsize = 10), 
                                                   legend_label_gp=gpar(fontsize=10), title = "z-score", 
                                                   labels_gp=gpar(fontsize=10), direction = "horizontal"), 
                       column_split = factor(c(rep("N", 6), rep("A", 6), rep("3hR", 6), rep("24hR", 6)), 
                                             levels = c("N", "A", "3hR", "24hR")), column_names_gp = gpar(fontsize = 10), 
                       cluster_column_slices = F, column_order = sort(colnames(scaled_all_Gill_degs)))

#render heatmap
Gill_heat_all <- draw(Gill_heat_all, heatmap_legend_side = "top")

#retrieve the order in which rows (genes) were arranged in the heatmap after any clustering operations had been applied
#Can be used in further analysis, ex clustering
clusters=row_order(Gill_heat_all)

#Create objects containing rows from different clusters
cluster1=all_Gill_degs_tmm[clusters$`1`,]
cluster2=all_Gill_degs_tmm[clusters$`2`,]
cluster3=all_Gill_degs_tmm[clusters$`3`,]
cluster4=all_Gill_degs_tmm[clusters$`4`,]

#Export heatmap
pdf("./Gill_heatmap_all.pdf", width = (17/2.56), height = (10/2.56))  # Adjust width and height as needed (unit in inches)
draw(Gill_heat_all)
dev.off()


tiff("./Gill_heatmap_all.tif", width = 4.5, height = 6,units = "in", res =300)
Gill_heat_all
dev.off()



#--2.5: Plot heatmap for the top 10 differentially expressed Gill genes--

#Create a matrix containing only genes that are specified in the toptags liver file
toptags_Gill_matrix <- all_toptags_Gill_tmm[all_toptags_Gill_tmm$gene_id %in% toptags_Gill, ]

#Assign the gene ids as row names
rownames(toptags_Gill_matrix) = as.vector(toptags_Gill_matrix$symbol.x) 

#Make a matrix with the top tags
toptags_Gill_matrix=as.matrix(toptags_Gill_matrix[,c(5:28)]) 

#Scale and transpose data to prepare for plotting
scaled_top_Gill_degs = scale(t(toptags_Gill_matrix))
scaled_top_Gill_degs = t(scaled_top_Gill_degs)

#Generate heatmap
Gill_heat_top = Heatmap(scaled_top_Gill_degs, col = mypalette, cluster_columns = F, 
                         cluster_rows = T, show_column_dend = F, show_row_dend = T, row_km = 5, 
                         show_row_names = T, row_title = NULL, 
                        column_title = "Gill",
                         column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                         row_names_gp = gpar(fontface = "italic", fontsize = 8), show_column_names = F, 
                         heatmap_legend_param = list(title_gp = gpar(fontface = "plain", fontsize = 10), 
                                                     legend_label_gp=gpar(fontsize=10), title = "z-score", 
                                                     labels_gp = gpar(fontsize = 10), 
                                                     direction = "horizontal"), 
                         column_split = factor(c(rep("N", 6), rep("A", 6),rep("3hR", 6),rep("24hR", 6)), 
                                               levels = c("N","A","3hR","24hR")), 
                         rect_gp = gpar(col = "grey", lwd = 0.5), 
                         border_gp = gpar(col = "white", lty = 0.5))



#render heatmap
Gill_heat_top <- draw(Gill_heat_top, heatmap_legend_side = "top")
Gill_heat_top

#Export heatmap
pdf("./Gill_heat_top.pdf", width = (17/2.56), height = (14/2.56))  # Adjust width and height as needed (unit in inches)
draw(Gill_heat_top)
dev.off()

tiff("./Gill_heat_top.tif", width = 3.7, height = 10.5,units = "in", res =300)
Gill_heat_top
dev.off()




