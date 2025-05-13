## Set Working directory
setwd("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/2_DESeq_analysis")


## Load packages
library(DESeq2)
library(ggplot2)
library(apeglm)
library(dplyr)
library(edgeR)
library(viridis)
library(gridExtra)


## Upload count matrices and sample information
Liver_file <- "./source_files_DESeq/Liver_Matrix_Count.txt"

Gill_file <- "./source_files_DESeq/Gill_Matrix_Count.txt"

data_file <- "./source_files_DESeq/sample_info_3hr.txt"


## Upload Annotation file
annotation <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/1_file_generation/output_file_generation/blastp_hits_swissprot_genes.txt")

## Read count matrices
Liver <- read.delim(Liver_file)
Liver <- Liver[,c(1,20:25,14:19,8:13,2:7)] #reorder 

Gill <- read.delim(Gill_file)
Gill <- Gill[,c(1,20:25,14:19,8:13,2:7)] #reorder


## Process count matrices
rownames(Liver) <- Liver[,1] #use the first column as rownames 
Liver <- Liver[,-1] #remove the first column

#filter counts so that there are at least 6 samples with non-zero counts for each gene
Liver_filt <- Liver[rowSums(cpm(Liver) > 0) >= 6, ] 


Gill_filt <- Gill

#filter the gill df to only include the rows where the gene ID matches the rownames in the liver_filt df
Gill_filt <- Gill[Gill$Geneid %in% row.names(Liver_filt),] 


rownames(Gill_filt) <- Gill_filt[,1] #use the first column as rownames
Gill_filt <- Gill_filt[,-1] #remove the first column

rownames(Gill) <- Gill[,1] #use the first column as rownames
Gill <- Gill[,-1] #remove the first column


## Merge count matrices
merge <- cbind(Liver_filt,Gill_filt)

head(merge) #show the first 20 lines of the merge df


## Upload sample information file
coldata <- read.delim(data_file)
#turn each column into a factor and then convert the results back into a df
coldata <- as.data.frame(apply(coldata,2,as.factor)) 

head(coldata) #show first 6 lines of df


#-------------------------------------
#Start DESeq2 analysis!
#-------------------------------------

#-----Step 1: Create DESeq2 Object for Liver -----

#Filter for rows with tissue = Liver (create index of row numbers)
ind_Liver = which(coldata$Tissue == "Liver")

#Keep the rows with the index numbers
coldata_Liver = coldata[ind_Liver,]

## Create DESeq object containing the liver tissue data 
#Liver_filt means there are at least 6 samples with non-zero counts for each gene
#coldata_liver means rows from sample info file with tissue = liver
deseq_object_Liver <- DESeqDataSetFromMatrix(countData = Liver_filt,
                                       colData = coldata_Liver, design =~ Condition)

#Normalization step (adjusts raw counts to make them comparable)
deseq_object_Liver <- estimateSizeFactors(deseq_object_Liver)


#Create matrix of normalised counts for comparison using TMM normalisation:
#normalize: look at patterns for each gene/row (up, down, etc) rather than exact values, since some genes have more counts than others
deseq_object_Liver_normalized <- counts(deseq_object_Liver, normalized=TRUE)

#Create a vector of gene ID names
Liver_geneid_vector=as.vector(row.names(deseq_object_Liver_normalized))

#Adding a column to the beginning containing the gene IDs from the vector
deseq_object_Liver_normalized=as.data.frame(cbind(Liver_geneid_vector, deseq_object_Liver_normalized))

#Remove rownames and do index numbers instead 
row.names(deseq_object_Liver_normalized)=c()

#Merge the normalized deseq liver object with the annotation file based on the gene IDs
deseq_object_Liver_normalized = merge(annotation, deseq_object_Liver_normalized, by.x="gene_id", 
                                      by.y="Liver_geneid_vector", all.y = T)

#Save the normalized Liver DESEq object
write.table(deseq_object_Liver_normalized, file="./output_DESeq/Liver_deseq2_normalised.matrix", row.names=F, sep="\t", quote=F, col.names=TRUE)



#-----Step 2: Run DESeq2 For Liver-----

#--2.1: Run DESeq using normoxia as refence level--

##"Make" a new DESeq object using the complete Liver DESeq object
deseq_object_Liver_vsN <- deseq_object_Liver

#Set reference level for Condition; we start with Normoxia
#Called A_Normoxia to set normoxia as first position in treatment order
deseq_object_Liver_vsN$Condition = relevel(deseq_object_Liver_vsN$Condition,"A_Normoxia") 

## Run DESeq function (comparing it to all other conditions)
deseq_object_Liver_vsN <- DESeq(deseq_object_Liver_vsN)

#Show comparisons from DESeq (what was compared against what)
resultsNames(deseq_object_Liver_vsN)


#--2.2: Working with output from DESeq (Anoxia vs Normoxia)--

#Filtering so that we only keep the results that show significant changes (p < 0.05)
#and logfold change >= log2(1.2) (increase/decrease of ~26%)
#Also only compare anoxia vs normoxia
#deg = differentially expressed genes
deg_Liver_AvsN <- results(deseq_object_Liver_vsN, alpha = 0.05, lfcThreshold = log2(1.2),  
                          contrast=c("Condition","B_Anoxia","A_Normoxia"))

#Shrinkage is applied to reduce the variance of log fold changes for genes with low count numbers
deg_Liver_AvsN <- lfcShrink(deseq_object_Liver_vsN, coef="Condition_B_Anoxia_vs_A_Normoxia",
                            res=deg_Liver_AvsN, type='apeglm')

#Convert the DESeq results into a dataframe 
deg_Liver_AvsN_df=as.data.frame(deg_Liver_AvsN)

#Create a column of gene IDs using the rownames
deg_Liver_AvsN_df$gene_id<-as.vector(row.names(deg_Liver_AvsN_df))

#Move gene ID column to first position
deg_Liver_AvsN_df=deg_Liver_AvsN_df[, c(6, 1:5)]

#Remove the rownames (index numbers instead, Maggy happy)
row.names(deg_Liver_AvsN_df)=c()

#Merge with annotation file based on gene ID (to add protein product and symbols)
deg_Liver_AvsN_df = merge(annotation, deg_Liver_AvsN_df, by.x="gene_id", by.y="gene_id", all.y =T)

#Sort the dataframe based on padj (smallest first)
deg_Liver_AvsN_df <- deg_Liver_AvsN_df[ order( deg_Liver_AvsN_df$padj ), ]

#Export dataframe
write.table(deg_Liver_AvsN_df,"./output_DESeq/deg_Liver_AvsN.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--2.3: Working with output from DESeq (3hr Reox vs Normoxia)--

#For descriptions, see 2.2

deg_Liver_3hRvsN <- results(deseq_object_Liver_vsN, alpha = 0.05, lfcThreshold = log2(1.2),  
                            contrast=c("Condition","C_Reoxygenation3h","A_Normoxia"))
deg_Liver_3hRvsN <- lfcShrink(deseq_object_Liver_vsN, coef="Condition_C_Reoxygenation3h_vs_A_Normoxia",
                              res=deg_Liver_3hRvsN, type='apeglm')
deg_Liver_3hRvsN_df=as.data.frame(deg_Liver_3hRvsN)
deg_Liver_3hRvsN_df$gene_id<-as.vector(row.names(deg_Liver_3hRvsN_df))
deg_Liver_3hRvsN_df=deg_Liver_3hRvsN_df[, c(6, 1:5)]
row.names(deg_Liver_3hRvsN_df)=c()
deg_Liver_3hRvsN_df = merge(annotation, deg_Liver_3hRvsN_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Liver_3hRvsN_df <- deg_Liver_3hRvsN_df[ order( deg_Liver_3hRvsN_df$padj ), ]

#Export dataframe
write.table(deg_Liver_3hRvsN_df,"./output_DESeq/deg_Liver_3hRvsN.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--2.4: Working with output from DESeq (24hr Reox vs Normoxia)--

#For descriptions, see 2.2

deg_Liver_24hRvsN <- results(deseq_object_Liver_vsN, alpha = 0.05, lfcThreshold = log2(1.2),  
                             contrast=c("Condition","D_Reoxygenation24h","A_Normoxia"))
deg_Liver_24hRvsN <- lfcShrink(deseq_object_Liver_vsN, coef="Condition_D_Reoxygenation24h_vs_A_Normoxia",
                               res=deg_Liver_24hRvsN, type='apeglm')
deg_Liver_24hRvsN_df=as.data.frame(deg_Liver_24hRvsN)
deg_Liver_24hRvsN_df$gene_id<-as.vector(row.names(deg_Liver_24hRvsN_df))
deg_Liver_24hRvsN_df=deg_Liver_24hRvsN_df[, c(6, 1:5)]
row.names(deg_Liver_24hRvsN_df)=c()
deg_Liver_24hRvsN_df = merge(annotation, deg_Liver_24hRvsN_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Liver_24hRvsN_df <- deg_Liver_24hRvsN_df[ order( deg_Liver_24hRvsN_df$padj ), ]

#Export dataframe
write.table(deg_Liver_24hRvsN_df,"./output_DESeq/deg_Liver_24hRvsN.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--2.5: Run DESeq using anoxia as refence level--

##"Make" a new DESeq object using the complete Liver DESeq object
deseq_object_Liver_vsA <- deseq_object_Liver

#Set reference level for Condition; now we're using Anoxia
deseq_object_Liver_vsA$Condition = relevel(deseq_object_Liver_vsA$Condition,"B_Anoxia")

## Run DESeq function (comparing it to all other conditions)
deseq_object_Liver_vsA <- DESeq(deseq_object_Liver_vsA)



#--2.6: Working with output from DESeq (3hr Reox vs Anoxia)--

#For descriptions, see 2.2

deg_Liver_3hRvsA <- results(deseq_object_Liver_vsA, alpha = 0.05, lfcThreshold = log2(1.2),  
                            contrast=c("Condition","C_Reoxygenation3h","B_Anoxia"))
deg_Liver_3hRvsA <- lfcShrink(deseq_object_Liver_vsA, coef="Condition_C_Reoxygenation3h_vs_B_Anoxia",
                              res=deg_Liver_3hRvsA, type='apeglm')
deg_Liver_3hRvsA_df=as.data.frame(deg_Liver_3hRvsA)
deg_Liver_3hRvsA_df$gene_id<-as.vector(row.names(deg_Liver_3hRvsA_df))
deg_Liver_3hRvsA_df=deg_Liver_3hRvsA_df[, c(6, 1:5)]
row.names(deg_Liver_3hRvsA_df)=c()
deg_Liver_3hRvsA_df = merge(annotation, deg_Liver_3hRvsA_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Liver_3hRvsA_df <- deg_Liver_3hRvsA_df[ order( deg_Liver_3hRvsA_df$padj ), ]

#Export dataframe
write.table(deg_Liver_3hRvsA_df,"./output_DESeq/deg_Liver_3hRvsA.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--2.7: Working with output from DESeq (24hr Reox vs Anoxia)--

#For descriptions, see 2.2

deg_Liver_24hRvsA <- results(deseq_object_Liver_vsA, alpha = 0.05, lfcThreshold = log2(1.2),  
                             contrast=c("Condition","D_Reoxygenation24h","B_Anoxia"))
deg_Liver_24hRvsA <- lfcShrink(deseq_object_Liver_vsA, coef="Condition_D_Reoxygenation24h_vs_B_Anoxia",
                               res=deg_Liver_24hRvsA, type='apeglm')
deg_Liver_24hRvsA_df=as.data.frame(deg_Liver_24hRvsA)
deg_Liver_24hRvsA_df$gene_id<-as.vector(row.names(deg_Liver_24hRvsA_df))
deg_Liver_24hRvsA_df=deg_Liver_24hRvsA_df[, c(6, 1:5)]
row.names(deg_Liver_24hRvsA_df)=c()
deg_Liver_24hRvsA_df = merge(annotation, deg_Liver_24hRvsA_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Liver_24hRvsA_df <- deg_Liver_24hRvsA_df[ order( deg_Liver_24hRvsA_df$padj ), ]

#Export dataframe
write.table(deg_Liver_24hRvsA_df,"./output_DESeq/deg_Liver_24hRvsA.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--2.8: Run DESeq using 3h reox as refence level--

##"Make" a new DESeq object using the complete Liver DESeq object
deseq_object_Liver_vs3hR <- deseq_object_Liver

#Set reference level for Condition; now we're using 3hr reox
deseq_object_Liver_vs3hR$Condition = relevel(deseq_object_Liver_vs3hR$Condition,"C_Reoxygenation3h")

## Run DESeq function (comparing it to all other conditions)
deseq_object_Liver_vs3hR <- DESeq(deseq_object_Liver_vs3hR)



#--2.9: Working with output from DESeq (24hr Reox vs 3hr Reox)--

#For descriptions, see 2.2

deg_Liver_24hRvs3hR <- results(deseq_object_Liver_vs3hR, alpha = 0.05, lfcThreshold = log2(1.2),  
                               contrast=c("Condition","D_Reoxygenation24h","C_Reoxygenation3h"))
deg_Liver_24hRvs3hR <- lfcShrink(deseq_object_Liver_vs3hR, coef="Condition_D_Reoxygenation24h_vs_C_Reoxygenation3h",
                                 res=deg_Liver_24hRvs3hR, type='apeglm')
deg_Liver_24hRvs3hR_df=as.data.frame(deg_Liver_24hRvs3hR)
deg_Liver_24hRvs3hR_df$gene_id<-as.vector(row.names(deg_Liver_24hRvs3hR_df))
deg_Liver_24hRvs3hR_df=deg_Liver_24hRvs3hR_df[, c(6, 1:5)]
row.names(deg_Liver_24hRvs3hR_df)=c()
deg_Liver_24hRvs3hR_df = merge(annotation, deg_Liver_24hRvs3hR_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Liver_24hRvs3hR_df <- deg_Liver_24hRvs3hR_df[ order( deg_Liver_24hRvs3hR_df$padj ), ]

#Export dataframe
write.table(deg_Liver_24hRvs3hR_df,"./output_DESeq/deg_Liver_24hRvs3hR.txt",quote=F,sep="\t",col.names = T,row.names = F)



#-----Step 3: Create DESeq2 Object for Gill-----

#Filter for rows with tissue = Gill (create index of row numbersw)
ind_Gill = which(coldata$Tissue == "Gill") 

#Keep the rows with the index numbers
coldata_Gill = coldata[ind_Gill,]

## Create DESeq object containing the gill tissue data 
#Gill_filt means there are at least 6 samples with non-zero counts for each gene
#coldata_liver means rows from sample info file with tissue = gill
deseq_object_Gill <- DESeqDataSetFromMatrix(countData = Gill_filt,
                                      colData = coldata_Gill, design =~ Condition)

#Normalization step (adjusts raw counts to make them comparable)
deseq_object_Gill <- estimateSizeFactors(deseq_object_Gill)



#-----Step 4: Run DESeq2 For Gill-----

#--4.1: Run DESeq using normoxia as refence level--

##"Make" a new DESeq object using the complete Gill DESeq object
deseq_object_Gill_vsN <- deseq_object_Gill

#Set reference level for Condition; first, we use Normoxia
deseq_object_Gill_vsN$Condition = relevel(deseq_object_Gill_vsN$Condition,"A_Normoxia") 

## Run DESeq function (comparing it to all other conditions)
deseq_object_Gill_vsN <- DESeq(deseq_object_Gill_vsN)

#Show comparisons from DESeq (what was compared against what)
resultsNames(deseq_object_Gill_vsN)



#--4.2: Working with output from DESeq (Anoxia vs Normoxia)--

#Filtering so that we only keep the results that show significant changes (p <= 0.05)
#and logfold change >= log2(1.2) (increase/decrease of ~26%)
#Also only compare anoxia vs normoxia
#deg = differentially expressed genes
deg_Gill_AvsN <- results(deseq_object_Gill_vsN, alpha = 0.05, lfcThreshold = log2(1.2), 
                         contrast=c("Condition","B_Anoxia","A_Normoxia"))

#Shrinkage is applied to reduce the variance of log fold changes for genes with low count numbers
deg_Gill_AvsN <- lfcShrink(deseq_object_Gill_vsN, coef="Condition_B_Anoxia_vs_A_Normoxia",
                           res=deg_Gill_AvsN, type='apeglm')

#Convert the DESeq results into a dataframe
deg_Gill_AvsN_df=as.data.frame(deg_Gill_AvsN)

#Create a column of gene IDs using the rownames
deg_Gill_AvsN_df$gene_id<-as.vector(row.names(deg_Gill_AvsN_df))

#Move gene ID column to first position
deg_Gill_AvsN_df=deg_Gill_AvsN_df[, c(6, 1:5)]

#Remove the rownames (index numbers instead)
row.names(deg_Gill_AvsN_df)=c()

#Merge with annotation file based on gene ID (to add protein product and symbols)
deg_Gill_AvsN_df = merge(annotation, deg_Gill_AvsN_df, by.x="gene_id", by.y="gene_id", all.y =T)

#Sort the dataframe based on padj (smallest first)
deg_Gill_AvsN_df <- deg_Gill_AvsN_df[ order( deg_Gill_AvsN_df$padj ), ]

#Export dataframe
write.table(deg_Gill_AvsN_df,"./output_DESeq/deg_Gill_AvsN.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--4.3: Working with output from DESeq (3hr reox vs Normoxia)--

#For descriptions, see 4.2

deg_Gill_3hRvsN <- results(deseq_object_Gill_vsN, alpha = 0.05, lfcThreshold = log2(1.2),  
                           contrast=c("Condition","C_Reoxygenation3h","A_Normoxia"))
deg_Gill_3hRvsN <- lfcShrink(deseq_object_Gill_vsN, coef="Condition_C_Reoxygenation3h_vs_A_Normoxia",
                             res=deg_Gill_3hRvsN, type='apeglm')
deg_Gill_3hRvsN_df=as.data.frame(deg_Gill_3hRvsN)
deg_Gill_3hRvsN_df$gene_id<-as.vector(row.names(deg_Gill_3hRvsN_df))
deg_Gill_3hRvsN_df=deg_Gill_3hRvsN_df[, c(6, 1:5)]
row.names(deg_Gill_3hRvsN_df)=c()
deg_Gill_3hRvsN_df = merge(annotation, deg_Gill_3hRvsN_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Gill_3hRvsN_df <- deg_Gill_3hRvsN_df[ order( deg_Gill_3hRvsN_df$padj ), ]

#Export dataframe
write.table(deg_Gill_3hRvsN_df,"./output_DESeq/deg_Gill_3hRvsN.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--4.4: Working with output from DESeq (24hr Reox vs Normoxia)--

#For descriptions, see 4.2

deg_Gill_24hRvsN <- results(deseq_object_Gill_vsN, alpha = 0.05, lfcThreshold = log2(1.2),  
                            contrast=c("Condition","D_Reoxygenation24h","A_Normoxia"))
deg_Gill_24hRvsN <- lfcShrink(deseq_object_Gill_vsN, coef="Condition_D_Reoxygenation24h_vs_A_Normoxia",
                              res=deg_Gill_24hRvsN, type='apeglm')
deg_Gill_24hRvsN_df=as.data.frame(deg_Gill_24hRvsN)
deg_Gill_24hRvsN_df$gene_id<-as.vector(row.names(deg_Gill_24hRvsN_df))
deg_Gill_24hRvsN_df=deg_Gill_24hRvsN_df[, c(6, 1:5)]
row.names(deg_Gill_24hRvsN_df)=c()
deg_Gill_24hRvsN_df = merge(annotation, deg_Gill_24hRvsN_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Gill_24hRvsN_df <- deg_Gill_24hRvsN_df[ order( deg_Gill_24hRvsN_df$padj ), ]

#Export dataframe
write.table(deg_Gill_24hRvsN_df,"./output_DESeq/deg_Gill_24hRvsN.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--4.5: Run DESeq using anoxia as refence level--

##"Make" a new DESeq object using the complete Gill DESeq object
deseq_object_Gill_vsA <- deseq_object_Gill

#Set reference level for Condition; now, we use Anoxia
deseq_object_Gill_vsA$Condition = relevel(deseq_object_Gill_vsA$Condition,"B_Anoxia")

## Run DESeq function (comparing it to all other conditions)
deseq_object_Gill_vsA <- DESeq(deseq_object_Gill_vsA)



#--4.6: Working with output from DESeq (3hr Reox vs Anoxia)--

#Filtering so that we only keep the results that show significant changes (p <= 0.05)
#and logfold change >= log2(1.2) (increase/decrease of ~26%)
#Also only compare 3hr Reox vs Anoxia
#deg = differentially expressed genes
deg_Gill_3hRvsA <- results(deseq_object_Gill_vsA, alpha = 0.05, lfcThreshold = log2(1.2),  
                           contrast=c("Condition","C_Reoxygenation3h","B_Anoxia"))
deg_Gill_3hRvsA <- lfcShrink(deseq_object_Gill_vsA, coef="Condition_C_Reoxygenation3h_vs_B_Anoxia",
                             res=deg_Gill_3hRvsA, type='apeglm')
deg_Gill_3hRvsA_df=as.data.frame(deg_Gill_3hRvsA)
deg_Gill_3hRvsA_df$gene_id<-as.vector(row.names(deg_Gill_3hRvsA_df))
deg_Gill_3hRvsA_df=deg_Gill_3hRvsA_df[, c(6, 1:5)]
row.names(deg_Gill_3hRvsA_df)=c()
deg_Gill_3hRvsA_df = merge(annotation, deg_Gill_3hRvsA_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Gill_3hRvsA_df <- deg_Gill_3hRvsA_df[ order( deg_Gill_3hRvsA_df$padj ), ]

#Export dataframe
write.table(deg_Gill_3hRvsA_df,"./output_DESeq/deg_Gill_3hRvsA.txt",quote=F,sep="\t",col.names = T,row.names = F)



#--4.7: Working with output from DESeq (24hr Reox vs Anoxia)--

#For descriptions, see 4.2

deg_Gill_24hRvsA <- results(deseq_object_Gill_vsA, alpha = 0.05, lfcThreshold = log2(1.2),  
                            contrast=c("Condition","D_Reoxygenation24h","B_Anoxia"))
deg_Gill_24hRvsA <- lfcShrink(deseq_object_Gill_vsA, coef="Condition_D_Reoxygenation24h_vs_B_Anoxia",
                              res=deg_Gill_24hRvsA, type='apeglm')
deg_Gill_24hRvsA_df=as.data.frame(deg_Gill_24hRvsA)
deg_Gill_24hRvsA_df$gene_id<-as.vector(row.names(deg_Gill_24hRvsA_df))
deg_Gill_24hRvsA_df=deg_Gill_24hRvsA_df[, c(6, 1:5)]
row.names(deg_Gill_24hRvsA_df)=c()
deg_Gill_24hRvsA_df = merge(annotation, deg_Gill_24hRvsA_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Gill_24hRvsA_df <- deg_Gill_24hRvsA_df[ order( deg_Gill_24hRvsA_df$padj ), ]

#Export dataframe
write.table(deg_Gill_24hRvsA_df,"./output_DESeq/deg_Gill_24hRvsA.txt",quote=F,sep="\t",col.names = T,row.names = F)


#--4.8: Run DESeq using 3hr Reox as refence level--

##"Make" a new DESeq object using the complete Gill DESeq object
deseq_object_Gill_vs3hR <- deseq_object_Gill

#Set reference level for Condition; now, we use 3hr Reox
deseq_object_Gill_vs3hR$Condition = relevel(deseq_object_Gill_vs3hR$Condition,"C_Reoxygenation3h")

## Run DESeq function (comparing it to all other conditions)
deseq_object_Gill_vs3hR <- DESeq(deseq_object_Gill_vs3hR)



#--4.9: Working with output from DESeq (24hr Reox vs 3hr Reox)--

#For descriptions, see 4.2

deg_Gill_24hRvs3hR <- results(deseq_object_Gill_vs3hR, alpha = 0.05, lfcThreshold = log2(1.2),  
                              contrast=c("Condition","D_Reoxygenation24h","C_Reoxygenation3h"))
deg_Gill_24hRvs3hR <- lfcShrink(deseq_object_Gill_vs3hR, coef="Condition_D_Reoxygenation24h_vs_C_Reoxygenation3h",
                                res=deg_Gill_24hRvs3hR, type='apeglm')
deg_Gill_24hRvs3hR_df=as.data.frame(deg_Gill_24hRvs3hR)
deg_Gill_24hRvs3hR_df$gene_id<-as.vector(row.names(deg_Gill_24hRvs3hR_df))
deg_Gill_24hRvs3hR_df=deg_Gill_24hRvs3hR_df[, c(6, 1:5)]
row.names(deg_Gill_24hRvs3hR_df)=c()
deg_Gill_24hRvs3hR_df = merge(annotation, deg_Gill_24hRvs3hR_df, by.x="gene_id", by.y="gene_id", all.y =T)
deg_Gill_24hRvs3hR_df <- deg_Gill_24hRvs3hR_df[ order( deg_Gill_24hRvs3hR_df$padj ), ]

#Export dataframe
write.table(deg_Gill_24hRvs3hR_df,"./output_DESeq/deg_Gill_24hRvs3hR.txt",quote=F,sep="\t",col.names = T,row.names = F)



#-----Step 5: PCA plots (separate)-----

#Variance stabilisation for PCA plots:
#normalize the counts in a way that makes them more suitable for downstream analysis, such as visualization or clustering
vsd_Liver <- vst(deseq_object_Liver)
vsd_Gill <- vst(deseq_object_Gill)

#Make color df based on treatment
color_info <- data.frame(
  treatment = c("Normoxia", "Anoxia", "3h Reox", "24h Reox"), 
  color = c("darkorange", "blue", "deeppink", "darkorchid")
)


#--Liver PCA plot--
#generate PCA plot, color/label by "condition"
#Return a dataframe with PCA results
#Use the top 5000 most variable genes
pcaData_Liver <- plotPCA(vsd_Liver, intgroup="Condition", returnData=TRUE, ntop = 5000)

#calculate the percentage of variance explained by the principal components derived from the PCA analysis
percentVar_Liver <- round(100 * attr(pcaData_Liver, "percentVar"))

#Create df containing sex of sample
sex_info_Liver <- data.frame(
  sample_name = c("L.N1","L.N2","L.N3","L.N4","L.N5","L.N6",
                  "L.A1","L.A2","L.A3","L.A4","L.A5","L.A6",
                  "L.R11","L.R12","L.R13","L.R14","L.R15","L.R16",
                  "L.R21","L.R22","L.R23","L.R24","L.R25","L.R26"), 
  Sex = c("F","M","M","F","M","F",
          "M","F","F","M","M","M",
          "F","F","F","F","F","F",
          "F","F","F","F","F","M")
)

#Add the sex of each sample to PCA data
pcaData_Liver <- merge(pcaData_Liver, sex_info_Liver, by.x = "name", by.y = "sample_name")

# Visualize PCA plot using GGplot
pca12Liverplot = ggplot(pcaData_Liver, aes(PC1, PC2, group = Condition)) +
  geom_point(aes(shape = Sex, colour = Condition, size = Condition)) +
  ggtitle("Liver") + 
  xlab(paste0("PC1: ", percentVar_Liver[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_Liver[2], "% variance")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(aspect.ratio = 1) +
  theme(legend.background = element_blank(), 
        legend.key = element_blank(), 
        axis.text = element_text(size = rel(1.5), colour = "black"),
        axis.title = element_text(size = rel(1.5)), 
        legend.title = element_text(size = rel(1.4)), 
        legend.text = element_text(size = rel(1.4)), 
        legend.position = c(1.2, 0.78)) +
  scale_shape_manual(name = "Sex", values = c(16, 17),
                     labels = c("F" = "Female", "M" = "Male"),
                     guide = guide_legend(override.aes = list(size = 3.5))) +
  scale_colour_manual(name = "Condition", 
                      values = color_info$color, 
                      labels = color_info$treatment) +
  scale_size_manual(values = c(3.5, 3.5, 3.5, 3.5), 
                    labels = c("Normoxia", "Anoxia", "3h Reox", "24h Reox"))

# Print the plot
print(pca12Liverplot)

#Save plot separately
cairo_pdf(file=paste0("./output_DESeq/pca_Liver.pdf"), width = 8, height = 5)
print(pca12Liverplot)
dev.off()

tiff(file=paste0("./output_DESeq/pca_Liver.tiff"), width = 19, height = 12, units = "cm", res = 300)
print(pca12Liverplot)
dev.off()


#--Gill PCA plot--
#generate PCA plot, color/label by "condition"
#Return a dataframe with PCA results
#Use the top 5000 most variable genes
pcaData_Gill <- plotPCA(vsd_Gill, intgroup="Condition", returnData=TRUE, ntop = 5000)

#calculate the percentage of variance explained by the principal components derived from the PCA analysis
percentVar_Gill <- round(100 * attr(pcaData_Gill, "percentVar"))

#Create df containing sex of sample
sex_info_Gill <- data.frame(
  sample_name = c("G.N2","G.N3","G.N4","G.N6","G.N7","G.N8",
                  "G.A1","G.A2","G.A3","G.A5","G.A7","G.A8",
                  "G.R11","G.R12","G.R13","G.R14","G.R17","G.R18",
                  "G.R22","G.R23","G.R24","G.R25","G.R27","G.R28"), 
  Sex = c("F","F","M","F","M","M",
          "M","F","F","M","M","M",
          "F","F","F","F","M","F",
          "M","M","M","F","M","F")
)

#Add the sex of each sample to PCA data
pcaData_Gill <- merge(pcaData_Gill, sex_info_Gill, by.x = "name", by.y = "sample_name")

# Visualize PCA plot using GGplot
pca12Gillplot = ggplot(pcaData_Gill, aes(PC1, PC2, group = Condition)) +
  geom_point(aes(shape = Sex, colour = Condition, size = Condition)) +
  ggtitle("Gill") + 
  xlab(paste0("PC1: ", percentVar_Gill[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_Gill[2], "% variance")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(aspect.ratio = 1) +
  theme(legend.background = element_blank(), 
        legend.key = element_blank(), 
        axis.text = element_text(size = rel(1.5), colour = "black"),
        axis.title = element_text(size = rel(1.5)), 
        legend.title = element_text(size = rel(1.4)), 
        legend.text = element_text(size = rel(1.4)), 
        legend.position = c(1.2, 0.78)) +
  scale_shape_manual(name = "Sex", values = c(16, 17),
                     labels = c("F" = "Female", "M" = "Male"),
                     guide = guide_legend(override.aes = list(size = 3.5))) +
  scale_colour_manual(name = "Condition", 
                      values = color_info$color, 
                      labels = color_info$treatment) +
  scale_size_manual(values = c(3.5, 3.5, 3.5, 3.5), 
                    labels = c("Normoxia", "Anoxia", "3h Reox", "24h Reox"))

# Print the plot
print(pca12Gillplot)

#Save plot separately
cairo_pdf(file=paste0("./output_DESeq/pca_Gill.pdf"), width = 8, height = 5)
print(pca12Gillplot)
dev.off()

tiff(file=paste0("./output_DESeq/pca_Gill.tiff"), width = 19, height = 12, units = "cm", res = 300)
print(pca12Gillplot)
dev.off()


#Combine and save the PCA plots!
cairo_pdf(file=paste0("./output_DESeq/pcaplots.pdf"), width = 14, height = 5)
print(grid.arrange(pca12Gillplot,pca12Liverplot, ncol=2, nrow=1))
dev.off()

tiff(file=paste0("./output_DESeq/pcaplots.tiff"), width = 36, height = 12, units = "cm", res = 300)
print(grid.arrange(pca12Gillplot,pca12Liverplot, ncol=2, nrow=1))
dev.off()

