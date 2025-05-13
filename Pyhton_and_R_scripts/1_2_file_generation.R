## Set Working directory
setwd("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/1_file_generation")

#load packages
library(ontologyIndex)


## Upload crucian carp transcript ID and protein ID from swissprot
ccar2swissprot=read.delim("./source_files/blastp.swissprotGO.e01bit30.outfmt6", header=F, stringsAsFactors=FALSE, sep="\t")
ccar2swissprot=ccar2swissprot[,c(1,2)] #switch column order
colnames(ccar2swissprot)=c("ccar_transcript_id","protein_id") #change column name


## Upload transcript ID and gene ID
ccar_gene2trans=read.delim("./source_files/carcar_annotation_v5.trans_gene_map.txt", header=F, stringsAsFactors=FALSE, sep="\t")
colnames(ccar_gene2trans)=c("transcript_id","gene_id") #change column name
ccar_gene2trans <- ccar_gene2trans[, c(2,1)] #change column order


## Combine based on transcript ID to make annotation file
annotation = merge(ccar_gene2trans, ccar2swissprot, by.x = "transcript_id", by.y="ccar_transcript_id", all.x=T)


## Upload protein names from swissprot, and add to annotation file based on transcript ID
prot2names=read.delim("./source_files/uniprotkb_reviewed_true_2024_07_08_wGO.tsv", header=T, stringsAsFactors=FALSE, sep="\t") # protein names and descriptions for swissprot protein entries
prot2names <- prot2names[, c(1:3)] #keep only the first 3 columns
colnames(prot2names)=c("protein_ID", "symbol", "product") #change column names


## Add protein names to annotation file based on the protein ID
annotation = merge(annotation, prot2names, by.x = "protein_id", by.y = "protein_ID") #There are now 5 columns with info


## Save annotation table to output files
annotation_transcript <- annotation[, c(3,2,1,4,5)] #change column order
write.table(annotation_transcript, "./output_file_generation/blastp_hits_swissprot_transcripts.txt", quote=F,sep="\t",col.names = T, row.names = F)


## Remove the transcript ID and protein ID columns and save to output files (???)
annotation=annotation[!duplicated(annotation$gene_id),c(3:5)]
write.table(annotation, "./output_file_generation/blastp_hits_swissprot_genes.txt", quote=F,sep="\t",col.names = T, row.names = F)


## Upload and format the "full" GO file
swissprot_go <- read.delim("./source_files/ccar2swissprotGO.txt", header=T, stringsAsFactors=FALSE, sep="\t") # table that contains GO terms for swissprot protein entries
colnames(swissprot_go) <- c("gene_id", "GO_id")


## Save reformatted GO file
write.table(swissprot_go, "./output_file_generation/swissprot_go.txt", quote=F,sep="\t",col.names = T, row.names = F)


## Upload gene ontology information
obo <- get_OBO("./source_files/go-basic_release20240617.obo", extract_tags = "everything") # read in gene ontology information
go_data <- as.data.frame(cbind(as.character(obo$id), as.character(obo$name), as.character(obo$namespace))) #turn the OBO list into a data frame
colnames(go_data) <- c("id", "name", "namespace") #change column names
rownames(go_data) <- NULL #remove any rownames if present


## Shorten names
go_data$namespace <- gsub("biological_process", "BP", go_data$namespace) #shorten biological process to BP
go_data$namespace <- gsub("molecular_function", "MF", go_data$namespace) #shorten molecular function to MF
go_data$namespace <- gsub("cellular_component", "CC", go_data$namespace) #shorten cellular component to CC


## Turn the "namespace" column into a factor (like "categories")
go_data$namespace <- as.factor(go_data$namespace)


## Remove any rows where the namespace is "external"
go_data <- droplevels(go_data[!go_data$namespace=="external",])


## Save the GO data
write.table(go_data, "./output_file_generation/go_data.txt", quote=F,sep="\t",col.names = T, row.names = F)
