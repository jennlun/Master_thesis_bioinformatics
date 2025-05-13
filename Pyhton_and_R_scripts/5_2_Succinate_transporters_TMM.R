## Set Working directory
setwd("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/5_Succinate_transporters")

##Load Packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)


##--Upload data files--
#Upload file with GO ids and gene ids
GO_gene_id_file <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/3_GO_analysis/GO_source_files/ccar2swissprotGO_long.txt", 
                              header=T, stringsAsFactors=FALSE, sep="\t")

#upload list of succinate transporter GO ids and GO terms
succinate_go_ids = read.delim("./Succinate_input_files/2025_01_17_go_terms_succinate_transport.txt", header=T, stringsAsFactors=FALSE, sep="\t")

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

#Upload tmm files
Liver_tmm <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/4_plot_venn_heat/output_files/Liver_genes_nonfiltered_TMM.matrix")
Gill_tmm <- read.delim("C:/Users/jenny/OneDrive - Universitetet i Oslo/Master/2025_05_03_RNAseq_3hr/4_plot_venn_heat/output_files/Gill_genes_nonfiltered_TMM.matrix")


#-----Step 1: Create a dataframe with only the relevant geneIDs-----

#Only keep the geneIDs corresponding to the relevant GO terms
relevant_genes <- GO_gene_id_file[GO_gene_id_file$GO_id %in% succinate_go_ids$go_id, ]

#Add the corresponding GO terms
relevant_genes <- merge(relevant_genes, succinate_go_ids, by.x = "GO_id", by.y = "go_id")



#-----Step 2: Go through the raw DEG files and pick out relevant genes-----

#--2.1: Add comparison label to DEG dfs--
#Add new column containing the comparison name on each row
#Liver files
deg_Liver_AvsN$comparison = "AvsN"
deg_Liver_3hRvsN$comparison = "3hRvsN"
deg_Liver_24hRvsN$comparison = "24hRvsN"
deg_Liver_3hRvsA$comparison = "3hRvsA"
deg_Liver_24hRvsA$comparison = "24hRvsA"
deg_Liver_24hRvs3hR$comparison = "24hRvs3hR"

#Gill Files
deg_Gill_AvsN$comparison = "AvsN"
deg_Gill_3hRvsN$comparison = "3hRvsN"
deg_Gill_24hRvsN$comparison = "24hRvsN"
deg_Gill_3hRvsA$comparison = "3hRvsA"
deg_Gill_24hRvsA$comparison = "24hRvsA"
deg_Gill_24hRvs3hR$comparison = "24hRvs3hR"



#--2.2: Filter out relevant genes from each comparison--
#Keep only rows that have the relevant gene IDs in each DEG file
#Liver data
relevant_deg_Liver_AvsN= deg_Liver_AvsN[deg_Liver_AvsN$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Liver_3hRvsN= deg_Liver_3hRvsN[deg_Liver_3hRvsN$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Liver_24hRvsN= deg_Liver_24hRvsN[deg_Liver_24hRvsN$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Liver_3hRvsA= deg_Liver_3hRvsA[deg_Liver_3hRvsA$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Liver_24hRvsA= deg_Liver_24hRvsA[deg_Liver_24hRvsA$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Liver_24hRvs3hR= deg_Liver_24hRvs3hR[deg_Liver_24hRvs3hR$gene_id %in% relevant_genes$gene_id, ]

#Gill data
relevant_deg_Gill_AvsN= deg_Gill_AvsN[deg_Gill_AvsN$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Gill_3hRvsN= deg_Gill_3hRvsN[deg_Gill_3hRvsN$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Gill_24hRvsN= deg_Gill_24hRvsN[deg_Gill_24hRvsN$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Gill_3hRvsA= deg_Gill_3hRvsA[deg_Gill_3hRvsA$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Gill_24hRvsA= deg_Gill_24hRvsA[deg_Gill_24hRvsA$gene_id %in% relevant_genes$gene_id, ]
relevant_deg_Gill_24hRvs3hR= deg_Gill_24hRvs3hR[deg_Gill_24hRvs3hR$gene_id %in% relevant_genes$gene_id, ]

#Combine all the relevant DEG data into one dataframe
relevant_Liver_degs <- rbind(relevant_deg_Liver_AvsN, relevant_deg_Liver_3hRvsN, relevant_deg_Liver_24hRvsN, 
                             relevant_deg_Liver_3hRvsA, relevant_deg_Liver_24hRvsA, relevant_deg_Liver_24hRvs3hR)

relevant_Gill_degs <- rbind(relevant_deg_Gill_AvsN, relevant_deg_Gill_3hRvsN, relevant_deg_Gill_24hRvsN, 
                            relevant_deg_Gill_3hRvsA, relevant_deg_Gill_24hRvsA, relevant_deg_Gill_24hRvs3hR)

#Sort the dataframes based on the gene IDs
relevant_Liver_degs <- relevant_Liver_degs[order(relevant_Liver_degs$gene_id), ]
relevant_Gill_degs <- relevant_Gill_degs[order(relevant_Gill_degs$gene_id), ]

#Only keep the Tmms corresponding to the relevant gene IDs
relevant_Liver_tmm <- Liver_tmm[ Liver_tmm$gene_id %in% relevant_Liver_degs$gene_id, ]
relevant_Gill_tmm <- Gill_tmm[ Gill_tmm$gene_id %in% relevant_Gill_degs$gene_id, ]



#-----Step 3: Calculate mean of treatment groups and add to original dataframe-----

#Start with Liver
#make df with only normoxic values
r_Liver_tmm_N = relevant_Liver_tmm %>% 
  dplyr::select(starts_with("L.N"))

#Calculate mean for each gene
r_Liver_tmm_N_means <- rowMeans(r_Liver_tmm_N)

#Add the means as a new column in the original df
relevant_Liver_tmm <- relevant_Liver_tmm %>%
  mutate(Normoxia_mean = r_Liver_tmm_N_means)


#Repeat for Anoxia
r_Liver_tmm_A = relevant_Liver_tmm %>% 
  dplyr::select(starts_with("L.A"))

r_Liver_tmm_A_means <- rowMeans(r_Liver_tmm_A)

#create new dataframe with all means
all_Liver_means <- relevant_Liver_tmm %>%
  mutate(Anoxia_mean = r_Liver_tmm_A_means)

#keep only gene_ID and mean
all_Liver_means = all_Liver_means %>% dplyr::select(gene_id, Normoxia_mean, Anoxia_mean)


#Repeat for 3h Reox
r_Liver_tmm_3h = relevant_Liver_tmm %>% 
  dplyr::select(starts_with("L.R1"))

r_Liver_tmm_3h_means <- rowMeans(r_Liver_tmm_3h)

all_Liver_means <- all_Liver_means %>%
  mutate("3hReox_mean" = r_Liver_tmm_3h_means)


#Repeat for 24h Reox
r_Liver_tmm_24h = relevant_Liver_tmm %>% 
  dplyr::select(starts_with("L.R2"))

r_Liver_tmm_24h_means <- rowMeans(r_Liver_tmm_24h)

all_Liver_means <- all_Liver_means %>%
  mutate("24hReox_mean" = r_Liver_tmm_24h_means)


#Repeat for Gill
r_Gill_tmm_N = relevant_Gill_tmm %>%
  dplyr::select(starts_with("G.N"))

r_Gill_tmm_N_means <- rowMeans(r_Gill_tmm_N)

relevant_Gill_tmm <- relevant_Gill_tmm %>%
  mutate(Normoxia_mean = r_Gill_tmm_N_means)


#Repeat for Anoxia
r_Gill_tmm_A = relevant_Gill_tmm %>% 
  dplyr::select(starts_with("G.A"))

r_Gill_tmm_A_means <- rowMeans(r_Gill_tmm_A)

#create new dataframe with all means
all_Gill_means <- relevant_Gill_tmm %>%
  mutate(Anoxia_mean = r_Gill_tmm_A_means)

#keep only gene_ID and mean
all_Gill_means = all_Gill_means %>% dplyr::select(gene_id, Normoxia_mean, Anoxia_mean)


#Repeat for 3h Reox
r_Gill_tmm_3h = relevant_Gill_tmm %>% 
  dplyr::select(starts_with("G.R1"))

r_Gill_tmm_3h_means <- rowMeans(r_Gill_tmm_3h)

all_Gill_means <- all_Gill_means %>%
  mutate("3hReox_mean" = r_Gill_tmm_3h_means)


#Repeat for 24h Reox
r_Gill_tmm_24h = relevant_Gill_tmm %>% 
  dplyr::select(starts_with("G.R2"))

r_Gill_tmm_24h_means <- rowMeans(r_Gill_tmm_24h)

all_Gill_means <- all_Gill_means %>%
  mutate("24hReox_mean" = r_Gill_tmm_24h_means)

#Remove "symbol" and "product" columns
relevant_Liver_tmm = relevant_Liver_tmm %>% dplyr::select(-symbol, -product)
relevant_Gill_tmm = relevant_Gill_tmm %>% dplyr::select(-symbol, -product)



#-----Step 4: Modify dataframes-----

###Liver samples 
#Flip df so that the columns are rows
relevant_Liver_tmm_treat <- as.data.frame(t(relevant_Liver_tmm))

#Make colnames gene IDs
colnames(relevant_Liver_tmm_treat) = relevant_Liver_tmm_treat[1,]

#remove the first and last row
relevant_Liver_tmm_treat <- relevant_Liver_tmm_treat[-1,] 
relevant_Liver_tmm_treat <- relevant_Liver_tmm_treat[-nrow(relevant_Liver_tmm_treat), ]

#Add column containing the treatment for each sample
relevant_Liver_tmm_treat$treatment <- rep(NA, nrow(relevant_Liver_tmm_treat))

#turn row names into a list
a_L_row <- row.names(relevant_Liver_tmm_treat)

#Go through the row names and assign the correct treatment to the "treatment" column
for (i in seq_along(a_L_row)) {
  row_name <- a_L_row[i]
  if (grepl("^L.N", row_name)) {
    relevant_Liver_tmm_treat$treatment[i] <- "Normoxia"
  } else if (grepl("^L.A", row_name)) {
    relevant_Liver_tmm_treat$treatment[i] <- "Anoxia"
  } else if (grepl("^L.R1", row_name)) {
    relevant_Liver_tmm_treat$treatment[i] <- "3h Reox"
  } else if (grepl("^L.R2", row_name)) {
    relevant_Liver_tmm_treat$treatment[i] <- "24h Reox"
  }
}


###Gill samples
#Flip df so that the columns are rows
relevant_Gill_tmm_treat <- as.data.frame(t(relevant_Gill_tmm))

#Make colnames gene IDs
colnames(relevant_Gill_tmm_treat) = relevant_Gill_tmm_treat[1,]

#remove the first and last row
relevant_Gill_tmm_treat <- relevant_Gill_tmm_treat[-1,] 
relevant_Gill_tmm_treat <- relevant_Gill_tmm_treat[-nrow(relevant_Gill_tmm_treat), ]

#Add column containing the treatment for each sample
relevant_Gill_tmm_treat$treatment <- rep(NA, nrow(relevant_Gill_tmm_treat))

#turn row names into a list
a_G_row <- row.names(relevant_Gill_tmm_treat)

#Go through the row names and assign the correct treatment to the "treatment" column
for (i in seq_along(a_G_row)) {
  row_name <- a_G_row[i]
  if (grepl("^G.N", row_name)) {
    relevant_Gill_tmm_treat$treatment[i] <- "Normoxia"
  } else if (grepl("^G.A", row_name)) {
    relevant_Gill_tmm_treat$treatment[i] <- "Anoxia"
  } else if (grepl("^G.R1", row_name)) {
    relevant_Gill_tmm_treat$treatment[i] <- "3h Reox"
  } else if (grepl("^G.R2", row_name)) {
    relevant_Gill_tmm_treat$treatment[i] <- "24h Reox"
  }
}



#-----Step 5: Make Bar plots!-----

#--5.1: Prepare data--

#Flip df so that the columns are rows
all_Liver_means <- as.data.frame(t(all_Liver_means))

#Make colnames gene IDs
colnames(all_Liver_means) = all_Liver_means[1,]

#remove the first row
all_Liver_means <- all_Liver_means[-1,] 

#Add column with treatments
all_Liver_means$treatment = c("Normoxia", "Anoxia", "3h Reox", "24h Reox")

# Convert treatments to a factor 
all_Liver_means$treatment <- as.factor(all_Liver_means$treatment)

# Reverse treatment order
all_Liver_means$treatment <- factor(all_Liver_means$treatment, 
                                    levels = rev(levels(all_Liver_means$treatment)))


#Gill

#Flip df so that the columns are rows
all_Gill_means <- as.data.frame(t(all_Gill_means))

#Make colnames gene IDs
colnames(all_Gill_means) = all_Gill_means[1,]

#remove the first row
all_Gill_means <- all_Gill_means[-1,] 

#Add column with treatments
all_Gill_means$treatment = c("Normoxia", "Anoxia", "3h Reox", "24h Reox")

# Convert treatments to a factor 
all_Gill_means$treatment <- as.factor(all_Gill_means$treatment)

# Reverse treatment order
all_Gill_means$treatment <- factor(all_Gill_means$treatment, 
                                   levels = rev(levels(all_Gill_means$treatment)))



#--5.2: Make barplots for genes with relevant changes in gene expression--

#Make function for bar plotting!
#gene = gene name as string
#tissue = tissue name as string
bar_plot_relevant = function(all_means, rel_dev_treat, gene, tissue){
  
  # Ensure numeric conversion
  all_means[[gene]] <- as.numeric(as.character(all_means[[gene]]))
  rel_dev_treat[[gene]] <- as.numeric(as.character(rel_dev_treat[[gene]]))
  
  # Create treatment information DataFrame
  treatment_info <- data.frame(
    treatment = c("Normoxia", "Anoxia", "3h Reox", "24h Reox"), 
    color1 = c("darkorange", "blue", "deeppink", "darkorchid"),
    color2 = c("darkorange3", "blue4", "deeppink3", "darkorchid4")
  )
  
  # Calculate min and max of relative deviation for each treatment
  max_rel_dev <- rel_dev_treat %>%
    group_by(treatment) %>% 
    summarize(max_dev = max(!!sym(gene), na.rm = TRUE))
  
  min_rel_dev <- rel_dev_treat %>%
    group_by(treatment) %>% 
    summarize(min_dev = min(!!sym(gene), na.rm = TRUE))
  
  # Add min and max to treatment info dataframe
  treatment_info <- merge(treatment_info, max_rel_dev, by = "treatment")
  treatment_info <- merge(treatment_info, min_rel_dev, by = "treatment")
  
  # Calculate max value for scaling
  max_value <- max(rel_dev_treat[[gene]], na.rm = TRUE) * 1.2
  y_breaks <- seq(0, max_value, by = (max_value/5))
  
  # Plot
  bar_plot_gene = ggplot(all_means, aes(x = treatment, y = !!sym(gene))) +  
    geom_col(aes(fill = treatment)) + 
    labs(fill = "Condition", title = gene, subtitle = tissue) +
    scale_y_continuous(limits = c(0, max_value), expand = c(0,0), breaks = y_breaks,
                       labels = number_format(accuracy = 0.01)) +
    scale_fill_manual(values = setNames(treatment_info$color1, treatment_info$treatment)) +
    geom_point(data = rel_dev_treat, aes(color = treatment), 
               position = position_jitter(width = 0.0025, height = 0), show.legend = FALSE) +
    scale_color_manual(values = setNames(treatment_info$color2, treatment_info$treatment)) + 
    ylab("Counts per Million") +
    xlab("Condition") +theme(
      panel.background = element_rect(fill = "white"))   # White background for the panel
  
  print(bar_plot_gene)
  
  # Save the plot as PDF
  
  file_path = "./TMM_plots/"
  
  # Save the plot as TIFF
  tiff(file = paste0(file_path, gene, "_bar_plot_", tissue, ".tiff"), width = 4, height = 3.2, units = "in", res = 300) 
  print(bar_plot_gene)
  dev.off()
}

#Run function on Liver genes
bar_plot_relevant(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua21-g19814", "Liver")
bar_plot_relevant(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua21-g20040", "Liver")
bar_plot_relevant(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ub18-g39864", "Liver")
bar_plot_relevant(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ub21-g41947", "Liver")

#Run function on Gill genes
bar_plot_relevant(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ua04-g4151", "Gill")
bar_plot_relevant(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ub18-g39831", "Gill")
bar_plot_relevant(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ub18-g39864", "Gill")



#--5.2: Make barplots for remaining genes--

#Make function for bar plotting!
#gene = gene name as string
#tissue = tissue name as string
bar_plot = function(all_means, rel_dev_treat, gene, tissue){
  
  # Ensure numeric conversion
  all_means[[gene]] <- as.numeric(as.character(all_means[[gene]]))
  rel_dev_treat[[gene]] <- as.numeric(as.character(rel_dev_treat[[gene]]))
  
  # Create treatment information DataFrame
  treatment_info <- data.frame(
    treatment = c("Normoxia", "Anoxia", "3h Reox", "24h Reox"), 
    color1 = c("darkorange", "blue", "deeppink", "darkorchid"),
    color2 = c("darkorange3", "blue4", "deeppink3", "darkorchid4")
  )
  
  # Calculate min and max of relative deviation for each treatment
  max_rel_dev <- rel_dev_treat %>%
    group_by(treatment) %>% 
    summarize(max_dev = max(!!sym(gene), na.rm = TRUE))
  
  min_rel_dev <- rel_dev_treat %>%
    group_by(treatment) %>% 
    summarize(min_dev = min(!!sym(gene), na.rm = TRUE))
  
  # Add min and max to treatment info dataframe
  treatment_info <- merge(treatment_info, max_rel_dev, by = "treatment")
  treatment_info <- merge(treatment_info, min_rel_dev, by = "treatment")
  
  # Calculate max value for scaling
  max_value <- max(rel_dev_treat[[gene]], na.rm = TRUE) * 1.2
  y_breaks <- seq(0, max_value, by = (max_value/5))
  
  # Plot
  bar_plot_gene = ggplot(all_means, aes(x = treatment, y = !!sym(gene))) +  
    geom_col(aes(fill = treatment)) + 
    labs(fill = "Condition", title = gene, subtitle = tissue) +
    scale_y_continuous(limits = c(0, max_value), expand = c(0,0), breaks = y_breaks,
                       labels = number_format(accuracy = 0.01)) +
    scale_fill_manual(values = setNames(treatment_info$color1, treatment_info$treatment)) +
    geom_point(data = rel_dev_treat, aes(color = treatment), 
               position = position_jitter(width = 0.0025, height = 0), show.legend = FALSE) +
    scale_color_manual(values = setNames(treatment_info$color2, treatment_info$treatment)) + 
    ylab("Counts per Million") +
    xlab("Condition") +theme(
      panel.background = element_rect(fill = "white"))   # White background for the panel
  
  print(bar_plot_gene)
  
  # Save the plot as PDF
  
  file_path = "./Non_relevant_TMM_plots/"
  
  # Save the plot as TIFF
  tiff(file = paste0(file_path, gene, "_bar_plot_", tissue, ".tiff"), width = 4, height = 3.2, units = "in", res = 300) 
  print(bar_plot_gene)
  dev.off()
}

#Run function on Liver genes
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua01-g1240", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua04-g4151", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua18-g17104", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua21-g19814", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ua21-g19855", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ub04-g27034", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ub18-g39831", "Liver")
bar_plot(all_Liver_means, relevant_Liver_tmm_treat, "ccar_ub21-g41947", "Liver")

#Run function on Gill genes
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ua01-g1240", "Gill")
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ua18-g17104", "Gill")
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ua21-g19814", "Gill")
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ua21-g19855", "Gill")
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ua21-g20040", "Gill")
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ub04-g27034", "Gill")
bar_plot(all_Gill_means, relevant_Gill_tmm_treat, "ccar_ub21-g41947", "Gill")






