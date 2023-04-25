"""
# Solution to Case Study about Differential Expression of a Cancer Dataset

For this case study, one had to suppose to be a scientist investigating which genes and pathways are differentially expressed in cancer.
The gene expression of 3 affected patients (disease 1 - 3) and 3 controls (controls 1 - 3) of 1000 genes were measured. 
The resulting count matrix is stored in the file “expression_counts.txt”.

Differential expression using the library DESeq2 was called between two sets of RNA-seq samples (control and disease).
A list of differential expressed genes at False Discovery Rate (i.e. padj) after 0.05 was obtained, being then defined which 
ones of them are up and which ones are down regulated (part 1 of the script).

Additionally, the list of Tumour Suppressor Genes (TSGs) and oncogenes from the COSMIC database (https://cancer.sanger.ac.uk/census) 
was used to unravel if TSGs and/or oncogenes over-represented among over- or under-expressed genes, 
while assessing statistical significance of the enrichment, if any (part 2 of the script).

Author: Raquel Romão
"""

# load the necessary libraries
library(DESeq2)
library(plotly)
library(ggplot2)
library(purrr)
library(ggrepel)
library(clusterProfiler)
library(cowplot)

####1####

# load the expression count matrix
expression <- read.table('data/expression_counts.txt', header=TRUE, row.names=1)

# transpose expression count matrix
expression <- t(expression)

# add a condition column in the expression dataframe
conditions <- factor(c('disease', 'disease', 'disease', 'control', 'control','control'))

# Create a DataFrame containing the sample names and condition labels
sample_info <- data.frame(Sample = colnames(expression), Condition = conditions)
sample_info


# create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = expression, colData=sample_info, design = ~ Condition)

# fit the model
# note: relevel since DSEQ2 orders by alphabetical order (in this case its not an issue)
dds$Condition <- factor(dds$Condition, levels = c("control","disease")) 
dds$Condition <- relevel(dds$Condition, ref = "control")

dds <- DESeq(dds)

# extract the results
res <- results(dds)
write.table(as.data.frame(res), 'Results_table.tsv', sep='\t', row.names = T, col.names = TRUE, quote = FALSE)

# filter the results to include only genes with padj <= 0.05
res_filtered <- res[res$padj <= 0.05, ]
head(res_filtered)

# determine which genes are upregulated and which are downregulated
upregulated_genes <- res_filtered[res_filtered$log2FoldChange > 1, ]
downregulated_genes <- res_filtered[res_filtered$log2FoldChange < -1, ]


# convert res_filtered to a data frame
res_filtered_df <- as.data.frame(res_filtered)
head(res_filtered_df)

# create a MA plot
ma_plot <- ggplot(data = res_filtered_df, aes(x = baseMean, y = log2FoldChange, color = log2FoldChange > 1)) +
  geom_point() +
  scale_x_log10()

# specify different colors for the up- and down-regulated groups
ma_plot + scale_color_manual(values = c("red", "blue"))


# create an interactive version of the MA plot
ggplotly(ma_plot)

# show the MA plot
ma_plot

#############################################################################################################
# VolcanoPlot

## res as a dataframe and row names as column
res_df <- data.frame(res) %>%
  tibble::rownames_to_column(., "GeneSymbol")

## 5 most up and down-regulated genes (padj < 0.05 AND abs(log2FC) > 1):
data <- res_df %>%
  dplyr::filter(abs(log2FoldChange) > 1) %>%
  dplyr::filter(padj < 0.05)

### up-regulated
up <- data %>%
  dplyr::filter(log2FoldChange > 0) %>%
  dplyr::arrange(desc(padj)) %>%
  dplyr::slice_tail(n = 5)

### down-regulated
down <- data %>%
  dplyr::filter(log2FoldChange < 0) %>%
  dplyr::arrange(desc(padj)) %>%
  dplyr::slice_tail(n = 5)

### bing down and up regulated genes in a single dataframe
top <- base::rbind(up, down)

## Preparation for the Volcano Plot
res_df$Sig <- "NS"
res_df$Sig[(res_df$padj < 0.05)] <- "AdjPval"
res_df$Sig[(res_df$log2FoldChange < -1) & (res_df$padj < 0.05)] <- "log2FC_down_AdjPval"
res_df$Sig[(res_df$log2FoldChange > 1) & (res_df$padj < 0.05)] <- "log2FC_up_AdjPval"
res_df$Sig <- factor(res_df$Sig, levels=c("NS", "AdjPval", "log2FC_down_AdjPval", "log2FC_up_AdjPval"))

## Resume table of Sig column
table(res_df$Sig)

## set size for text
size = 16

## do volcano plot
volcano_plot <- ggplot2::ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  #Add points:
  #	Colour based on factors set a few lines up
  #	'alpha' provides gradual shading of colour
  #	Set size of points
  ggplot2::geom_point(aes(color=factor(Sig)), alpha=1/2, size=2) +
  ## show in plot 5 most up and down-regulated genes:
  ggrepel::geom_text_repel(data= top, aes(label=GeneSymbol)) +
  
  
  #Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
  ggplot2::scale_color_manual(values=c(NS="black", AdjPval ="grey30", log2FC_down_AdjPval="royalblue", log2FC_up_AdjPval="red2"),
                              labels=c(NS="AdjPval > 0.05", AdjPval = "AdjPval < 0.05", log2FC_down_AdjPval="down-regulated in disease vs control",
                                       log2FC_up_AdjPval="up-regulated in disease vs control")) +
  #Set the size of the plotting window
  ggplot2::theme_bw(base_size=24) +
  #Modify various aspects of the plot text and legend
  ggplot2::theme(legend.background=element_rect(),
                 plot.title=element_text(angle=0, size=size, face="bold", vjust=1),
                 panel.grid.major=element_blank(),	#Remove gridlines
                 panel.grid.minor=element_blank(),	#Remove gridlines
                 axis.text.x=element_text(angle=0, size=size, vjust=1),
                 axis.text.y=element_text(angle=0, size=size, vjust=1),
                 axis.title=element_text(size=size),
                 
                 #Legend
                 legend.position="top",			#Moves the legend to the top of the plot
                 legend.key=element_blank(),		#removes the border
                 legend.key.size=unit(0.5, "cm"),	#Sets overall area/size of the legend
                 legend.text=element_text(size=size),	#Text size
                 title=element_text(size=size),		#Title text size
                 legend.title=element_blank()) +		#Remove the title
  
  #Change the size of the icons/symbols in the legend
  ggplot2::guides(colour = guide_legend(override.aes=list(size=2.5))) +
  #Set the axis limits
  ggplot2::xlim(-abs(max(res_df$log2FoldChange)), abs(max(res_df$log2FoldChange))) +
  ggplot2::xlab("log2 Fold Change in disease vs control") +
  ggplot2::ylab("-log10(padj)") +
  #Set title
  ggplot2::ggtitle("disease vs control")+
  #Add a vertical line for fold change cut-offs
  ggplot2::geom_vline(xintercept=c(-1, 1), linetype="longdash", colour="black", size=0.4) +
  #Add a horizontal line for P-value cut-off
  ggplot2::geom_hline(yintercept=1.3, linetype="longdash", colour="black", size=0.4)

volcano_plot

####2####

# load the list of TSGs and oncogenes
cosmic_genes <- read.csv("data/Census.csv")
head(cosmic_genes)

cosmic_genes$Role.in.Cancer <- as.character(cosmic_genes$Role.in.Cancer)
cosmic_genes$Role.in.Cancer <- strsplit(cosmic_genes$Role.in.Cancer, ', ')


# merge the list of cosmic genes with the list of differential expressed genes
merged_df <- merge(res_filtered_df, cosmic_genes, by.x = "row.names", by.y = "Gene.Symbol")

merge_upregulated_genes <- merge(as.data.frame(upregulated_genes), cosmic_genes, by.x = "row.names", by.y = "Gene.Symbol")

merge_downregulated_genes <- merge(as.data.frame(downregulated_genes), cosmic_genes, by.x = "row.names", by.y = "Gene.Symbol")



# # determine entries for contingency table to assess statistical significance of TSG in upregulated genes 
TSG_up <- sum(grepl("TSG", merge_upregulated_genes$Role.in.Cancer))
TSG_notup <- sum(grepl("TSG", merged_df$Role.in.Cancer)) - sum(grepl("TSG", merge_upregulated_genes$Role.in.Cancer))
notTSG_up <- nrow(upregulated_genes) - sum(grepl("TSG", merge_upregulated_genes$Role.in.Cancer))
notTSG_notup <- nrow(res_filtered_df) - TSG_up - TSG_notup - notTSG_up



# Create the contingency table
contingency_table_TSG_up <- matrix(c(TSG_up, TSG_notup, notTSG_up, notTSG_notup), nrow = 2, byrow = TRUE)
contingency_table_TSG_up

# Perform the Fisher's exact test on the contingency table
fisher_TSG_up <- fisher.test(contingency_table_TSG_up)
fisher.test(contingency_table_TSG_up)

# # determine entries for contingency table to assess statistical significance of TSG in downregulated genes 
TSG_down <- sum(grepl("TSG", merge_downregulated_genes$Role.in.Cancer))
TSG_notdown <- sum(grepl("TSG", merged_df$Role.in.Cancer)) - sum(grepl("TSG", merge_downregulated_genes$Role.in.Cancer))
notTSG_down <- nrow(downregulated_genes) - sum(grepl("TSG", merge_downregulated_genes$Role.in.Cancer))
notTSG_notdown <- nrow(res_filtered_df) - TSG_down - TSG_notdown - notTSG_down



# Create the contingency table
contingency_table_TSG_down <- matrix(c(TSG_down, TSG_notdown, notTSG_down, notTSG_notdown), nrow = 2, byrow = TRUE)
contingency_table_TSG_down

# Perform the Fisher's exact test on the contingency table
fisher_TSG_down <- fisher.test(contingency_table_TSG_down)
fisher.test(contingency_table_TSG_down)

# # determine entries for contingency table to assess statistical significance of oncogenes in upregulated genes 
Onc_up <- sum(grepl("oncogene", merge_upregulated_genes$Role.in.Cancer))
Onc_notup <- sum(grepl("oncogene", merged_df$Role.in.Cancer)) - sum(grepl("oncogene", merge_upregulated_genes$Role.in.Cancer))
notOnc_up <- nrow(upregulated_genes) - sum(grepl("oncogene", merge_upregulated_genes$Role.in.Cancer))
notOnc_notup <- nrow(res_filtered_df) - Onc_up - Onc_notup - notOnc_up



# Create the contingency table
contingency_table_Onc_up <- matrix(c(Onc_up, Onc_notup, notOnc_up, notOnc_notup), nrow = 2, byrow = TRUE)
contingency_table_Onc_up

# Perform the Fisher's exact test on the contingency table
fisher_Onc_up <- fisher.test(contingency_table_Onc_up)
fisher.test(contingency_table_Onc_up)


# # determine entries for contingency table to assess statistical significance of oncogenes in downregulated genes 
Onc_down <- sum(grepl("oncogene", merge_downregulated_genes$Role.in.Cancer))
Onc_notdown <- sum(grepl("oncogene", merged_df$Role.in.Cancer)) - sum(grepl("oncogene", merge_downregulated_genes$Role.in.Cancer))
notOnc_down <- nrow(downregulated_genes) - sum(grepl("oncogene", merge_downregulated_genes$Role.in.Cancer))
notOnc_notdown <- nrow(res_filtered_df) - Onc_down - Onc_notdown - notOnc_down



# Create the contingency table
contingency_table_Onc_down <- matrix(c(Onc_down, Onc_notdown, notOnc_down, notOnc_notdown), nrow = 2, byrow = TRUE)
contingency_table_Onc_down

# Perform the Fisher's exact test on the contingency table
fisher_Onc_down <- fisher.test(contingency_table_Onc_down)
fisher.test(contingency_table_Onc_down)


# # determine entries for contingency table to assess statistical significance of TSG&oncogenes in upregulated genes 
TSG_Onc_up <- sum(grepl("oncogene", merge_upregulated_genes$Role.in.Cancer)) & grepl("TSG", merge_upregulated_genes$Role.in.Cancer)
TSG_Onc_notup <- (sum(grepl("oncogene", merged_df$Role.in.Cancer)) & sum(grepl("TSG", merged_df$Role.in.Cancer))) - TSG_Onc_up
notTSG_Onc_up <- nrow(upregulated_genes) - TSG_Onc_up
notTSG_Onc_notup <- nrow(res_filtered_df) - TSG_Onc_up - TSG_Onc_notup - notTSG_Onc_up



# Create the contingency table
contingency_table_TSG_Onc_up <- matrix(c(TSG_Onc_up, TSG_Onc_notup, notTSG_Onc_up, notTSG_Onc_notup), nrow = 2, byrow = TRUE)
contingency_table_TSG_Onc_up

# Perform the Fisher's exact test on the contingency table
fisher_TSG_Onc_up <- fisher.test(contingency_table_TSG_Onc_up)
fisher.test(contingency_table_TSG_Onc_up)


# # determine entries for contingency table to assess statistical significance of TSG&oncogenes in upregulated genes 
TSG_Onc_down <- sum(grepl("oncogene", merge_downregulated_genes$Role.in.Cancer)) & grepl("TSG", merge_downregulated_genes$Role.in.Cancer)
TSG_Onc_notdown <- (sum(grepl("oncogene", merged_df$Role.in.Cancer)) & sum(grepl("TSG", merged_df$Role.in.Cancer))) - TSG_Onc_down
notTSG_Onc_down <- nrow(downregulated_genes) - TSG_Onc_down
notTSG_Onc_notdown <- nrow(res_filtered_df) - TSG_Onc_down - TSG_Onc_notdown - notTSG_Onc_down



# Create the contingency table
contingency_table_TSG_Onc_down <- matrix(c(TSG_Onc_down, TSG_Onc_notdown, notTSG_Onc_down, notTSG_Onc_notdown), nrow = 2, byrow = TRUE)
contingency_table_TSG_Onc_down

# Perform the Fisher's exact test on the contingency table
fisher_TSG_Onc_down <- fisher.test(contingency_table_TSG_Onc_down)
fisher.test(contingency_table_TSG_Onc_down)

#######################################################################################################################################
## Visualization of results

# Create data frames for upregulated and downregulated genes
upregulated_data <- data.frame(
  Category = c("TSG", "Oncogene", "TSG&Oncogene"),
  Direction = rep("Upregulated", 3),
  P.value = c(fisher_TSG_up$p.value, fisher_Onc_up$p.value, fisher_TSG_Onc_up$p.value),
  Significant = c(fisher_TSG_up$p.value <= 0.05, fisher_Onc_up$p.value <= 0.05, fisher_TSG_Onc_up$p.value <= 0.05)
)

downregulated_data <- data.frame(
  Category = c("TSG", "Oncogene", "TSG&Oncogene"),
  Direction = rep("Downregulated", 3),
  P.value = c(fisher_TSG_down$p.value, fisher_Onc_down$p.value, fisher_TSG_Onc_down$p.value),
  Significant = c(fisher_TSG_down$p.value <= 0.05, fisher_Onc_down$p.value <= 0.05, fisher_TSG_Onc_down<= 0.05)
)

# Combine upregulated and downregulated data into a single data frame
combined_data <- rbind(upregulated_data, downregulated_data)

# Define bar plot to access significant of TSG and Oncogenes in Upregulated and Downregulated Genes
ggplot(na.omit(combined_data), aes(x = Category, y = P.value, fill = Significant)) +
  geom_point(position = position_dodge(width = 0.75), aes(color = Significant)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ Direction, scales = "free") +
  ylab("P-value") +
  ggtitle("Significance of TSG and Oncogenes in Upregulated and Downregulated Genes") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(aes(label = ifelse(Significant, sprintf("%.2e", P.value), ""), color = Significant), 
            position = position_dodge(width = 0.75), 
            vjust = -0.5, size = 3.5)