# -------------------------------
# Snakemake-compatible DESeq2 script
# -------------------------------
#' @title RNA-seq Differential Expression Pipeline (DESeq2)
#' @description 
#' This script automates the differential expression analysis for RNA-seq data.
#' It processes raw counts, standardizes metadata genotype labels (SNAI1 knockout 
#' vs Wildtype), performs statistical testing, and maps results to genomic 
#' coordinates using biomaRt. It also generates MA and Volcano plots for QC.
#' 
#' @author Amna Asghar
#' @usage Rscript Analyzing_with_DESeq.R <counts_file> <metadata_file> <output_file>
#' @param counts_file Input path to the gene count matrix (tab-delimited).
#' @param metadata_file Input path to the sample metadata (tab-delimited).
#' @param output_file Output path for the final annotated CSV result.

# arguments from Snakemake
args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
metadata_file <- args[2]
output_file <- args[3]

# -------------------------------
# Load libraries
# -------------------------------
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(biomaRt)

# -------------------------------
# Read input files
# -------------------------------
counts <- read.delim(counts_file)
metadata <- read.delim(metadata_file)

# -------------------------------
# Data wrangling
# -------------------------------

# Set gene IDs as rownames
rownames(counts) <- counts$Gene.ID

# Save gene info separately
genes <- counts[, c("Gene.ID", "Gene.Name")]

# Remove unnecessary columns
counts <- counts[, -c(1,2)]

# Set metadata rownames
rownames(metadata) <- metadata$Run

# Keep only genotype column
metadata <- metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]

# Rename column
colnames(metadata) <- c("genotype")

# Clean genotype labels
metadata$genotype[metadata$genotype == "wild type genotype"] <- "wildtype"
metadata$genotype[metadata$genotype == "Snai1 knockout"] <- "knockout"

# Convert to factor
metadata$genotype <- factor(metadata$genotype, levels = c("wildtype", "knockout"))

# -------------------------------
# DESeq2 analysis
# -------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ genotype
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10,]

# Run DESeq
dds <- DESeq(dds)

# Get results
res <- results(dds,
               contrast = c("genotype", "knockout", "wildtype"),
               alpha = 1e-5)

# Convert to dataframe
res_df <- as.data.frame(res)

# Merge with gene names
res_df$Gene.ID <- rownames(res_df)
res_df <- merge(res_df, genes, by = "Gene.ID")

# -------------------------------
# Visualization (optional but kept)
# -------------------------------

# MA plot
png("results/MA_plot.png")
plotMA(res)
dev.off()

# Volcano plot
png("results/volcano.png")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'SNAI1 Knockout vs Wildtype',
                pCutoff = 10e-6,
                FCcutoff = 1.5)
dev.off()

# -------------------------------
# biomaRt annotation
# -------------------------------

options(timeout = 300)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

attributes <- c('ensembl_gene_id', 'chromosome_name',
                'start_position', 'end_position')

all.genes <- getBM(attributes = attributes, mart = ensembl)

# Clean IDs
all.genes$ensembl_gene_id <- trimws(all.genes$ensembl_gene_id)
res_df$Gene.ID_clean <- gsub("\\..*", "", res_df$Gene.ID)

# Merge annotation
merged_data <- merge(all.genes,
                     res_df,
                     by.x = "ensembl_gene_id",
                     by.y = "Gene.ID_clean")

# -------------------------------
# Save outputs
# -------------------------------

if (nrow(merged_data) > 0) {
  
  merged_data$chromosome_name <- paste("chr", merged_data$chromosome_name, sep="")
  
  write.csv(merged_data, output_file, row.names = FALSE)
  
  # Save subset
  genes_to_check <- c('THY1','SFMBT2','PASD1')
  subset_data <- merged_data[merged_data$Gene.Name %in% genes_to_check, ]
  
  write.csv(subset_data, "results/deseq_subset.csv", row.names = FALSE)
  
  message("Analysis completed successfully.")
  
} else {
  stop("Merge failed! No matching IDs found.")
}