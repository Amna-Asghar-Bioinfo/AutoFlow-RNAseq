#' @title RNA-seq Biological Significance Filtering
#' @description 
#' This script performs post-statistical filtering on DESeq2 results. It extracts 
#' biologically relevant gene sets by applying strict thresholds on adjusted 
#' p-values and fold-change magnitudes. The script bifurcates the results into 
#' upregulated and downregulated gene lists for downstream pathway analysis.
#' 
#' @author Amna Asghar
#' @usage Rscript Filtration.R <input_file> <up_output> <down_output>
#' @param input_file Path to the annotated DESeq2 results (CSV).
#' @param up_output Path to save genes with Log2FC >= 2.0 and padj < 0.05.
#' @param down_output Path to save genes with Log2FC <= -2.0 and padj < 0.05.
# ------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------
library(dplyr)

# ------------------------------------------------------------
# 2. READ INPUT FROM COMMAND LINE
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
up_output <- args[2]
down_output <- args[3]

# ------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

data <- read.csv(input_file)

# ------------------------------------------------------------
# 4. FILTER DATA (BIOLOGICAL SIGNIFICANCE FILTERS)
# ------------------------------------------------------------

filtered_data <- data %>% 
  filter(
    !is.na(padj) &
      !is.na(log2FoldChange) &
      !is.na(Gene.Name) &
      padj < 0.05 &
      abs(log2FoldChange) >= 2.0
  )

# ------------------------------------------------------------
# 5. SPLIT UP AND DOWN REGULATED GENES
# ------------------------------------------------------------

up_regulated <- filtered_data %>% 
  filter(log2FoldChange >= 2.0) %>%
  arrange(desc(log2FoldChange))

down_regulated <- filtered_data %>% 
  filter(log2FoldChange <= -2.0) %>%
  arrange(log2FoldChange)

# ------------------------------------------------------------
# 6. SAVE OUTPUT FILES
# ------------------------------------------------------------

write.csv(up_regulated, up_output, row.names = FALSE)
write.csv(down_regulated, down_output, row.names = FALSE)

# ------------------------------------------------------------
# 7. SUMMARY OUTPUT (FOR TERMINAL LOGGING)
# ------------------------------------------------------------

cat("\n--- Filtering Summary ---\n")
cat("Total significant genes: ", nrow(filtered_data), "\n")
cat("Upregulated genes: ", nrow(up_regulated), "\n")
cat("Downregulated genes: ", nrow(down_regulated), "\n")
cat("Outputs saved successfully.\n")