
# =========================
# Venn Diagram and enrichment 
# =========================
# Load libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(readxl)
library(writexl)

# Input files
file_DT3 <- "data/HCC38_DT3_DEG_filtered_results.xlsx"
file_GT3 <- "data/HCC38_GT3_DEG_filtered_results.xlsx"

# Output folder
output_dir <- "DEG_pipeline_output"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Column names in your DEG table
# Modify these according to your file
gene_col   <- "Gene"        # gene symbol column
fc_col     <- "FC"      
pval_col   <- "PValue"     # raw p-value column

# Thresholds
fc_cutoff <- 2              
pval_cutoff <- 0.05

# Organism for enrichment
organism_db <- org.Hs.eg.db
kegg_org <- "hsa"           

# Read data
deg_DT3 <- read_excel(file_DT3) %>% as.data.frame()
deg_GT3 <- read_excel(file_GT3) %>% as.data.frame()

# Standardize columns (same format as your pipeline)
process_deg <- function(df) {
  df %>%
    dplyr::select(all_of(c(gene_col, fc_col, pval_col))) %>%
    dplyr::rename(
      Gene = all_of(gene_col),
      FC = all_of(fc_col),
      PValue = all_of(pval_col)
    ) %>%
    filter(!is.na(Gene), !is.na(FC), !is.na(PValue)) %>%
    mutate(
      Regulation = case_when(
        FC >= fc_cutoff & PValue < pval_cutoff  ~ "Up",
        FC <= -fc_cutoff & PValue < pval_cutoff ~ "Down",
        TRUE ~ "NotSig"
      )
    )
}

deg_DT3 <- process_deg(deg_DT3)
deg_GT3 <- process_deg(deg_GT3)

# Upregulated
up_DT3 <- deg_DT3 %>% filter(Regulation == "Up") %>% pull(Gene) %>% unique()
up_GT3 <- deg_GT3 %>% filter(Regulation == "Up") %>% pull(Gene) %>% unique()

# Downregulated
down_DT3 <- deg_DT3 %>% filter(Regulation == "Down") %>% pull(Gene) %>% unique()
down_GT3 <- deg_GT3 %>% filter(Regulation == "Down") %>% pull(Gene) %>% unique()

# Overlaps
overlap_up <- intersect(up_DT3, up_GT3)
overlap_down <- intersect(down_DT3, down_GT3)

cat("Overlap Up:", length(overlap_up), "\n")
cat("Overlap Down:", length(overlap_down), "\n")

# =========================
# Venn Diagram
# =========================

# Install if needed
if (!requireNamespace("eulerr", quietly = TRUE)) install.packages("eulerr")
library(eulerr)

# Downregulated Venn
# Counts for Euler/Venn
fit_down <- euler(c(
  "δT3&γT3" = length(overlap_down),
  δT3 = length(down_DT3) -length(overlap_down),
  γT3 = length(down_GT3)-length(overlap_down)
))

# Plot
png(file.path(output_dir, "HCC_Venn_Downregulated_aesthetic_1.png"),
    width = 2200, height = 2200, res = 300)

plot(
  fit_down,
  fills = list(fill = c("pink", "darkgreen"), alpha = 0.2),
  edges = list(col = c("pink", "darkgreen"), lwd = 8),
  labels = list(font = 4, cex = 2.8),
  quantities = list(cex = 2.5),
  legend = FALSE
)

dev.off()


# Upregulated Venn
# Counts for Euler/Venn
fit_up <- euler(c(
  "δT3&γT3" = length(overlap_up),
  δT3 = length(up_DT3) -length(overlap_up),
  γT3 = length(up_GT3)-length(overlap_up)
))

# Plot
png(file.path(output_dir, "HCC_Venn_Upregulated_aesthetic_1.png"),
    width = 2200, height = 2200, res = 300)

plot(
  fit_up,
  fills = list(fill = c("pink", "darkgreen"), alpha = 0.2),
  edges = list(col = c("pink", "darkgreen"), lwd = 8),
  labels = list(font = 4, cex = 2.5),
  quantities = list(cex = 2.5),
  legend = FALSE
)

dev.off()

# =========================
# Enrichment analysis of common genes
# =========================


convert_entrez <- function(genes) {
  bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = organism_db
  )$ENTREZID %>% unique()
}

entrez_up <- convert_entrez(overlap_up)
entrez_down <- convert_entrez(overlap_down)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
library(stringr)

run_go_kegg_singleplot <- function(gene_list, label) {
  
  # -------------------------
  # GO enrichment
  # -------------------------
  go_res <- enrichGO(
    gene = gene_list,
    OrgDb = organism_db,
    ont = "ALL",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  go_df <- as.data.frame(go_res)
  if (nrow(go_df) > 0) {
    go_df <- go_df %>%
      arrange(pvalue) %>%
      slice_head(n = 10) %>%
      mutate(
        Database = "GO",
        negLogP = -log10(pvalue)
      )
  } else {
    go_df <- data.frame()
  }
  
  # -------------------------
  # KEGG enrichment
  # -------------------------
  kegg_res <- enrichKEGG(
    gene = gene_list,
    organism = kegg_org,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  )
  
  kegg_df <- as.data.frame(kegg_res)
  if (nrow(kegg_df) > 0) {
    kegg_df <- kegg_df %>%
      arrange(pvalue) %>%
      slice_head(n = 10) %>%
      mutate(
        Database = "KEGG",
        negLogP = -log10(pvalue)
      )
  } else {
    kegg_df <- data.frame()
  }
  
  # -------------------------
  # Combine GO + KEGG
  # -------------------------
  combined_df <- bind_rows(go_df, kegg_df)
  
  if (nrow(combined_df) == 0) {
    message("No significant GO or KEGG terms found for ", label)
    return(NULL)
  }
  
  combined_df <- combined_df %>%
    mutate(
      Description_wrap = str_wrap(Description, width = 50)
    ) %>%
    arrange(negLogP) %>%
    mutate(
      Description_wrap = factor(Description_wrap, levels = Description_wrap)
    )
  
  # Save combined table
  write_xlsx(
    list(
      GO_top10 = go_df,
      KEGG_top10 = kegg_df,
      Combined_plot_data = combined_df
    ),
    file.path(output_dir, paste0("GO_KEGG_", label, "_singleplot_data_1.xlsx"))
  )
  
  # -------------------------
  # Single combined plot
  # -------------------------
  combined_plot <- ggplot(
    combined_df,
    aes(
      x = negLogP,
      y = Description_wrap,
      fill = Database
    )
  ) +
    geom_col(width = 0.75) +
    
    # Add labels directly on bars
    geom_text(
      aes(
        x = 0,                     # anchor at start of bar
        label = Description_wrap
      ),
      hjust = 0,                  # left aligned
      nudge_x = 0.05,             # small offset into bar
      color = "black",
      size = 3.5
    ) +
    
    scale_fill_manual(
      values = c("GO" = "#FFF090", "KEGG" = "#7FC7FF")
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.15))  # extra space for labels
    ) +
    
    labs(
      title = paste("GO and KEGG Enrichment -", label),
      x = expression(-log[10](Pvalue)),
      y = NULL,
      fill = "Database"
    ) +
    theme_bw(base_size = 15) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 14),
      
      # Remove y-axis text completely
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      axis.text.x = element_text(size = 12, color = "black"),
      
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) 
  
  ggsave(
    file.path(output_dir, paste0("GO_KEGG_", label, "_singleplot_1.png")),
    combined_plot,
    width = 5,
    height = 4,
    dpi = 600
  )
  
  return(combined_plot)
}

plot_up <- run_go_kegg_singleplot(entrez_up, "Overlap_Up_1")
plot_down <- run_go_kegg_singleplot(entrez_down, "Overlap_Down_1")

plot_up
plot_down

