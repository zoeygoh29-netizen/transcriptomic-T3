
#Load libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # Change based on organism: org.Mm.eg.db for mouse
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(readxl)
library(writexl)
library(stringr)


# Output folder
output_dir <- "DEG_pipeline_output"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Import TCGA BRCA data
tcga_basal <- read.table("data/tcga-basal.txt", header = FALSE, sep = "\t")


# Column names in your DEG table
# Modify these according to your file
gene_col   <- "Gene"        # gene symbol column
fc_col     <- "FC"      
pval_col   <- "PValue"     # raw p-value column


# Organism for enrichment
organism_db <- org.Hs.eg.db
kegg_org <- "hsa"  

deg_DT3 <- MDAMB231_DT3_DEG_filtered_results %>% pull(gene_col) %>% unique()
deg_GT3 <- MDAMB231_GT3_DEG_filtered_results %>% pull(gene_col) %>% unique()
deg_TCGA <- tcga_basal %>% pull(V1) %>% unique()

# =========================
# Venn Diagram
# =========================

# Install if needed
# Install if needed
if (!requireNamespace("eulerr", quietly = TRUE)) install.packages("eulerr")
library(eulerr)

# Build Euler object (FIXED TCGA line)
fit_ <- euler(list(
  "δT3" = deg_DT3,
  "γT3" = deg_GT3,
  "TCGA" = deg_TCGA
))

# Output file
outfile <- file.path(output_dir, "MDAMB231_Venn_TCGA_aesthetic.png")

png(outfile, width = 2200, height = 2200, res = 300)

plot(
  fit_,
  fills = list(fill = c("pink", "darkgreen", "lightblue"), alpha = 0.3),
  edges = list(col = c("pink", "darkgreen", "lightblue"), lwd = 4),
  labels = list(font = 3, cex = 1.5),
  quantities = list(cex = 1.5),
  legend = FALSE
)

dev.off()

#Overlap
overlap_up <- intersect((MDAMB231_DT3_DEG_filtered_results %>% filter(Regulation == "Up") %>% pull(gene_col) %>% unique()), deg_TCGA)
overlap_down <- intersect((MDAMB231_DT3_DEG_filtered_results %>% filter(Regulation == "Down") %>% pull(gene_col) %>% unique()), deg_TCGA)

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

run_enrichment_combined <- function(gene_list, label) {
  
  # ---------- GO ----------
  go_res <- enrichGO(
    gene = gene_list,
    OrgDb = organism_db,
    ont = "ALL",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  go_df <- as.data.frame(go_res)
  
  if (nrow(go_df) > 0) {
    go_df <- go_df %>%
      arrange(p.adjust, pvalue) %>%
      slice_head(n = 10) %>%
      mutate(
        Database = "GO",
        Regulation = label,
        negLogP = -log10(pvalue)
      )
  } else {
    go_df <- data.frame()
  }
  
  # ---------- KEGG ----------
  kegg_res <- enrichKEGG(
    gene = gene_list,
    organism = kegg_org,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  kegg_df <- as.data.frame(kegg_res)
  
  if (nrow(kegg_df) > 0) {
    kegg_df <- kegg_df %>%
      arrange(p.adjust, pvalue) %>%
      slice_head(n = 10) %>%
      mutate(
        Database = "KEGG",
        Regulation = label,
        negLogP = -log10(pvalue)
      )
  } else {
    kegg_df <- data.frame()
  }
  
  bind_rows(go_df, kegg_df)
}

df_up <- run_enrichment_combined(entrez_up, "Up")
df_down <- run_enrichment_combined(entrez_down, "Down")

combined_df <- bind_rows(df_up, df_down)

combined_df <- combined_df %>%
  mutate(
    Description_wrap = str_wrap(Description, width = 100)
  )

# Order within each facet
combined_df <- combined_df %>%
  group_by(Regulation) %>%
  arrange(negLogP, .by_group = TRUE) %>%
  mutate(
    Description_wrap = factor(Description_wrap, levels = unique(Description_wrap))
  ) %>%
  ungroup()

combined_plot <- ggplot(
  combined_df,
  aes(
    x = negLogP,
    y = Description_wrap,
    fill = Database
  )
) +
  geom_col(width = 0.75) +
  
  # Labels INSIDE bars (left aligned)
  geom_text(
    aes(x = 0, label = Description_wrap),
    hjust = 0,
    nudge_x = 0.05,
    color = "black",
    size = 3.5
  ) +
  
  facet_grid(
    Regulation ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_manual(
    values = c("GO" = "#FFF090", "KEGG" = "#7FC7FF")
  ) +
  
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    x = expression(-log[10](Pvalue)),
    y = NULL,
    fill = "Database"
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

combined_plot

ggsave(
  file.path(output_dir, "MDAMB231_DT3_GO_KEGG_Overlap_TCGA_combined.png"),
  combined_plot,
  width = 8,
  height = 5,
  dpi = 600
)
