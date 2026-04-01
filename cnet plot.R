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

# Column names
gene_col <- "User ID"
fc_col   <- "FC"
pval_col <- "PValue"

# Organism
organism_db <- org.Hs.eg.db
kegg_org <- "hsa"

# -----------------------------------
# 1. Get all overlap genes (both Up and Down together)
# -----------------------------------
deg_all <- HCC38_GT3_TCGA %>%
  filter(Regulation %in% c("Up", "Down")) %>%
  pull(.data[[gene_col]]) %>%
  unique()

# -----------------------------------
# 2. Prepare fold-change vector for all genes
# names(fc_vector) must be ENTREZID
# values remain signed FC, so Up/Down are colored properly
# -----------------------------------
prepare_fc_vector <- function(df, genes) {
  
  mapping <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = organism_db
  )
  
  df_fc <- df %>%
    filter(.data[[gene_col]] %in% genes) %>%
    inner_join(mapping, by = setNames("SYMBOL", gene_col)) %>%
    distinct(.data[[gene_col]], .keep_all = TRUE)
  
  fc_vector <- df_fc[[fc_col]]
  names(fc_vector) <- df_fc[[gene_col]]   # use gene symbols, not ENTREZID
  
  return(fc_vector)
}

fc_all <- prepare_fc_vector(HCC38_GT3_TCGA, deg_all)

# -----------------------------------
# 3. Convert all genes to ENTREZ IDs
# -----------------------------------
entrez_all <- bitr(
  deg_all,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = organism_db
)$ENTREZID %>%
  unique()

# -----------------------------------
# 4. Sanity checks
# -----------------------------------
cat("Total genes used:", length(deg_all), "\n")
cat("Mapped ENTREZ IDs:", length(entrez_all), "\n")
cat("Fold-change vector length:", length(fc_all), "\n")


# -----------------------------------
# 5. Run GO and KEGG enrichment
# -----------------------------------
go_res <- enrichGO(
  gene = entrez_all,
  OrgDb = organism_db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

kegg_res <- enrichKEGG(
  gene = entrez_all,
  organism = kegg_org,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_res <- setReadable(kegg_res, OrgDb = organism_db, keyType = "ENTREZID")


# -----------------------------------
# 6. Plot combined cnet plot
# -----------------------------------
p_go <- enrichplot::cnetplot(
  go_res,
  showCategory = 4,
  foldChange = fc_all,
  node_label = "all",
  color_category = "#F0E442"
)

p_go

p_kegg <- enrichplot::cnetplot(
  kegg_res,
  showCategory = 4,
  foldChange = fc_all,
  node_label = "all",
  color_category = "#F0E442"
)

p_kegg

ggsave(
  filename = file.path(output_dir, "HCC38_GT3_DAVID_GO_cnetplot.png"),
  plot = p_go,
  width = 4,
  height = 4,
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "HCC38_GT3_DAVID_KEGG_cnetplot.png"),
  plot = p_kegg,
  width = 4,
  height = 4,
  dpi = 600
)
