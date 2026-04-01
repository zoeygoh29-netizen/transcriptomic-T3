# ============================================
# RNA-Seq Data Analysis Pipeline
# Volcano Plot
# ============================================

# Install required packages (run once)
install_packages <- function() {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", 
                         "enrichplot", "pathview", "DOSE"))
  install.packages(c("ggplot2", "dplyr", "tidyr", "writexl", "readr", "ggrepel"))
}

# Load libraries
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

input_file <- "data/HCC_GT3_ALL.xlsx" #replace file name here

# Sheet name or index
sheet_to_read <- 1

# Output folder
output_dir <- "DEG_pipeline_output"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Column names in your DEG table
# Modify these according to your file
gene_col   <- "external_gene_name"        # gene symbol column
fc_col     <- "Fold.Change"      
pval_col   <- "P.val"     # raw p-value column

# Thresholds
fc_cutoff <- 2              
pval_cutoff <- 0.05

# Organism for enrichment
organism_db <- org.Hs.eg.db
kegg_org <- "hsa"           

deg_data <- read_excel(input_file, sheet = sheet_to_read) %>%
  as.data.frame()

# Basic check
cat("Imported rows:", nrow(deg_data), "\n")
cat("Imported columns:\n")
print(colnames(deg_data))

# Keep only required columns and remove NA
deg_data <- deg_data %>%
  dplyr::select(all_of(c(gene_col, fc_col, pval_col))) %>%
  dplyr::rename(
    Gene = all_of(gene_col),
    FC = all_of(fc_col),
    PValue = all_of(pval_col)
  ) %>%
  filter(!is.na(Gene), !is.na(FC), !is.na(PValue))

# =========================
# 4. Filter significant DEGs
# =========================
deg_data <- deg_data %>%
  mutate(
    Regulation = case_when(
      FC >= fc_cutoff & PValue < pval_cutoff  ~ "Up",
      FC <= -fc_cutoff & PValue < pval_cutoff ~ "Down",
      TRUE ~ "NotSig"
    )
  )

sig_deg <- deg_data %>%
  filter(Regulation != "NotSig")

up_deg <- sig_deg %>% filter(Regulation == "Up")
down_deg <- sig_deg %>% filter(Regulation == "Down")

cat("Significant DEGs:", nrow(sig_deg), "\n")
cat("Upregulated:", nrow(up_deg), "\n")
cat("Downregulated:", nrow(down_deg), "\n")

# Save DEG tables
write_xlsx(
  list(
    All_DEGs = deg_data,
    Significant_DEGs = sig_deg,
    Upregulated = up_deg,
    Downregulated = down_deg
  ),
  file.path(output_dir, "HCC38_GT3_DEG_filtered_results.xlsx")
)

# =========================
# 5. Add rank score
# rank score = abs(log2FC) * -log10(PValue)
# =========================
deg_data <- deg_data %>%
  mutate(
    RankScore = abs(FC) * (-log10(PValue))
  )

# =========================
# 6. Select top 10 up/down genes for labeling
# =========================
top_up <- deg_data %>%
  filter(Regulation == "Up") %>%
  arrange(desc(RankScore)) %>%
  slice_head(n = 10)

top_down <- deg_data %>%
  filter(Regulation == "Down") %>%
  arrange(desc(RankScore)) %>%
  slice_head(n = 10)

top_labeled_genes <- bind_rows(top_up, top_down)

# =========================
# 7. Volcano plot
# =========================
volcano_plot <- ggplot(deg_data, aes(x = FC, y = -log10(PValue), color = Regulation)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey") +
  geom_point(
    data = top_labeled_genes,
    aes(x = FC, y = -log10(PValue)),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    data = top_labeled_genes,
    aes(label = Gene),
    size = 4,
    color = "black",              # label text color
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  coord_cartesian(xlim = c(-30, 12)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",      # move legend to top
    text = element_text(color = "black"),   # all text black
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black")
  ) +
  labs(
    x = "Fold Change",
    y = "-log10(p-value)",
    color = "Regulation"
  )

volcano_plot

ggsave(
  filename = file.path(output_dir, "HCC38_GT3_Volcano_plot.png"),
  plot = volcano_plot,
  width = 9,
  height = 7,
  dpi = 600
)
