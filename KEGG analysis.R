
# =========================
# KEGG pathway clusterprofiler - Run the DEG_pipeline R script beforehand till step 4
# =========================

#Load library
library(clusterProfiler)
library(org.Hs.eg.db)

# Enrichment analysis
##Convert gene symbol to entrez id for enrichment.
##simply replace the value of gene sets for separate enrichment analysis
up_deg_entrez <- bitr( 
  up_deg$Gene, 
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db) #replace value name
up_deg_entrez <- up_deg_entrez$ENTREZID

down_deg_entrez <- bitr( 
  down_deg$Gene, 
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db) #replace value name
down_deg_entrez <- down_deg_entrez$ENTREZID

# =========================
# 6. Run KEGG enrichment
# =========================
kegg_up <- NULL
kegg_down <- NULL

if (length(up_deg_entrez) > 0) {
  kegg_up <- enrichKEGG(
    gene = up_deg_entrez,
    organism = kegg_org,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}

if (length(down_deg_entrez) > 0) {
  kegg_down <- enrichKEGG(
    gene = down_deg_entrez,
    organism = kegg_org,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}

# =========================
# 7. Convert KEGG results to data frames
# =========================
kegg_up_df <- if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
  as.data.frame(kegg_up) %>% mutate(Regulation = "Up")
} else {
  data.frame()
}

kegg_down_df <- if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
  as.data.frame(kegg_down) %>% mutate(Regulation = "Down")
} else {
  data.frame()
}

# Save full KEGG results
write_xlsx(
  list(
    KEGG_Up = kegg_up_df,
    KEGG_Down = kegg_down_df
  ),
  file.path(output_dir, "HCC38_GT3_KEGG_results.xlsx")
)

# =========================
# 8. Select top 10 significant pathways for each group
# by raw p-value; if you prefer adjusted p-value, change arrange(pvalue)
# =========================
top_kegg_up <- kegg_up_df %>%
  arrange(pvalue) %>%
  slice_head(n = 10)

top_kegg_down <- kegg_down_df %>%
  arrange(pvalue) %>%
  slice_head(n = 10)

kegg_plot_df <- bind_rows(top_kegg_up, top_kegg_down)

# Stop if nothing to plot
if (nrow(kegg_plot_df) == 0) {
  stop("No KEGG pathways found for either upregulated or downregulated genes.")
}

# =========================
# 9. Prepare plotting columns
# dot size   = Count
# dot shape  = FDR category
# color      = -log10(pvalue)
# =========================

kegg_plot_df <- kegg_plot_df %>%
  mutate(
    negLogP = -log10(pvalue),
    FDR_shape = ifelse(p.adjust < 0.05, "FDR < 0.05", "FDR >= 0.05"),
    Description_wrap = stringr::str_wrap(Description, width = 40),
    
    # Convert GeneRatio (e.g. "5/120") to numeric
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    })
  )

# Order terms within each facet
top_kegg_up <- kegg_plot_df %>%
  filter(Regulation == "Up") %>%
  arrange(GeneRatio_num) %>%
  mutate(Description_wrap = factor(Description_wrap, levels = Description_wrap))

top_kegg_down <- kegg_plot_df %>%
  filter(Regulation == "Down") %>%
  arrange(GeneRatio_num) %>%
  mutate(Description_wrap = factor(Description_wrap, levels = Description_wrap))

kegg_plot_df <- bind_rows(top_kegg_up, top_kegg_down)

# =========================
# 10. Combined KEGG dotplot
# =========================
kegg_combined_dotplot <- ggplot(
  kegg_plot_df,
  aes(
    x = GeneRatio_num,
    y = Description_wrap,
    size = Count,
    shape = FDR_shape,
    color = negLogP
  )
) +
  geom_point(alpha = 0.9) +
  facet_grid(
    Regulation ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_size_continuous(
    name = "Gene count",
    range = c(4, 10)
  ) +
  scale_shape_manual(
    name = "FDR",
    values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 15)
  ) +
  scale_color_gradient(
    name = expression(-log[10](Pvalue)),
    low = "skyblue3",
    high = "red"
  ) +
  labs(
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", lineheight = 0.95),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 13),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  guides(
    size = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    color = guide_colorbar(order = 3)
  )

kegg_combined_dotplot

ggsave(
  filename = file.path(output_dir, "HCC38_GT3_KEGG_top10_dotplot.png"),
  plot = kegg_combined_dotplot,
  width = 9,
  height = 7,
  dpi = 600
)