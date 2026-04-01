
# =========================
# 8. GO Enrichment Analysis - Run the DEG_pipeline R script beforehand till step 4
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

ego_up <- enrichGO(
  gene          = up_deg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",          # BP + MF + CC
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_up <- as.data.frame(ego_up)

down_deg_entrez <- bitr( 
  down_deg$Gene, 
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db) #replace value name
down_deg_entrez <- down_deg_entrez$ENTREZID

ego_down <- enrichGO(
  gene          = down_deg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",          # BP + MF + CC
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_down <- as.data.frame(ego_down)

write_xlsx(
  list(
    GO_Up = ego_up,
    GO_down = ego_down
  ),
  file.path(output_dir, "MDA_GT3_GO.xlsx")
)

#Prepare GO data for plotting
ego_up_top <- ego_up %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

ego_up_top <- ego_up_top %>%
  group_by(ONTOLOGY) %>%
  mutate(GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
           as.numeric(sub(".*/", "", GeneRatio))) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

ego_down_top <- ego_down %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

ego_down_top <- ego_down_top %>%
  group_by(ONTOLOGY) %>%
  mutate(GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
           as.numeric(sub(".*/", "", GeneRatio))) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

#GO dotplot
go_up_dotplot <- ggplot(
  ego_up_top,
  aes(
    x = GeneRatio_num,
    y = Description,
    size = Count,
    color = p.adjust
  )
) +
  geom_point() +
  facet_grid(
    ONTOLOGY ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    name = "Adjusted p-value"
  ) +
  scale_size_continuous(
    name = "Count"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

go_up_dotplot

ggsave(file.path(output_dir, "GO_MDA_GT3_UP.png"), go_up_dotplot, width = 7, height = 4, dpi = 600)

go_down_dotplot <- ggplot(
  ego_down_top,
  aes(
    x = GeneRatio_num,
    y = Description,
    size = Count,
    color = p.adjust
  )
) +
  geom_point() +
  facet_grid(
    ONTOLOGY ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    name = "Adjusted p-value"
  ) +
  scale_size_continuous(
    name = "Count"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

go_down_dotplot

ggsave(file.path(output_dir, "GO_MDA_GT3_DOWN.png"), go_down_dotplot, width = 8, height = 4, dpi = 600)

