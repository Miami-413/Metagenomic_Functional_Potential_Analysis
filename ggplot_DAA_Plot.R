# Load required libraries
library(tidyverse)
library(patchwork)
library(ggrepel)
# 1. Colorblind-Friendly Palette -------------------------------------------
color_map <- c(
  "DATE_UP"   = "#E69F00", "DATE_DOWN" = "#56B4E9",
  "LOC_UP"    = "#009E73", "LOC_DOWN"  = "#D55E00", "NS" = "#999999"
)
label_map <- c(
  "DATE_UP"   = "2021", "DATE_DOWN" = "2019",
  "LOC_UP"    = "Cooler", "LOC_DOWN"  = "Hotbox", "NS" = "Not significant"
)
# 2. Helper Function for Row Titles ---------------------------------------
make_row_title <- function(text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = text, fontface = "bold", size = 6) +
    theme_void() +
    theme(plot.margin = margin(t = 20, b = -10))
}
# 3. Updated Plotting Function (Strict 1 Label Per Direction) -------------
create_volcano <- function(filepath, subtitle, type = c("date", "location"),
                           x_limits = c(-25, 25), y_limits = c(0, 25),
                           x_breaks = seq(-25, 25, 5)) {
  
  type <- match.arg(type)
  
  df_base <- read.csv(filepath) %>%
    filter(!is.na(padj)) %>%
    mutate(
      # Create a unique ID for every single point to prevent duplicate naming
      RowID = row_number(),
      Enriched = case_when(
        type == "date" & log2FoldChange > 1  & padj < 0.05 ~ "DATE_UP",
        type == "date" & log2FoldChange < -1 & padj < 0.05 ~ "DATE_DOWN",
        type == "location" & log2FoldChange > 1  & padj < 0.05 ~ "LOC_UP",
        type == "location" & log2FoldChange < -1 & padj < 0.05 ~ "LOC_DOWN",
        TRUE ~ "NS"
      ),
      Enriched = factor(Enriched, levels = names(color_map))
    )
  
  # Identify the top upregulated gene by combined significance × magnitude score
  top_up_id <- df_base %>%
    filter(log2FoldChange > 0, padj < 0.05) %>%
    mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    pull(RowID)
  
  # Identify the top downregulated gene by combined significance × magnitude score
  top_down_id <- df_base %>%
    filter(log2FoldChange < 0, padj < 0.05) %>%
    mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    pull(RowID)
  
  # Only add the label if the point matches those specific Row IDs
  df <- df_base %>%
    mutate(plot_label = ifelse(RowID %in% c(top_up_id, top_down_id), Preferred_name, ""))
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Enriched)) +
    geom_vline(xintercept = c(-1, 1), color = "gray85", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "gray85", linetype = "dashed") +
    geom_point(alpha = 0.7, size = 1.3) +
    
    # REPEL SETTINGS FOR EXACTLY 2 LABELS
    geom_text_repel(
      aes(label = plot_label),
      size = 3.5,               
      fontface = "italic",
      box.padding = 1.5,        
      point.padding = 0.5,      
      min.segment.length = 0,   
      segment.color = "grey30",  
      segment.size = 0.5,
      max.overlaps = Inf,       
      force = 30,               
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = color_map, labels = label_map, drop = FALSE, name = "Enrichment") +
    scale_x_continuous(breaks = x_breaks) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    labs(subtitle = subtitle, x = expression(log[2]~"Fold Change"), y = expression("-log"[10]*"p-value"[adj])) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
      legend.position = "right",
      plot.margin = margin(20, 40, 20, 20)
    )
}
# 4. Generate Individual Plots --------------------------------------------
t1 <- make_row_title("ABC Transporters");
p1 <- create_volcano('ABC_transporters_location.csv', 'Hotbox vs. Cooler', type = "location");
p2 <- create_volcano('ABC_transporters_date.csv', '2019 vs. 2021', type = "date", x_limits = c(-10, 10), x_breaks = seq(-10, 10, 5), y_limits = c(0, 5))
t2 <- make_row_title("Bacterial Secretion Systems");
p3 <- create_volcano('bacterial_secretion_location.csv', 'Hotbox vs. Cooler', type = "location");
p4 <- create_volcano('bacterial_secretion_date.csv', '2019 vs. 2021', type = "date", x_limits = c(-10, 10), x_breaks = seq(-10, 10, 5), y_limits = c(0, 5))
t4 <- make_row_title("Quorum Sensing");
p7 <- create_volcano('quorum_sensing_location.csv', 'Hotbox vs. Cooler', type = "location");
p8 <- create_volcano('quorum_sensing_date.csv', '2019 vs. 2021', type = "date", x_limits = c(-10, 10), x_breaks = seq(-10, 10, 5), y_limits = c(0, 5))
t5 <- make_row_title("Biofilm Formation");
p9 <- create_volcano('biofilm_formation_location.csv', 'Hotbox vs. Cooler', type = "location");
p10 <- create_volcano('biofilm_formation_date.csv', '2019 vs. 2021', type = "date", x_limits = c(-10, 10), x_breaks = seq(-10, 10, 5), y_limits = c(0, 5))
t6 <- make_row_title("Bacterial Chemotaxis");
p11 <- create_volcano('bacterial_chemotaxis_location.csv', 'Hotbox vs. Cooler', type = "location");
p12 <- create_volcano('bacterial_chemotaxis_date.csv', '2019 vs. 2021', type = "date", x_limits = c(-10, 10), x_breaks = seq(-10, 10, 5), y_limits = c(0, 5))
t7 <- make_row_title("Flagellar Assembly");
p13 <- create_volcano('flagellar_assembly_location.csv', 'Hotbox vs. Cooler', type = "location");
p14 <- create_volcano('flagellar_assembly_date.csv', '2019 vs. 2021', type = "date", x_limits = c(-10, 10), x_breaks = seq(-10, 10, 5), y_limits = c(0, 5))
# 5. Unified Master Assembly (All 6 Included Pathways) --------------------
plot_master <- (
  t1 / (p1 | p2) /
    t2 / (p3 | p4) /
    t4 / (p7 | p8) /
    t5 / (p9 | p10) /
    t6 / (p11 | p12) /
    t7 / (p13 | p14)
) +
  plot_layout(guides = "collect", heights = rep(c(0.2, 1), 6)) +
  plot_annotation(
    title = "Differential Abundance Analysis - Selected KEGG Pathway Genes",
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )
# 6. Save Single Master File -----------------------------------------------
ggsave("Manuscript_Figure_Master_Combined.png", plot_master, width = 14, height = 26, dpi = 600)
message("Success! Check your working directory for Manuscript_Figure_Master_Combined.png")