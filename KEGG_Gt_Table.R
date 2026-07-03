library(tidyverse)
library(gt)

# 1. Framework Setup
pathway_metadata <- tibble(
  Pathway_ID   = c("map02010", "map03070", "map02024", "map05111, map02025, map02026", "map02030", "map02040"),
  Pathway_Name = c("ABC Transporters", "Bacterial Secretion Systems", "Quorum Sensing", 
                   "Biofilm Formation", "Bacterial Chemotaxis", "Flagellar Assembly")
)

# 2. Direct Quantification Logic
get_pathway_counts <- function(date_file, loc_file) {
  df_loc  <- read.csv(loc_file, check.names = FALSE)  %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  df_date <- read.csv(date_file, check.names = FALSE) %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  
  tibble(
    Hotbox    = sum(df_loc$log2FoldChange < 0),
    Cooler    = sum(df_loc$log2FoldChange > 0),
    Year_2019 = sum(df_date$log2FoldChange < 0),
    Year_2021 = sum(df_date$log2FoldChange > 0)
  )
}

# 3. File Mapping & Matrix Compilation
pathway_mapping <- list(
  "map02010" = list(date = "abc_transporters_date.csv", loc = "abc_transporters_location.csv"),
  "map03070" = list(date = "bacterial_secretion_date.csv", loc = "bacterial_secretion_location.csv"),
  "map02024" = list(date = "quorum_sensing_date.csv", loc = "quorum_sensing_location.csv"),
  "map05111, map02025, map02026" = list(date = "biofilm_formation_date.csv", loc = "biofilm_formation_location.csv"),
  "map02030" = list(date = "bacterial_chemotaxis_date.csv", loc = "bacterial_chemotaxis_location.csv"),
  "map02040" = list(date = "flagellar_assembly_date.csv", loc = "flagellar_assembly_location.csv")
)

compiled_counts <- map2_df(pathway_mapping, names(pathway_mapping), ~ {
  get_pathway_counts(.x$date, .x$loc) %>% mutate(Pathway_ID = .y)
})

supplemental_table_1 <- pathway_metadata %>% 
  left_join(compiled_counts, by = "Pathway_ID") %>% 
  select(`Pathway ID` = Pathway_ID, `Pathway Name` = Pathway_Name, Hotbox, Cooler, `2019` = Year_2019, `2021` = Year_2021)

write.csv(supplemental_table_1, file = "Supplemental_Table_1_Generated.csv", row.names = FALSE)

# 4. GT Presentation Layout
gt_table <- supplemental_table_1 %>% 
  gt() %>% 
  tab_header(
    title = md("**Summary of Differentially Abundant Genes**"),
    subtitle = "Gene counts categorized by KEGG Pathway and experimental grouping"
  ) %>% 
  tab_spanner(label = md("**Location**"), columns = c(Hotbox, Cooler)) %>% 
  tab_spanner(label = md("**Timepoint**"), columns = c(`2019`, `2021`)) %>% 
  fmt_number(columns = c(Hotbox, Cooler, `2019`, `2021`), decimals = 0, use_seps = TRUE) %>% 
  cols_align(align = "left", columns = c(`Pathway ID`, `Pathway Name`)) %>% 
  cols_align(align = "right", columns = c(Hotbox, Cooler, `2019`, `2021`)) %>% 
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = `Pathway Name`)) %>% 
  tab_options(
    heading.align = "center", table.border.top.color = "white",
    table.border.bottom.color = "darkgrey", table_body.border.bottom.color = "darkgrey",
    heading.border.bottom.color = "lightgrey", column_labels.border.bottom.color = "darkgrey",
    column_labels.border.bottom.width = px(2), data_row.padding = px(6)
  )

print(gt_table)