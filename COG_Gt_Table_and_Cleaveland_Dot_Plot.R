###Disclaimer###
#The following code assumes CSV outputs created using eggNOG-mapper.


library(tidyverse)
library(gt)
library(scales)

working_dir   <- "./" 
file_manifest <- c(
  "merged_eggnog_tax_2019.csv"   = "2019",
  "merged_eggnog_tax_2021.csv"   = "2021",
  "merged_eggnog_tax_hotbox.csv" = "Hotbox",
  "merged_eggnog_tax_cooler.csv" = "Cooler"
)


# STEP 1: COG Mapping Dictionary

cog_ref <- tibble(
  COG_letter = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","Z"),
  Broad_Class = c(
    "Information Storage & Processing", "Information Storage & Processing", "Metabolism",
    "Cellular Processes & Signaling",   "Metabolism",                       "Metabolism",
    "Metabolism",                       "Metabolism",                       "Metabolism",
    "Information Storage & Processing", "Information Storage & Processing", "Information Storage & Processing",
    "Cellular Processes & Signaling",   "Cellular Processes & Signaling",   "Cellular Processes & Signaling",
    "Metabolism",                       "Metabolism",                       "Poorly Characterized",
    "Poorly Characterized",             "Cellular Processes & Signaling",   "Cellular Processes & Signaling",
    "Cellular Processes & Signaling",   "Cellular Processes & Signaling",   "Cellular Processes & Signaling"
  ),
  COG_Name = c(
    "RNA processing and modification", "Chromatin structure and dynamics", "Energy production and conversion",
    "Cell cycle control, cell division, chromosome partitioning", "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism", "Carbohydrate transport and metabolism", "Coenzyme transport and metabolism",
    "Lipid transport and metabolism", "Translation, ribosomal structure and biogenesis", "Transcription",
    "Replication, recombination and repair", "Cell wall/membrane/envelope biogenesis", "Cell motility",
    "Post-translational modification, protein turnover, chaperones", "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism", "General function prediction only",
    "Function unknown", "Signal transduction mechanisms", "Intracellular trafficking, secretion, and vesicular transport",
    "Defense mechanisms", "Extracellular structures", "Cytoskeleton"
  )
)


# STEP 2: Parser Function (Tracks total contigs before duplicating annotations with multiple letters)

parse_eggnog_csv <- function(filepath, sample_name) {
  
  # Remove unannotated and missing contigs
  raw_filtered <- read_csv(filepath) %>%
    select(contig, COG_category) %>%
    filter(!is.na(COG_category) & COG_category != "-")
  
  # Total row count represents the denominator (all valid annotated contigs)
  total_contig_count <- nrow(raw_filtered)
  
  # Duplicate multi-letter annotations (e.g., "KL" creates an entry for K and L)
  counts <- raw_filtered %>%
    separate_rows(COG_category, sep = "") %>% 
    filter(COG_category != "") %>% 
    group_by(COG_category) %>%
    summarise(!!sym(sample_name) := n(), .groups = "drop")
  
  # Inject the denominator value into data framework for downstream math extraction
  bind_rows(
    counts,
    tibble(COG_category = "Total_Contigs", !!sym(sample_name) := total_contig_count)
  )
}


# STEP 3: Compile Master Data Frame via Full Outer Joins

cog_master_table <- cog_ref
active_samples   <- unname(file_manifest)

for (file_inst in names(file_manifest)) {
  full_path <- file.path(working_dir, file_inst)
  parsed_df <- parse_eggnog_csv(full_path, file_manifest[file_inst])
  
  # full_join ensures the customized "Total_Contigs" metric row isn't dropped
  cog_master_table <- full_join(cog_master_table, parsed_df, by = c("COG_letter" = "COG_category"))
}

# Replace any structurally unaligned NAs with zero
cog_master_table <- cog_master_table %>%
  mutate(across(all_of(active_samples), ~replace_na(., 0)))


# STEP 4: Build COG Category Gt Table

table_data <- cog_master_table %>%
  filter(!COG_letter %in% c("Total_Contigs")) %>%
  filter(Broad_Class != "Poorly Characterized")

cog_gt_table <- table_data %>%
  group_by(Broad_Class) %>%
  gt(rowname_col = "COG_letter") %>%
  tab_header(
    title = md("**Functional Profile of Metagenomic Assemblies by COG Category**"),
    subtitle = "Number of Annotations per Assembly"
  ) %>%
  tab_spanner(label = "Year", columns = c(`2019`, `2021`)) %>%
  tab_spanner(label = "Location", columns = c(Hotbox, Cooler)) %>%
  cols_label(
    COG_Name = "Category Description",
    `2019` = "2019", `2021` = "2021", Hotbox = "Hotbox", Cooler = "Cooler"
  ) %>%
  fmt_number(columns = c(`2019`, `2021`, Hotbox, Cooler), decimals = 0, use_seps = TRUE) %>%
  tab_options(
    row_group.font.weight = "bold", column_labels.font.weight = "bold",
    table.width = pct(100), data_row.padding = px(6)
  )

print(cog_gt_table)
gtsave(cog_gt_table, "COG_Category_Gt_Table.png", vwidth = 1000, vheight = 1500, zoom = 2)


# STEP 5: Reshape and Normalize (Proportion of Contigs Calculations)

# Extract true sample denominators from the data framework
denominators <- cog_master_table %>%
  filter(COG_letter == "Total_Contigs") %>%
  pivot_longer(cols = c(Hotbox, Cooler, `2019`, `2021`), names_to = "Condition", values_to = "Total_Contigs") %>%
  select(Condition, Total_Contigs)

# Format visual matrices and strip out poorly characterized categories
dot_data <- cog_master_table %>%
  filter(COG_letter != "Total_Contigs") %>% 
  pivot_longer(cols = c(Hotbox, Cooler, `2019`, `2021`), names_to = "Condition", values_to = "Counts") %>%
  left_join(denominators, by = "Condition") %>%
  group_by(Condition) %>%
  mutate(
    Percentage = (Counts / Total_Contigs) * 100,
    COG_Label = paste0(COG_Name, " (", COG_letter, ")")
  ) %>%
  ungroup() %>%
  filter(Broad_Class != "Poorly Characterized")


# STEP 6: Cleveland Dot Plot Configuration

cog_plot <- ggplot(dot_data, aes(x = Percentage, y = reorder(COG_Label, Percentage), color = Condition)) +
  geom_line(aes(group = COG_Label), color = "grey90") + 
  geom_point(size = 5) +
  scale_color_manual(
    values = c("Hotbox" = "#D55E00", "Cooler" = "#009E73", "2019" = "#56B4E9", "2021" = "#E69F00"),
    breaks = c("Hotbox", "Cooler", "2019", "2021")
  ) +
  scale_x_continuous(breaks = 1:10) + 
  facet_wrap(~Broad_Class, scales = "free_y", ncol = 1) +
  labs(
    title = "Relative Abundance of COG Functional Profiles Across Experimental Conditions",
    x = "Proportion of Annotated Contigs (%)", 
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(cog_plot)
ggsave("COG_Dot_Plot.png", plot = cog_plot, width = 11, height = 9, dpi = 600)