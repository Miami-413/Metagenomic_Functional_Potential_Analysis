library(dplyr)
library(tidyr)
library(janitor)
library(gt)

keep_cols <- c("Predicted ARG Class")


# STEP 1: Parse and Combine Raw DeepARG Datasets

load_and_clean <- function(file_path, group_name) {
  read.csv(file_path, check.names = FALSE) %>%
    select(any_of(keep_cols)) %>%
    mutate(Group = group_name)
}

df_2019   <- load_and_clean("DeepARG_2019.csv", "2019")
df_2021   <- load_and_clean("DeepARG_2021.csv", "2021")
df_hotbox <- load_and_clean("DeepARG_hotbox.csv", "Hotbox")
df_cooler <- load_and_clean("DeepARG_cooler.csv", "Cooler")

all_data <- bind_rows(df_2019, df_2021, df_hotbox, df_cooler)


# STEP 2: Aggregate Counts and Format Summary Framework

summary_table <- all_data %>%
  group_by(`Predicted ARG Class`, Group) %>% 
  summarise(Gene_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = Gene_Count, values_fill = 0) %>%
  adorn_totals(where = c("row", "col")) %>%
  arrange(desc(Total))

write.csv(summary_table, "Full_ARG_Comparison_Table.csv", row.names = FALSE)


# STEP 3: Format and Export Publication-Ready Report Table (gt)

summary_table_csv <- read.csv("Full_ARG_Comparison_Table.csv", check.names = FALSE)

complete_table <- summary_table_csv %>%
  # Replace database underscores with spaces for presentation polish
  mutate(across(where(is.character), ~ gsub("_", " ", .))) %>%
  # Drop row total to allow gt to dynamically recalculate summary values
  filter(`Predicted ARG Class` != "Total") %>%
  select(`Predicted ARG Class`, `2019`, `2021`, Hotbox, Cooler, Total) %>%
  gt() %>%
  # Re-inject professional, dynamic bottom totals
  grand_summary_rows(
    columns = c(`2019`, `2021`, Hotbox, Cooler), 
    fns = list(label = "Total", id = "totals", fn = "sum"),
    fmt = list(~ fmt_number(., decimals = 0))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_grand_summary()
  ) %>%
  tab_header(
    title = md("**ARG Comparison Summary**")
  ) %>%
  # Assign footnotes to target index rows
  tab_footnote(
    footnote = "Multidrug ARG class comprises of various resistance mechanisms that do not belong to a single, specific drug class (e.g., efflux pump, ATP-binding cassette, etc.).",
    locations = cells_body(columns = `Predicted ARG Class`, rows = 1)
  ) %>%
  tab_footnote(
    footnote = "Macrolides, Lincosamides, and Streptogramins.",
    locations = cells_body(columns = `Predicted ARG Class`, rows = 3)
  ) %>%
  tab_footnote(
    footnote = "Unclassified genes represent putative resistance factors that do not yet have a confirmed specific antibiotic or antimicrobial target in the current database version.",
    locations = cells_body(columns = `Predicted ARG Class`, rows = 5)
  )

print(complete_table)
gtsave(complete_table, "Experimental_factors_deepARG_comparison_table.png", vheight = 2000)