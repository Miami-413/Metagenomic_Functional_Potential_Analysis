# ==============================================================================
# AMR METAGENOMIC ALLUVIAL VISUALIZATION PIPELINE
# ==============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(RColorBrewer)
library(tidyr)

# ------------------------------------------------------------------------------
# STEP 1: Load and Process AMRFinder Data
# ------------------------------------------------------------------------------
AMR_data <- read.csv("merged_beef_AMRFinder_Report_4_1_2026.csv", header = TRUE, stringsAsFactors = FALSE)

plot_data <- AMR_data %>%
  # Standardize textual data by trimming whitespace
  mutate(
    Genus = trimws(Genus),
    Subclass = trimws(Subclass)
  ) %>%
  # Aggregate raw counts for each Genus-Subclass pair
  group_by(Genus, Subclass) %>%
  summarise(Frequency = n(), .groups = "drop") 

# ------------------------------------------------------------------------------
# STEP 2: Map Custom Dynamic Color Palette
# ------------------------------------------------------------------------------
# Extract all distinct values across both variables for consistent mapping
unique_labels <- unique(c(plot_data$Genus, plot_data$Subclass))
getPalette    <- colorRampPalette(brewer.pal(12, "Paired"))
custom_colors <- getPalette(length(unique_labels))
names(custom_colors) <- unique_labels

# ------------------------------------------------------------------------------
# STEP 3: Configure and Output Alluvial Diagram
# ------------------------------------------------------------------------------
beef_alluvial_plot <- ggplot(data = plot_data,
                             aes(axis1 = Genus, axis2 = Subclass, y = Frequency)) +
  
  # Draw flow links shaded by the source organism
  geom_alluvium(aes(fill = Genus), width = 0.4, alpha = 0.4, knot.pos = 0.5) +
  
  # Draw categorical strata blocks colored by their specific identity
  geom_stratum(aes(fill = after_stat(stratum)), 
               width = 0.4, color = "grey40", linewidth = 0.3) +
  
  # Label blocks; applies italics only to column 1 (Genus names)
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum),
                fontface = ifelse(after_stat(x) == 1, "italic", "plain")), 
            size = 2.8) +
  
  # Scale and coordinate definitions
  scale_x_discrete(limits = c("Genus", "Resistance/Mechanism"), expand = c(.1, .1)) +
  scale_fill_manual(values = custom_colors) + 
  
  # Fine-tune visual layout 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12)
  ) +
  labs(title = "AMR Distribution by Bacterial Genus Across Metagenome",
       y = "Count of Hits")

print(beef_alluvial_plot)

# Export publication-ready image asset
ggsave("Beef_AMR_Alluvial.png", 
       plot = beef_alluvial_plot, 
       width = 8.5, 
       height = 12, 
       units = "in", 
       dpi = 300)