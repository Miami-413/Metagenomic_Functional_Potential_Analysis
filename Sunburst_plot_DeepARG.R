library(plotly)
library(dplyr)
library(tidyr)
library(RColorBrewer)

sunburst_2019   <- read.csv("DeepARG_2019.csv", header = TRUE, sep = ",")
sunburst_2021   <- read.csv("DeepARG_2021.csv", header = TRUE, sep = ",")
sunburst_hotbox <- read.csv("DeepARG_hotbox.csv", header = TRUE, sep = ",")
sunburst_cooler <- read.csv("DeepARG_cooler.csv", header = TRUE, sep = ",")


# STEP 1: Data Preparation Function and Hierarchical Formatting

prep_sunburst_data <- function(df, label) {
  counts <- df %>%
    group_by(Predicted.ARG.Class) %>%
    summarise(Frequency = n(), .groups = "drop") %>%
    arrange(desc(Frequency))
  
  bind_rows(
    # Center Hub
    data.frame(ids = label, 
               labels = paste0("<b>", label, "</b>"), 
               parents = "", 
               values = sum(counts$Frequency)),
    # Outer Ring
    counts %>% mutate(ids = paste(label, Predicted.ARG.Class, sep = " - "), 
                      labels = Predicted.ARG.Class, 
                      parents = label, 
                      values = Frequency)
  )
}

data_2019   <- prep_sunburst_data(sunburst_2019, "2019")
data_2021   <- prep_sunburst_data(sunburst_2021, "2021")
data_hotbox <- prep_sunburst_data(sunburst_hotbox, "Hotbox")
data_cooler <- prep_sunburst_data(sunburst_cooler, "Cooler")


# STEP 2: Harmonize Global Color Map Across Subplots

all_labels <- unique(c(data_2019$labels, data_2021$labels, data_hotbox$labels, data_cooler$labels))
centers    <- c("<b>2019</b>", "<b>2021</b>", "<b>Hotbox</b>", "<b>Cooler</b>")
class_only <- all_labels[!all_labels %in% centers]

color_count <- length(class_only)
palette     <- colorRampPalette(brewer.pal(min(color_count, 12), "Set3"))(color_count)
color_map   <- setNames(palette, class_only)

# Assign a consistent neutral grey to center hub categories
for(node in centers) { color_map[node] <- "#D3D3D3" }


# STEP 3: Construct and Render 2x2 Layout Grid

fig_2x2 <- plot_ly() %>%
  # Top Left: 2019
  add_trace(
    ids = data_2019$ids, labels = data_2019$labels, parents = data_2019$parents, values = data_2019$values,
    type = 'sunburst', branchvalues = 'total',
    domain = list(x = c(0, 0.48), y = c(0.55, 1)),
    marker = list(colors = unname(color_map[data_2019$labels])),
    insidetextfont = list(size = 18, family = "Arial"),
    hovertemplate = "<b>%{label}</b><br>Count: %{value}<extra></extra>"
  ) %>%
  # Top Right: 2021
  add_trace(
    ids = data_2021$ids, labels = data_2021$labels, parents = data_2021$parents, values = data_2021$values,
    type = 'sunburst', branchvalues = 'total',
    domain = list(x = c(0.52, 1), y = c(0.55, 1)),
    marker = list(colors = unname(color_map[data_2021$labels])),
    insidetextfont = list(size = 18, family = "Arial"),
    hovertemplate = "<b>%{label}</b><br>Count: %{value}<extra></extra>"
  ) %>%
  # Bottom Left: Hotbox
  add_trace(
    ids = data_hotbox$ids, labels = data_hotbox$labels, parents = data_hotbox$parents, values = data_hotbox$values,
    type = 'sunburst', branchvalues = 'total',
    domain = list(x = c(0, 0.48), y = c(0, 0.45)),
    marker = list(colors = unname(color_map[data_hotbox$labels])),
    insidetextfont = list(size = 18, family = "Arial"),
    hovertemplate = "<b>%{label}</b><br>Count: %{value}<extra></extra>"
  ) %>%
  # Bottom Right: Cooler
  add_trace(
    ids = data_cooler$ids, labels = data_cooler$labels, parents = data_cooler$parents, values = data_cooler$values,
    type = 'sunburst', branchvalues = 'total',
    domain = list(x = c(0.52, 1), y = c(0, 0.45)),
    marker = list(colors = unname(color_map[data_cooler$labels])),
    insidetextfont = list(size = 18, family = "Arial"),
    hovertemplate = "<b>%{label}</b><br>Count: %{value}<extra></extra>"
  ) %>%
  layout(
    title = list(text = "<b>ARG Class Distribution Comparison</b>", y = 0.98, font = list(size = 30)),
    margin = list(l = 10, r = 10, b = 20, t = 100),
    uniformtext = list(minsize = 10, mode = 'hide')
  )

fig_2x2