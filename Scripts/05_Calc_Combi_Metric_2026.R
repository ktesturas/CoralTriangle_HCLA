# 05_Calc_Combi_Metrics: This R script calculates the combined resilience metrics 

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "05_Calc_Combi_Metrics"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)

# Set filepaths for RDS
rds_dir1 <- file.path(out_dir, "RDS", "01_Setup_PUs")  # to retrieve PUs

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Load prerequisites
# Load libraries
library(dplyr)
library(purrr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(cmocean)
library(patchwork)
library(cowplot)
library(viridis)
library(scales)
library(biscale)
library(grid)

# Load datasets
PUs <- readRDS(file = file.path(rds_dir1, "PUs.rds")) # PUs at whole planning region
landmass <- ne_countries(scale = "medium", returnclass = "sf") 

# Load Climate Metrics
scenarios <- c("ssp126", "ssp245", "ssp585")
variables <- c("tos", "o2os", "phos", "mhw")

# Create loop to store metrics in a list 
dfs_list <- list()

for (variable in variables) {
  for (scenario in scenarios) {
    df <- readRDS(file.path("Data", "Metrics", paste0(variable, "_", scenario, "_roc.rds"))) %>%
      sf::st_drop_geometry()
    colnames(df)[2] <- paste0(variable, "_", scenario) # rename second column to include the variable and scenario
    dfs_list[[paste0(variable, "_", scenario)]] <- df
  }
}

# Join everything as one df
metrics_tbl <- reduce(dfs_list, full_join, by = "cellID") %>%
  dplyr::filter(cellID %in% PUs$cellID)
metrics <- PUs %>%
  dplyr::left_join(metrics_tbl, by = "cellID")

# 2) Calculate Combined Metrics
# 2.1 Adaptation Hotspots
# Create function to calculate climate_smart metric
combined_metric_AH <- function(metrics, scenario, PUs){
  tos_col <- paste0("tos_", scenario)
  mhw_col <- paste0("mhw_", scenario)
  
  # Calculate combined metric
  combined <- metrics %>%
    mutate(
      !!paste0(tos_col, "_scaled") := scales::rescale(.data[[tos_col]], to = c(0, 100)), # the higher the metric values for tos, the more climate-smart
      !!paste0(mhw_col, "_scaled") := scales::rescale(-.data[[mhw_col]], to = c(0, 100)) # the lower the metric values for mhw, the more climate-smart
    ) %>%
    mutate(!!paste0("combined_", scenario) := (.data[[paste0(tos_col, "_scaled")]] * .data[[paste0(mhw_col, "_scaled")]])^(1/2)) %>% # geometric mean
    select(cellID, !!paste0("combined_", scenario), geometry) %>%
    sf::st_as_sf(crs = st_crs(PUs))
  
  return(combined)
}

# Loop through scenarios to calculate, save, and load 
combined_metrics <- lapply(scenarios, function(scenario) {
  combined_metric <- combined_metric_AH(metrics, scenario, PUs)
  file_path <- file.path(rds_dir, paste0("combined_metric_", scenario, ".rds"))
  saveRDS(combined_metric, file = file_path)
  readRDS(file_path)
})
names(combined_metrics) <- scenarios

# Access loaded objects
combined_metric_ssp126 <- combined_metrics[["ssp126"]]
combined_metric_ssp245 <- combined_metrics[["ssp245"]]
combined_metric_ssp585 <- combined_metrics[["ssp585"]]

# 2.2 Climate Refugia 
# Create function to calculate climate_smart metric
combined_metric_R <- function(metrics, scenario, PUs){
  tos_col <- paste0("tos_", scenario)
  mhw_col <- paste0("mhw_", scenario)
  
  # Calculate combined metric
  combined <- metrics %>%
    mutate(
      !!paste0(tos_col, "_scaled") := scales::rescale(-.data[[tos_col]], to = c(0, 100)), # the lower the metric values for tos, the more climate-smart
      !!paste0(mhw_col, "_scaled") := scales::rescale(-.data[[mhw_col]], to = c(0, 100)) # the lower the metric values for mhw, the more climate-smart
    ) %>%
    mutate(!!paste0("combined_", scenario) := (.data[[paste0(tos_col, "_scaled")]] * .data[[paste0(mhw_col, "_scaled")]])^(1/2)) %>% # geometric mean
    select(cellID, !!paste0("combined_", scenario), geometry) %>%
    sf::st_as_sf(crs = st_crs(PUs))
  
  return(combined)
}

# Loop through scenarios to calculate, save, and load 
combined_metrics_refugia <- lapply(scenarios, function(scenario) {
  combined_metric <- combined_metric_R(metrics, scenario, PUs)
  file_path <- file.path(rds_dir, paste0("combined_metric_", scenario, "_Refugia.rds"))
  saveRDS(combined_metric, file = file_path)
  readRDS(file_path)
})
names(combined_metrics_refugia) <- scenarios

# Access loaded objects
combined_metric_ssp126_refugia <- combined_metrics_refugia[["ssp126"]]
combined_metric_ssp245_refugia <- combined_metrics_refugia[["ssp245"]]
combined_metric_ssp585_refugia <- combined_metrics_refugia[["ssp585"]]


# 3) Plot resilience metric layers 
# 3.1 Adaptation Hotspots
# Create a function to plot adaptation hotspots 
fPlot_AH <- function(combined_metric){
  names(combined_metric)[2] <- "metric" 
  breaks <- c(0, 20, 40, 60, 80)
  
  gg_combi245_AH <- ggplot() +
    geom_sf(data = combined_metric, aes(fill = metric, color = metric), size = 0.1) +
    scale_fill_cmocean(
      name = "matter",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = c(0, 90),
      breaks = breaks,
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 1, barheight = 20)
    ) +
    scale_color_cmocean(
      name = "matter",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = c(0, 90),
      breaks = breaks,
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 1, barheight = 20)
    ) +
    geom_sf(data = landmass, fill = "grey80", color = "black") +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
    scale_x_continuous(breaks = seq(100, 160, by = 20)) +
    labs(fill = "Score",
         color = "Score") +
    ggtitle("Adaptation Hotspots") +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black", size = 18),
      axis.title = element_text(size = 27),
      plot.title = element_text(size = 27, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),
      legend.text = element_text(size = 22, hjust = 0.5),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(2.00, "cm"),
      legend.box.margin = margin(l = 10),
      strip.text = element_text(size = 27)
    ) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5),
           color = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  return(gg_combi245_AH)
}

# Save plots 
AH_SSP126 <- fPlot_AH(combined_metric_ssp126)
saveRDS(AH_SSP126, file = file.path(rds_dir, "AH_SSP126.rds"))
AH_SSP245 <- fPlot_AH(combined_metric_ssp245)
AH_SSP585 <- fPlot_AH(combined_metric_ssp585)

# Combine plots with spacing
combined_diff_plot_AH <- (AH_SSP126 + theme(plot.title = element_blank()) |
                            AH_SSP245 + theme(plot.title = element_blank()) |
                            AH_SSP585 + theme(plot.title = element_blank()) |
                            plot_spacer()) + 
  plot_layout(guides = "collect", widths = c(1, 1, 1, 0.05)) & 
  theme(legend.position = "right")

# Save the combined plot
ggsave(
  filename = file.path(fig_dir, "all_Hotspot_plots.png"),
  plot = combined_diff_plot_AH,
  width = 27,
  height = 15,
  dpi = 300
)


# 3.2 Climate Refugia
# Create a function to plot climate refugia 
fPlot_R <- function(combined_metric){
  names(combined_metric)[2] <- "metric" 
  breaks <- c(0, 20, 40, 60, 80)
  
  gg_combi245_R <- ggplot() +
    geom_sf(data = combined_metric, aes(fill = metric, color = metric), size = 0.1) +
    scale_fill_cmocean(
      name = "deep",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = c(0, 90),
      breaks = breaks,
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 1, barheight = 20)
    ) +
    scale_color_cmocean(
      name = "deep",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = c(0, 90),
      breaks = breaks,
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 1, barheight = 20)
    ) +
    geom_sf(data = landmass, fill = "grey80", color = "black") +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
    scale_x_continuous(breaks = seq(100, 160, by = 20)) +
    labs(fill = "Score",
         color = "Score") +
    ggtitle("Climate Refugia") +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black", size = 18),
      axis.title = element_text(size = 27),
      plot.title = element_text(size = 27, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),
      legend.text = element_text(size = 22, hjust = 0.5),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(2.00, "cm"),
      legend.box.margin = margin(l = 10),
      strip.text = element_text(size = 27)
    ) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5),
           color = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  return(gg_combi245_R)
}

# Save plots 
R_SSP126 <- fPlot_R(combined_metric_ssp126_refugia)
R_SSP245 <- fPlot_R(combined_metric_ssp245_refugia)
R_SSP585 <- fPlot_R(combined_metric_ssp585_refugia)

# Combine plots
combined_diff_plot_R <- (R_SSP126 + theme(plot.title = element_blank()) |
                           R_SSP245 + theme(plot.title = element_blank()) |
                           R_SSP585 + theme(plot.title = element_blank()) |
                           plot_spacer()) +
  plot_layout(guides = "collect", widths = c(1, 1, 1, 0.05)) & 
  theme(legend.position = "right")

# Save the combined plot
ggsave(
  filename = file.path(fig_dir, "all_Refugia_plots.png"),
  plot = combined_diff_plot_R,
  width = 27,
  height = 15,
  dpi = 300
)

# 4) Plot everything together
# Create function to modify plots
modify_plot <- function(plot) {
  plot + ggtitle(NULL) + theme(legend.position = "none")
}

modify_plot2 <- function(plot) {
  plot + ggtitle(NULL) + theme(legend.position = "right")
}

ah_plots <- list(AH_SSP126, AH_SSP245)
r_plots <- list(R_SSP126, R_SSP245)

ah_plots2 <- list(AH_SSP585)
r_plots2 <- list(R_SSP585)

# Apply function 
ah_plots <- lapply(ah_plots, modify_plot)
r_plots <- lapply(r_plots, modify_plot)

ah_plots2 <- lapply(ah_plots2, modify_plot2)
r_plots2 <- lapply(r_plots2, modify_plot2)

# Reassign modified plots
AH_SSP126 <- ah_plots[[1]]
AH_SSP245 <- ah_plots[[2]]
AH_SSP585 <- ah_plots2[[1]]
R_SSP126 <- r_plots[[1]]
R_SSP245 <- r_plots[[2]]
R_SSP585 <- r_plots2[[1]]

# Combine plots
combined_plot <- (AH_SSP126 | AH_SSP245 | AH_SSP585) /
  (R_SSP126 | R_SSP245 | R_SSP585)

ggsave(
  filename = file.path(fig_dir, "all_resilience_metrics.png"),
  plot = combined_plot,
  width = 27,
  height = 15,
  dpi = 300
)

# Plot SSP 2-4.5
gg_combined <- AH_SSP245 + R_SSP245 + 
  plot_layout(ncol = 2, nrow = 1)

ggsave(
  filename = file.path(fig_dir, "SSP245_resiience_metrics.png"),
  plot = gg_combined,
  width = 27,
  height = 15,
  dpi = 300
)

# 5) Calculate normalized differences 
# Create function to calculate and plot normalized differences
fCalc_and_plot_diff <- function(scenario, AH_data, R_data) {
  
  names(AH_data)[2] <- "metric" 
  names(R_data)[2] <- "metric" 
  
  # Calc the difference
  diff_data <- AH_data
  diff_data$combined_diff <- AH_data$metric - R_data$metric
  
  max_diff <- max(abs(diff_data$combined_diff), na.rm = TRUE)
  
  # Normalize the difference
  diff_data$normalized_diff <- diff_data$combined_diff / max_diff
  
  saveRDS(diff_data, file = file.path(rds_dir, paste0("Combined_", scenario, "_AH_R.rds")))
  
  # Create plot
  gg_diff <- ggplot() +
    geom_sf(data = diff_data, aes(fill = normalized_diff, color = normalized_diff), size = 0.1) +
    scale_fill_cmocean(
      name = "balance",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.5),
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 20, barheight = 1)
    ) +
    scale_color_cmocean(
      name = "balance",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.5),
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 20, barheight = 1)
    ) +
    geom_sf(data = landmass, fill = "grey40", color = "black") +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
    labs(fill = "Normalized Differences",
         color = "Normalized Differences") +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black", size = 10),
      legend.position = "right",
      legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),
      legend.text = element_text(size = 22, hjust = 0.5),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(2.00, "cm"),
      legend.box.margin = margin(l = 10),
      strip.text = element_text(size = 27)
    ) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5),
           color = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  return(gg_diff)
}

# Save plot 
Diff_126 <- fCalc_and_plot_diff("ssp126", combined_metric_ssp126, combined_metric_ssp126_refugia)
Diff_245 <- fCalc_and_plot_diff("ssp245", combined_metric_ssp245, combined_metric_ssp245_refugia)
Diff_585 <- fCalc_and_plot_diff("ssp585", combined_metric_ssp585, combined_metric_ssp585_refugia)

# Combine plots
combined_diff_plot <- (Diff_126 | Diff_245 | Diff_585) + 
  plot_layout(guides = "collect")

# Save the combined plot
ggsave(
  filename = file.path(fig_dir, "all_combined_difference_plots.png"),
  plot = combined_diff_plot,
  width = 27,
  height = 15,
  dpi = 300
)

# 6) Create Bivariate Plots 

# Define your function to calculate bivariate data and plot
fCalc_and_plot_bivariate <- function(scenario, AH_data, R_data, show_legend = TRUE) {
  names(AH_data)[2] <- "metric_AH"
  names(R_data)[2] <- "metric_R"
  
  bivariate_data <- st_join(AH_data, R_data)
  
  bivariate_data <- bivariate_data %>%
    mutate(
      norm_AH = scales::rescale(metric_AH),
      norm_R = scales::rescale(metric_R)
    )
  
  # Create the bivariate color scale
  bivariate_data <- bi_class(bivariate_data, x = norm_AH, y = norm_R, style = "quantile", dim = 4)
  
  map <- ggplot() +
    geom_sf(data = bivariate_data, aes(fill = bi_class), color = NA, size = 0.1) +
    bi_scale_fill(pal = "DkViolet2", dim = 4) +
    geom_sf(data = landmass, fill = "grey90", color = "black") +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
    scale_x_continuous(breaks = seq(100, 160, by = 20)) +
    bi_theme() +
    labs(fill = "Bivariate Choropleth") +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black", size = 18),
      axis.title = element_text(size = 27),
      plot.title = element_text(size = 27, hjust = 0.5),
      legend.position = "none",
      strip.text = element_text(size = 27)
    )
  
  if (show_legend) {
    legend <- bi_legend(
      pal = "DkViolet2",
      dim = 4,
      xlab = "Adaptation hotspot",
      ylab = "Climate refugia",
      size = 20
    )
    
    finalPlot <- ggdraw() +
      draw_plot(map, 0, 0, 1, 1) +
      draw_plot(legend, 0.1, 0.1, 0.4, 0.4)
    
    return(finalPlot)
  } else {
    return(map)
  }
}

# Generate the plot with legend for saving on its own
Biv_126_legend <- fCalc_and_plot_bivariate(
  "ssp126",
  combined_metric_ssp126,
  combined_metric_ssp126_refugia,
  show_legend = TRUE
)
ggsave(
  filename = file.path(fig_dir, "all_bivarlegend.png"),
  plot = Biv_126_legend,
  width = 27,
  height = 15,
  dpi = 300
)

# Generate plots without legend for combining
Biv_126 <- fCalc_and_plot_bivariate(
  "ssp126",
  combined_metric_ssp126,
  combined_metric_ssp126_refugia,
  show_legend = FALSE
)

Biv_245 <- fCalc_and_plot_bivariate(
  "ssp245",
  combined_metric_ssp245,
  combined_metric_ssp245_refugia,
  show_legend = FALSE
)

Biv_585 <- fCalc_and_plot_bivariate(
  "ssp585",
  combined_metric_ssp585,
  combined_metric_ssp585_refugia,
  show_legend = FALSE
)

# Combine plots
combined_bivariate_plot <- (Biv_126 + theme(plot.title = element_blank()) |
                              Biv_245 + theme(plot.title = element_blank()) |
                              Biv_585 + theme(plot.title = element_blank()))

# Save the combined plot
ggsave(
  filename = file.path(fig_dir, "all_combined_bivariate_plots.png"),
  plot = combined_bivariate_plot,
  width = 27,
  height = 15,
  dpi = 300
)
