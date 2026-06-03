# 04_Plot_Metrics: This R script plots the climate metrics 

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "04_Plot_Metrics"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)

# Set filepaths for RDS
rds_dir1 <- file.path(out_dir, "RDS", "01_Setup_PUs")  # to retrieve PUs

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Load packages
pacman::p_load(tidyverse, sf, terra, here, purrr, rnaturalearth, parallel, patchwork, cmocean)
library(dplyr)

# Load datasets
PUs <- readRDS(file = file.path(rds_dir1, "PUs.rds"))
scenarios <- c("ssp126", "ssp245", "ssp585")
variables <- c("tos", "mhw")

for(variable in variables){
  for (scenario in scenarios){
    df <- readRDS(file = file.path("Data", "Metrics", paste0(variable, "_", scenario, "_roc.rds")))
    colnames(df)[2] <- paste0(names(df)[2], "_", scenario)
    assign(paste0(variable,"_",scenario), df, envir = .GlobalEnv)
  }
}

landmass <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") 

# --------------------------------------- #
# 1) Rate of Climate Warming 
# Create function for plotting 
fPlot_CW <- function(ssp_metric){
  
  names(ssp_metric)[2] <- "metric"
  
  gg_tos <- ggplot() +
    geom_sf(data = ssp_metric, aes(fill = metric, color = metric), size = 0.1) +
    scale_fill_cmocean(
      name = "thermal",
      direction = 1,
      aesthetics = c("fill", "color"),
      limits = range(ssp_metric$metric, na.rm = TRUE),
      labels = scales::label_number(accuracy = 0.001), 
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 20, barheight = 1)
    ) +
    geom_sf(data = landmass, fill = "grey80", color = "black") +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
    labs(
      fill = expression(''^"o"*'C yr'^"-1"*''),
      color = expression(''^"o"*'C yr'^"-1"*'')
    ) +
    ggtitle("Rate of Climate Warming") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 18),
          axis.title = element_text(size = 27),
          plot.title = element_text(size = 27, hjust = 0.5), # Center and enlarge the plot title
          legend.position = "bottom",  # Move the legend to the right
          legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 5)),  # Bold legend title and add margin
          legend.text = element_text(size = 22, hjust = 0.5),   # Center legend text
          #     legend.direction = "vertical",
          #     legend.box = "vertical",
          legend.key.width = unit(3.00, "cm"),
          legend.key.height = unit(0.50, "cm"),
          legend.box.margin = margin(l = 10),  # Move the legend box to the right
          strip.text = element_text(size = 27)) +  # Increase facet label text size
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5),
           color = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  return(gg_tos)
}

# Save plots 
CW_SSP126 <- fPlot_CW(tos_ssp126)
CW_SSP245 <- fPlot_CW(tos_ssp245)
CW_SSP585 <- fPlot_CW(tos_ssp585)

# 2) Sum of Cumulative MHW Intensities 
# Create function for plotting 
fPlot_MHW <- function(ssp_metric){
  
  names(ssp_metric)[2] <- "metric"
  
  gg_mhw <- ggplot() +
    geom_sf(data = ssp_metric, aes(fill = metric, color = metric), size = 0.1) +
    scale_fill_cmocean(
      name = "solar",
      direction = -1,
      aesthetics = c("fill", "color"),
      limits = range(ssp_metric$metric, na.rm = TRUE),
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 20, barheight = 1),
      labels = scales::comma  
    ) +
    geom_sf(data = landmass, fill = "grey80", color = "black") +
    coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
    labs(
      fill = "total degree days",
      color = "total degree days"
    ) +
    ggtitle("Sum of Cumulative MHW Intensities") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 18),
          axis.title = element_text(size = 27),
          plot.title = element_text(size = 27, hjust = 0.5), # Center and enlarge the plot title
          legend.position = "bottom",  # Move the legend to the right
          legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 5)),  # Bold legend title and add margin
          legend.text = element_text(size = 22, hjust = 0.5),   # Center legend text
          #     legend.direction = "vertical",
          #     legend.box = "vertical",
          legend.key.width = unit(3.00, "cm"),
          legend.key.height = unit(0.50, "cm"),
          legend.box.margin = margin(l = 10),  # Move the legend box to the right
          strip.text = element_text(size = 27)) +  # Increase facet label text size
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5),
           color = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  return(gg_mhw)
}

# Save plots
MHW_SSP126 <- fPlot_MHW(mhw_ssp126)
MHW_SSP245 <- fPlot_MHW(mhw_ssp245)
MHW_SSP585 <- fPlot_MHW(mhw_ssp585)

# --------------------------------------- #
# Create function to modify plots
modify_plot <- function(plot) {
  plot + ggtitle(NULL) + theme(legend.position = "bottom")
}

plots <- list(CW_SSP126, CW_SSP245, CW_SSP585, MHW_SSP126, MHW_SSP245, MHW_SSP585)
plots <- lapply(plots, modify_plot)

CW_SSP126 <- plots[[1]]
CW_SSP245 <- plots[[2]]
CW_SSP585 <- plots[[3]]
MHW_SSP126 <- plots[[4]]
MHW_SSP245 <- plots[[5]]
MHW_SSP585 <- plots[[6]]

# Plot everything together
combined_plot <- (CW_SSP126 | CW_SSP245 | CW_SSP585) /
  (MHW_SSP126 | MHW_SSP245 | MHW_SSP585)

ggsave(
  filename = file.path(fig_dir, "all_metrics.png"),
  plot = combined_plot,
  width = 27,
  height = 15,
  dpi = 300
)

# Plot SSP 2-4.5 
gg_combined <- CW_SSP245 + MHW_SSP245 + 
  plot_layout(ncol = 2, nrow = 1) 

ggsave(
  filename = file.path(fig_dir, "SSP245_metrics.png"),
  plot = gg_combined,
  width = 27,
  height = 15,
  dpi = 300
)
