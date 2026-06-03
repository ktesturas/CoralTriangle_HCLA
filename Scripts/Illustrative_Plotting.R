# Illustrative Plotting

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

script_name <- "Illustrative_Plotting"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)
sol_dir <- file.path(out_dir, "Solutions", script_name)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)


# Illustrative figure climate priority area appraoaches 

# 1) Plot for Adaptation Hotspot metric
# Load libraries
library(dplyr)
library(purrr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(cmocean)
library(here)
library(patchwork)


#===============================================================================
# Load data
#===============================================================================

rds_dir3 <- file.path(out_dir, "RDS", "03_Create_Features")

PUs <- readRDS(file.path(rds_dir3, "out_Coral.rds"))

landmass <- ne_countries(scale = "medium", returnclass = "sf")




combined_metric_ssp245 <- readRDS(
  file.path(out_dir, "RDS", "combined_metric_ssp245.rds")
)

# Prepare the metric data for SSP245
names(combined_metric_ssp245)[2] <- "metric"

# Make sure layers are in lon/lat
combined_metric_ssp245 <- st_transform(combined_metric_ssp245, 4326)
landmass <- st_transform(landmass, 4326)

#===============================================================================
# Define crop extent directly
#===============================================================================

plot_bbox <- st_bbox(
  c(
    xmin = 118.5,
    xmax = 132.0,
    ymin = -10.5,
    ymax = 1.5
  ),
  crs = st_crs(4326)
)

plot_box <- st_as_sfc(plot_bbox)

# Crop layers to plot extent
combined_metric_ssp245_crop <- st_crop(combined_metric_ssp245, plot_bbox)
landmass_crop <- st_crop(landmass, plot_bbox)

# Define breaks for the color scale
breaks <- c(0, 20, 40, 60, 80)

#===============================================================================
# Plot adaptation hotspots for SSP245
#===============================================================================

gg_combi245_AH <- ggplot() +
  geom_sf(
    data = combined_metric_ssp245_crop,
    aes(fill = metric),
    colour = NA,
    size = 0.1
  ) +
  scale_fill_cmocean(
    name = "matter",
    direction = 1,
    limits = c(0, 90),
    breaks = breaks,
    guide = guide_colourbar(
      title = "Score",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1,
      barheight = 20
    )
  ) +
  geom_sf(
    data = landmass_crop,
    fill = "grey80",
    color = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = c(118.5, 132.0),
    ylim = c(-10.5, 1.5),
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = c(120, 130)
  ) +
  scale_y_continuous(
    breaks = c(-10, -5, 0)
  ) +
  labs(
    fill = "Score"
  ) +
  ggtitle("Adaptation Hotspots") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.25, "cm"),
    
    plot.title = element_text(size = 27, hjust = 0.5),
    
    legend.position = "right",
    legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),
    legend.text = element_text(size = 22, hjust = 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(2.00, "cm"),
    legend.box.margin = margin(l = 10),
    
    strip.text = element_text(size = 27),
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Show plot
gg_combi245_AH

# Save plot
ggsave(
  plot = gg_combi245_AH,
  filename = file.path(fig_dir, "AH_245_illustrative.png"),
  width = 9,
  height = 5,
  dpi = 300
)




#===============================================================================
# 2) Plot for Climate Refugia metric
#===============================================================================

combined_metric_ssp245_Refugia <- readRDS(
  file.path(out_dir, "RDS", "combined_metric_ssp245_Refugia.rds")
)

# Prepare the metric data for SSP245 Refugia
names(combined_metric_ssp245_Refugia)[2] <- "metric"

# Make sure layers are in lon/lat
combined_metric_ssp245_Refugia <- st_transform(combined_metric_ssp245_Refugia, 4326)

# Crop layer to plot extent
combined_metric_ssp245_Refugia_crop <- st_crop(combined_metric_ssp245_Refugia, plot_bbox)

#===============================================================================
# Plot climate refugia for SSP245
#===============================================================================

gg_combi245_Refugia <- ggplot() +
  geom_sf(
    data = combined_metric_ssp245_Refugia_crop,
    aes(fill = metric),
    colour = NA,
    size = 0.1
  ) +
  scale_fill_cmocean(
    name = "deep",
    direction = 1,
    limits = c(0, 90),
    breaks = breaks,
    guide = guide_colourbar(
      title = "Score",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1,
      barheight = 20
    )
  ) +
  geom_sf(
    data = landmass_crop,
    fill = "grey80",
    color = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = c(118.5, 132.0),
    ylim = c(-10.5, 1.5),
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = c(120, 130)
  ) +
  scale_y_continuous(
    breaks = c(-10, -5, 0)
  ) +
  labs(
    fill = "Score"
  ) +
  ggtitle("Climate Refugia") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.25, "cm"),
    
    plot.title = element_text(size = 27, hjust = 0.5),
    
    legend.position = "right",
    legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),
    legend.text = element_text(size = 22, hjust = 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(2.00, "cm"),
    legend.box.margin = margin(l = 10),
    
    strip.text = element_text(size = 27),
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Show plot
gg_combi245_Refugia

# Save plot
ggsave(
  plot = gg_combi245_Refugia,
  filename = file.path(fig_dir, "Refugia_245_illustrative.png"),
  width = 9,
  height = 5,
  dpi = 300
)






#===============================================================================
# 3) Plot for Combined AH-Refugia metric
#===============================================================================

Combined_ssp245_AH_R <- readRDS(
  file.path(out_dir, "RDS", "Combined_ssp245_AH_R.rds")
)

# Prepare the metric data for SSP245 Combined AH-Refugia
names(Combined_ssp245_AH_R)[2] <- "metric"

# Make sure layer is in lon/lat
Combined_ssp245_AH_R <- st_transform(Combined_ssp245_AH_R, 4326)

# Crop layer to plot extent
Combined_ssp245_AH_R_crop <- st_crop(Combined_ssp245_AH_R, plot_bbox)

# Define breaks for the combined difference scale
breaks_combined <- c(-1, -0.5, 0, 0.5, 1)

##===============================================================================
# 3) Plot for Combined AH-Refugia metric
#===============================================================================

Combined_ssp245_AH_R <- readRDS(
  file.path(out_dir, "RDS", "Combined_ssp245_AH_R.rds")
)

# Make sure layer is in lon/lat
Combined_ssp245_AH_R <- st_transform(Combined_ssp245_AH_R, 4326)

# Use the normalized difference as the plotting metric
Combined_ssp245_AH_R <- Combined_ssp245_AH_R %>%
  mutate(metric = normalized_diff)

# Crop layer to plot extent
Combined_ssp245_AH_R_crop <- st_crop(Combined_ssp245_AH_R, plot_bbox)

# Define breaks for the combined difference scale
breaks_combined <- c(-1, -0.5, 0, 0.5, 1)

#===============================================================================
# Plot combined AH-Refugia metric for SSP245
#===============================================================================

gg_combi245_AH_R <- ggplot() +
  geom_sf(
    data = Combined_ssp245_AH_R_crop,
    aes(fill = metric),
    colour = NA,
    size = 0.1
  ) +
  scale_fill_cmocean(
    name = "balance",
    direction = 1,
    limits = c(-1, 1),
    breaks = breaks_combined,
    na.value = "grey80",
    guide = guide_colourbar(
      title = "Score",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1,
      barheight = 20
    )
  ) +
  geom_sf(
    data = landmass_crop,
    fill = "grey40",
    color = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = c(118.5, 132.0),
    ylim = c(-10.5, 1.5),
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = c(120, 130)
  ) +
  scale_y_continuous(
    breaks = c(-10, -5, 0)
  ) +
  labs(
    fill = "Score"
  ) +
  ggtitle("Combined AH-Refugia") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.25, "cm"),
    
    plot.title = element_text(size = 27, hjust = 0.5),
    
    legend.position = "right",
    legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),
    legend.text = element_text(size = 22, hjust = 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(2.00, "cm"),
    legend.box.margin = margin(l = 10),
    
    strip.text = element_text(size = 27),
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Show plot
gg_combi245_AH_R

# Save plot
ggsave(
  plot = gg_combi245_AH_R,
  filename = file.path(fig_dir, "Combined_AH_R_245_illustrative.png"),
  width = 9,
  height = 5,
  dpi = 300
)





#===============================================================================
# Define total threshold N
# Example: if N = 10, this selects top 5% and bottom 5%
#===============================================================================

N <- 10
tail_prop <- (N / 100) / 2

# Get lower and upper cutoffs
lower_cut <- quantile(combined_diff_ssp245$metric, probs = tail_prop, na.rm = TRUE)
upper_cut <- quantile(combined_diff_ssp245$metric, probs = 1 - tail_prop, na.rm = TRUE)

# Create binary selected / not selected column
combined_diff_ssp245_bin <- combined_diff_ssp245 %>%
  mutate(
    selected = ifelse(metric <= lower_cut | metric >= upper_cut, 1, 0)
  )

#===============================================================================
# Crop directly to the map window you want
#===============================================================================

plot_bbox <- st_bbox(
  c(
    xmin = 118.5,
    xmax = 132.0,
    ymin = -10.5,
    ymax = 1.5
  ),
  crs = st_crs(4326)
)

combined_diff_ssp245_bin_crop <- st_crop(combined_diff_ssp245_bin, plot_bbox)
landmass_crop <- st_crop(landmass, plot_bbox)

#===============================================================================
# Plot binary selected cells
#===============================================================================

gg_diff_binary <- ggplot() +
  
  geom_sf(
    data = landmass_crop,
    fill = "grey80",
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_sf(
    data = combined_diff_ssp245_bin_crop,
    aes(fill = factor(selected), color = factor(selected)),
    size = 0.1
  ) +
  
  scale_fill_manual(
    values = c(
      "0" = "#E6F2DA",
      "1" = "#D24757"
    ),
    name = "Cells",
    labels = c(
      "0" = "No",
      "1" = "Yes"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "0" = "#E6F2DA",
      "1" = "#D24757"
    ),
    name = "Cells",
    labels = c(
      "0" = "No",
      "1" = "Yes"
    )
  ) +
  
  geom_sf(
    data = landmass_crop,
    fill = NA,
    color = "grey70",
    linewidth = 0.4
  ) +
  
  coord_sf(
    xlim = c(118.5, 132.0),
    ylim = c(-10.5, 1.5),
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(120, 130)
  ) +
  
  scale_y_continuous(
    breaks = c(-10, -5, 0)
  ) +
  
  ggtitle("Selected cells") +
  
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.25, "cm"),
    
    plot.title = element_text(size = 22, hjust = 0.5),
    
    legend.position = "right",
    legend.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.9, "cm"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    fill = guide_legend(title.position = "top", title.hjust = 0.5),
    color = "none"
  )
# Show plot
gg_diff_binary

# Save plot
ggsave(
  plot = gg_diff_binary,
  filename = file.path(fig_dir, "Diff_top_bottom_binary_ssp245.png"),
  width = 6,
  height = 4.5,
  dpi = 300
)







#===============================================================================
# Top N/2 percentile maps for climate adaptation hotspots and climate refugia
#===============================================================================

N <- 10
top_prop <- (N / 100) / 2

plot_bbox <- st_bbox(
  c(
    xmin = 118.5,
    xmax = 132.0,
    ymin = -10.5,
    ymax = 1.5
  ),
  crs = st_crs(4326)
)

landmass <- st_transform(landmass, 4326)
landmass_crop <- st_crop(landmass, plot_bbox)


AH_cut <- quantile(
  combined_metric_ssp245$metric,
  probs = 1 - top_prop,
  na.rm = TRUE
)

R_cut <- quantile(
  combined_metric_ssp245_Refugia$metric,
  probs = 1 - top_prop,
  na.rm = TRUE
)

combined_metric_ssp245_AH_bin_crop <- combined_metric_ssp245 %>%
  mutate(selected = ifelse(metric >= AH_cut, 1, 0)) %>%
  st_crop(plot_bbox)

combined_metric_ssp245_R_bin_crop <- combined_metric_ssp245_Refugia %>%
  mutate(selected = ifelse(metric >= R_cut, 1, 0)) %>%
  st_crop(plot_bbox)

make_top_binary_map <- function(data, plot_title) {
  
  ggplot() +
    geom_sf(
      data = landmass_crop,
      fill = "grey80",
      color = "black",
      linewidth = 0.3
    ) +
    geom_sf(
      data = data,
      aes(fill = factor(selected), color = factor(selected)),
      size = 0.1
    ) +
    scale_fill_manual(
      values = c(
        "0" = "#E6F2DA",
        "1" = "#D24757"
      ),
      name = "Cells",
      labels = c(
        "0" = "No",
        "1" = "Yes"
      )
    ) +
    scale_color_manual(
      values = c(
        "0" = "#E6F2DA",
        "1" = "#D24757"
      ),
      name = "Cells",
      labels = c(
        "0" = "No",
        "1" = "Yes"
      )
    ) +
    geom_sf(
      data = landmass_crop,
      fill = NA,
      color = "grey70",
      linewidth = 0.4
    ) +
    coord_sf(
      xlim = c(118.5, 132.0),
      ylim = c(-10.5, 1.5),
      expand = FALSE
    ) +
    scale_x_continuous(
      breaks = c(120, 130)
    ) +
    scale_y_continuous(
      breaks = c(-10, -5, 0)
    ) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black", size = 18),
      axis.title = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.25, "cm"),
      
      plot.title = element_text(size = 22, hjust = 0.5),
      
      legend.position = "right",
      legend.title = element_text(size = 18, hjust = 0.5),
      legend.text = element_text(size = 16),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(0.9, "cm"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    guides(
      fill = guide_legend(title.position = "top", title.hjust = 0.5),
      color = "none"
    )
}

gg_AH_binary <- make_top_binary_map(
  combined_metric_ssp245_AH_bin_crop,
  "Hotspots"
)

gg_R_binary <- make_top_binary_map(
  combined_metric_ssp245_R_bin_crop,
  "Refugia"
)

gg_AH_binary
gg_R_binary

ggsave(
  plot = gg_AH_binary,
  filename = file.path(fig_dir, "AH_top_binary_ssp245.png"),
  width = 6,
  height = 4.5,
  dpi = 300
)

ggsave(
  plot = gg_R_binary,
  filename = file.path(fig_dir, "R_top_binary_ssp245.png"),
  width = 6,
  height = 4.5,
  dpi = 300
)




