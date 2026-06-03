# 07_SpatialPlans: This R script produces the spatial plans comparing climate refugia and adaptation hotspots planning approaches.

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "07_SpatialPlans"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)
sol_dir <- file.path(out_dir, "Solutions", script_name)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)

# Set filepaths for RDS
rds_dir1 <- file.path(out_dir, "RDS", "01_Setup_PUs")            # to retrieve PUs
rds_dir3 <- file.path(out_dir, "RDS", "03_Create_Features")     # to retrieve biodiversity data
rds_dir5 <- file.path(out_dir, "RDS", "05_Calc_Combi_Metrics")  # to retrieve combined metrics

# Load required libraries
# install.packages("devtools")
# devtools::install_github("https://github.com/SpatialPlanning/spatialplanr")
pacman::p_load(
  spatialplanr, prioritizr, dplyr, magrittr, stringr, sf, here, ggplot2, ggridges, forcats, purrr, patchwork, corrplot, Rsymphony)

# Load PUs and calculate median coral area
PUs <- readRDS(file.path(rds_dir3, "out_Coral.rds")) %>%
  dplyr::select(-any_of(c("EcoRegion", "Province", "Realm", "LatZone", "Area_MEOW")))
median_coral_area <- median(PUs$Coral_area_km2, na.rm = TRUE)

# Keep only PUs with coral areas >= median coral area
PUs <- PUs %>%
  filter(Coral_area_km2 >= median_coral_area)

# Load coral species occurrence data
out_spp_C <- readRDS(file.path(rds_dir3, "PUs.fin_spp_C.rds"))

# Load spatial prioritization functions
source("Utils/fClimatePriorityArea_CSapproach.R")        # CPA method

# Convert species occurrence probabilities to binary using 0.5 cutoff
out_spp_C_bin <- out_spp_C %>%
  as_tibble() %>%
  mutate(across(-any_of(c("cellID", "geometry")), ~ case_when(. >= 0.5 ~ 1, TRUE ~ 0))) %>%
  select(-which(colSums(select(., -geometry)) == 0)) %>%  # Drop species columns with no presence
  sf::st_as_sf()

# Keep only the filtered PUs and join additional PU info
out_spp_C_bin <- out_spp_C_bin %>%
  filter(cellID %in% PUs$cellID) %>%
  left_join(st_drop_geometry(PUs), by = "cellID") %>%
  select(-PU_area_km2, -Coral_area_km2) %>%
  st_as_sf()

# Load chronic (tos) and acute (mhw) climate data for each SSP scenario
scenarios <- c("ssp126", "ssp245", "ssp585")
variables <- c("tos", "mhw")

for(variable in variables){
  for (scenario in scenarios){
    df <- readRDS(file = paste0("Data/Metrics/", variable, "_", scenario, "_roc.rds"))
    colnames(df)[2] <- paste0(names(df)[2])  # Standardize second column name
    assign(paste0(variable,"_",scenario), df, envir = .GlobalEnv)  # Save to global environment
  }
}

# Load combined adaptation hotspot, refugia, and combined metrics for each scenario
for (scenario in scenarios) {
  combined_ah <- readRDS(file = paste0(rds_dir5,"/combined_metric_", scenario, ".rds")) %>%
    dplyr::rename(metric = names(.)[2])
  assign(paste0("Combined_", scenario, "_AH"), combined_ah, envir = .GlobalEnv)
  
  combined_r <- readRDS(file = paste0(rds_dir5,"/combined_metric_", scenario, "_Refugia.rds")) %>%
    dplyr::rename(metric = names(.)[2])
  assign(paste0("Combined_", scenario, "_R"), combined_r, envir = .GlobalEnv)
  
  combined_norm <- readRDS(file = paste0(rds_dir5,"/Combined_", scenario, "_AH_R.rds")) %>%
    dplyr::select(cellID, normalized_diff, geometry) %>%
    dplyr::rename(metric = "normalized_diff")
  assign(paste0("Combined_", scenario, "_norm"), combined_norm, envir = .GlobalEnv)
}

# Load world landmass for plotting
landmass <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Define plotting theme
splnr_theme <- list(
  ggplot2::theme_bw(),
  ggplot2::theme(
    legend.position = "right",
    legend.direction = "vertical",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    plot.title = ggplot2::element_text(size = 9),
    axis.title = ggplot2::element_blank()
  )
)

# Create species targets: 30% of each feature's distribution
C_target_df <- out_spp_C_bin %>%
  as_tibble() %>%
  select(-geometry, -cellID) %>%
  names() %>%
  tibble::enframe(name = NULL, value = "feature") %>%
  mutate(target = 0.3)

# Prepare cost layer
UniformCost <- PUs %>%
  select(-Coral_area_km2) %>%
  rename(cost = PU_area_km2)

# Function to merge climate data, cost into spatial layer
prepare_out_sf <- function(CPA, scenario) {
  tos_var <- get(paste0("tos_ssp", scenario))
  mhw_var <- get(paste0("mhw_ssp", scenario))
  combined_ah_var <- get(paste0("Combined_ssp", scenario, "_AH"))
  combined_r_var <- get(paste0("Combined_ssp", scenario, "_R"))
  combined_norm_var <- get(paste0("Combined_ssp", scenario, "_norm"))
  
  CPA %>%
    st_drop_geometry() %>%
    left_join(UniformCost %>% st_drop_geometry(), by = "cellID") %>%
    left_join(tos_var %>% st_drop_geometry(), by = "cellID") %>%
    left_join(mhw_var %>% st_drop_geometry(), by = "cellID") %>%
    left_join(combined_ah_var %>% st_drop_geometry(), by = "cellID") %>%
    rename(metric_AH = metric) %>%
    left_join(combined_r_var %>% st_drop_geometry(), by = "cellID") %>%
    rename(metric_R = metric) %>%
    left_join(combined_norm_var %>% st_drop_geometry(), by = "cellID") %>%
    rename(metric_norm = metric) %>%
    st_as_sf(geometry = CPA$geometry)
}

# Set up prioritizr spatial problem
prepare_problem <- function(out_sf, targets) {
  prioritizr::problem(out_sf, features = targets$feature, cost_column = "cost") %>%
    add_min_set_objective() %>%
    add_relative_targets(targets$target) %>%
    add_binary_decisions() %>%
    add_rsymphony_solver(gap = 0.1, time_limit = 600, verbose = TRUE)
}

# Solve prioritization problem
solve_problem <- function(problem, approach, file_name) {
  solution <- prioritizr::solve.ConservationProblem(problem, force = TRUE) %>%
    mutate(approach = approach)
  saveRDS(solution, file = file_name)
  solution
}

# Plot solution using custom plotting functions
plot_solution <- function(solution) {
  splnr_plot_solution(solution) +
    splnr_gg_add(overlay = landmass, cropOverlay = PUs, ggtheme = splnr_theme)
}

# Calculate % each species distribution protected in the solution
calculate_percent_protection <- function(data, name) {
  sp_cols_with_suffix <- colnames(data)[grepl("_CS|_NCS", colnames(data))]
  sp_cols <- unique(sub("_CS|_NCS", "", sp_cols_with_suffix))
  sp_cols <- unique(c(sp_cols, colnames(data)[grepl("^C_", colnames(data)) & !grepl("_CS|_NCS", colnames(data))]))
  protect_r <- data.frame(Species = character(), PercentProtection = numeric(), stringsAsFactors = FALSE)
  
  for(sp in sp_cols) {
    if(paste0(sp, "_CS") %in% colnames(data) & paste0(sp, "_NCS") %in% colnames(data)) {
      cs_col <- paste0(sp, "_CS")
      ncs_col <- paste0(sp, "_NCS")
      
      total <- sum(data[[cs_col]], data[[ncs_col]], na.rm = TRUE)
      protected <- sum((data[[cs_col]] + data[[ncs_col]]) * data$solution_1, na.rm = TRUE)
    } else {
      total <- sum(data[[sp]], na.rm = TRUE)
      protected <- sum(data[[sp]] * data$solution_1, na.rm = TRUE)
    }
    
    percent <- (protected / total) * 100
    protect_r <- rbind(protect_r, data.frame(Species = sp, PercentProtection = percent))
  }
  
  protect_r <- protect_r %>%
    dplyr::mutate(name = name)
  return(protect_r)
}

# Clean targets by removing features that are missing from out_sf or contain only 0/NA values
clean_targets <- function(out_sf, targets, id_col = "cellID", verbose = TRUE) {
  # checks
  if (!inherits(out_sf, "sf")) {
    stop("out_sf must be an sf object.")
  }
  
  if (!all(c("feature", "target") %in% names(targets))) {
    stop("targets must have columns named 'feature' and 'target'.")
  }
  
  # requested target features
  target_features <- targets$feature
  
  # features that exist in out_sf
  existing_features <- intersect(target_features, names(out_sf))
  
  # features missing from out_sf
  missing_features <- setdiff(target_features, names(out_sf))
  
  # subset existing feature columns only
  feat_df <- sf::st_drop_geometry(out_sf)[, existing_features, drop = FALSE]
  
  # features that are entirely 0 / NA
  zero_features <- names(feat_df)[
    vapply(feat_df, function(x) all(is.na(x) | x == 0), logical(1))
  ]
  
  # all bad features
  bad_features <- union(missing_features, zero_features)
  
  # cleaned targets
  targets_fixed <- targets %>%
    dplyr::filter(!feature %in% bad_features)
  
  if (verbose) {
    cat("Original target features:", nrow(targets), "\n")
    cat("Missing features removed:", length(missing_features), "\n")
    if (length(missing_features) > 0) print(missing_features)
    
    cat("Zero/NA-only features removed:", length(zero_features), "\n")
    if (length(zero_features) > 0) print(zero_features)
    
    cat("Remaining target features:", nrow(targets_fixed), "\n")
  }
  
  if (nrow(targets_fixed) == 0) {
    stop("All target features were removed by clean_targets().")
  }
  
  return(list(
    targets_fixed = targets_fixed,
    missing_features = missing_features,
    zero_features = zero_features,
    bad_features = bad_features
  ))
}


#===============================================================================
# LOCATOR MAP ONLY
#===============================================================================

# Load archived solution used for manuscript figures
archived_solution_path <- file.path("Data", "Solutions", "s-CPA-10-245.rds")
solution_AH <- readRDS(archived_solution_path)

#PUs <- readRDS(file = file.path(rds_dir1, "PUs.rds"))

# Make sure PUs matches the CRS of the solution / land layers
PUs <- st_transform(PUs, st_crs(solution_AH))

# Define tighter full-region bounding box in lon/lat
bbox_full <- st_as_sfc(
  st_bbox(c(xmin = 98, xmax = 166, ymin = -12, ymax = 19), crs = 4326)
)

# Transform bounding box to the CRS of the plotting layers
bbox_full <- st_transform(bbox_full, st_crs(solution_AH))
bbox_full <- st_bbox(bbox_full)

# Add breathing room around the map
xpad <- (bbox_full["xmax"] - bbox_full["xmin"]) * 0.05
ypad <- (bbox_full["ymax"] - bbox_full["ymin"]) * 0.05

# Create locator boxes for the three panel windows
window_df <- data.frame(
  xmin = c(113, 119, 147),
  xmax = c(126, 132, 160),
  ymin = c(5, -8, -11.5),
  ymax = c(14, 1, -2.5)
)

make_bbox_polygon <- function(xmin, xmax, ymin, ymax) {
  st_polygon(
    list(
      matrix(
        c(
          xmin, ymin,
          xmax, ymin,
          xmax, ymax,
          xmin, ymax,
          xmin, ymin
        ),
        ncol = 2,
        byrow = TRUE
      )
    )
  )
}

locator_boxes <- st_sf(
  geometry = st_sfc(
    mapply(
      make_bbox_polygon,
      xmin = window_df$xmin,
      xmax = window_df$xmax,
      ymin = window_df$ymin,
      ymax = window_df$ymax,
      SIMPLIFY = FALSE
    ),
    crs = 4326
  )
)

locator_boxes <- st_transform(locator_boxes, st_crs(solution_AH))

Locator_Map_Full_Region <- ggplot() +
  geom_sf(
    data = PUs,
    fill = "grey25",
    colour = NA,
    alpha = 0.9
  ) +
  geom_sf(
    data = landmass,
    colour = "grey55",
    fill = "grey90",
    alpha = 1,
    linewidth = 0.1
  ) +
  geom_sf(
    data = locator_boxes,
    fill = NA,
    colour = "blue",
    linewidth = 0.35,
    linetype = "solid",
    inherit.aes = FALSE
  ) +
  coord_sf(
    xlim = c(bbox_full["xmin"] - xpad, bbox_full["xmax"] + xpad),
    ylim = c(bbox_full["ymin"] - ypad, bbox_full["ymax"] + ypad),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text    = element_text(color = "black", size = 10),
    plot.title   = element_text(size = 16)
  ) +
  labs(
    title = "Locator Map",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(
  filename = file.path(fig_dir, "Locator_Map_Full_Region.png"),
  plot     = Locator_Map_Full_Region,
  width    = 9,
  height   = 4.5,
  dpi      = 300
)


#===============================================================================
# A) ADAPTATION HOTSPOTS
# Run SpatPrio for SSP 2-4.5 at 10% ONLY
#===============================================================================

scenario <- "245"
percentile <- 10

# Load adaptation hotspot metric
metricDF <- get(paste0("Combined_ssp", scenario, "_AH"))

# Run Climate Priority Area approach
CPA <- fClimatePriorityArea_CSapproach_bin(
  featuresDF = out_spp_C_bin,
  percentile = percentile,
  metricDF   = metricDF,
  direction  = 1
)

# Prepare spatial data and targets
out_sf <- prepare_out_sf(CPA, scenario)

targets <- fAssignTargets_CPA(
  climatesmartDF = CPA,
  targetsDF      = C_target_df,
  refugiaTarget  = 1
)

target_check  <- clean_targets(out_sf, targets)
targets_fixed <- target_check$targets_fixed

# Build and solve spatial prioritisation problem
p <- prepare_problem(out_sf, targets_fixed)

set.seed(123)

solution_AH <- solve_problem(
  p,
  paste0("CPA-", percentile),
  file.path(sol_dir, paste0("s-CPA-", percentile, "-", scenario, ".rds"))
)

#===============================================================================
# QUICK FULL-REGION SOLUTION PLOT — ADAPTATION HOTSPOTS
#===============================================================================

# Load archived solution used for manuscript figures.
# This avoids re-solving the prioritisation problem when only regenerating HD plots.
archived_solution_path <- file.path("Data", "Solutions", "s-CPA-10-245.rds")
solution_AH <- readRDS(archived_solution_path)

# Define tighter full-region bounding box in lon/lat
bbox_full <- st_as_sfc(
  st_bbox(c(xmin = 98, xmax = 166, ymin = -12, ymax = 19), crs = 4326)
)

# Transform bounding box to the CRS of the solution
bbox_full <- st_transform(bbox_full, st_crs(solution_AH))
bbox_full <- st_bbox(bbox_full)

# Add breathing room around the map
xpad <- (bbox_full["xmax"] - bbox_full["xmin"]) * 0.05
ypad <- (bbox_full["ymax"] - bbox_full["ymin"]) * 0.05

Hotspots_Full_Region_245_10 <- ggplot() +
  geom_sf(
    data = solution_AH,
    aes(fill = as.factor(solution_1)),
    colour = NA,
    size = 0.1
  ) +
  geom_sf(
    data = landmass,
    colour = "grey50",
    fill = "grey90",
    alpha = 0.9,
    size = 0.1
  ) +
  coord_sf(
    xlim = c(bbox_full["xmin"] - xpad, bbox_full["xmax"] + xpad),
    ylim = c(bbox_full["ymin"] - ypad, bbox_full["ymax"] + ypad),
    expand = FALSE
  ) +
  scale_fill_manual(
    values = c("1" = "#295E11", "0" = "#FED674"),
    breaks = c("0", "1"),
    labels = c("Not Selected", "Selected"),
    name   = "Planning Units"
  ) +
  theme_bw() +
  theme(
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text    = element_text(color = "black", size = 10),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 16)
  ) +
  labs(
    title = "Climate Adaptation Hotspots",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(
  filename = file.path(fig_dir, "Hotspots_Full_Region_245_10.png"),
  plot     = Hotspots_Full_Region_245_10,
  width    = 9,
  height   = 4.5,
  dpi      = 300
)
#===============================================================================
# REGIONAL PLOTS
#===============================================================================
# Helper function to crop planning units using bounding box
crop_pus <- function(bbox_values) {
  PUs_region <- readRDS(file = file.path(rds_dir1, "PUs.rds"))
  
  bbox <- st_as_sfc(st_bbox(bbox_values, crs = 4326))
  bbox <- st_transform(bbox, st_crs(PUs_region))
  
  PUs_region %>%
    filter(st_intersects(geometry, bbox, sparse = FALSE))
}


# Helper function to plot regional solution
plot_region <- function(PUs_region, y_breaks, x_breaks = waiver()) {
  ggplot() +
    geom_sf(
      data = solution_AH,
      aes(fill = as.factor(solution_1)),
      colour = NA,
      size = 0.1,
      show.legend = FALSE
    ) +
    geom_sf(
      data = landmass,
      colour = "grey30",
      fill = "grey90",
      alpha = 0.9,
      size = 0.1,
      show.legend = FALSE
    ) +
    coord_sf(
      xlim = c(st_bbox(PUs_region)$xmin, st_bbox(PUs_region)$xmax),
      ylim = c(st_bbox(PUs_region)$ymin, st_bbox(PUs_region)$ymax)
    ) +
    scale_fill_manual(
      values = c("1" = "#295E11", "0" = "#FED674"),
      breaks = c("0", "1")
    ) +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_continuous(breaks = x_breaks) +
    theme_bw() +
    theme(
      axis.ticks   = element_line(color = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.text    = element_text(color = "black", size = 14),
      axis.title   = element_blank(),
      plot.title   = element_blank()
    )
}

# Crop planning units
PUs_ph <- crop_pus(c(xmin = 113, xmax = 126, ymin = 5, ymax = 14))
PUs_west <- crop_pus(c(xmin = 119, xmax = 132, ymin = -8, ymax = 1))
PUs_east <- crop_pus(c(xmin = 147, xmax = 160, ymin = -11.5, ymax = -2.5))

# Create regional plots
plot_ph <- plot_region(
  PUs_region = PUs_ph,
  y_breaks   = c(5, 10, 15, 20)
)

plot_west <- plot_region(
  PUs_region = PUs_west,
  y_breaks   = seq(-10, 5, by = 5),
  x_breaks   = seq(95, 130, by = 10)
)

plot_east <- plot_region(
  PUs_region = PUs_east,
  y_breaks   = seq(-10, 5, by = 5),
  x_breaks   = seq(130, 165, by = 10)
)

# Combine and save plots
combined_plot_PH <- (
  plot_ph +
    theme(plot.margin = margin(0, 0, 10, 0))
) / (
  plot_west +
    theme(plot.margin = margin(10, 0, 10, 0))
) / (
  plot_east +
    theme(plot.margin = margin(10, 0, 0, 0))
)

ggsave(
  filename = file.path(fig_dir, "Hotspots_Plan_245_10perc_PH.png"),
  plot     = combined_plot_PH,
  width    = 6.8,
  height   = 14,
  dpi      = 300
)

#===============================================================================
# B) CLIMATE REFUGIA
# Run SpatPrio for SSP 2-4.5 at 10% ONLY
#===============================================================================

scenario <- "245"
percentile <- 10

# Load climate refugia metric
metricDF <- get(paste0("Combined_ssp", scenario, "_R"))

# Run Climate Priority Area approach
CPA <- fClimatePriorityArea_CSapproach_bin(
  featuresDF = out_spp_C_bin,
  percentile = percentile,
  metricDF   = metricDF,
  direction  = 1
)

# Prepare spatial data and targets
out_sf <- prepare_out_sf(CPA, scenario)

targets <- fAssignTargets_CPA(
  climatesmartDF = CPA,
  targetsDF      = C_target_df,
  refugiaTarget  = 1
)

target_check  <- clean_targets(out_sf, targets)
targets_fixed <- target_check$targets_fixed

# Build and solve spatial prioritisation problem
p <- prepare_problem(out_sf, targets_fixed)

set.seed(123)

solution_AH <- solve_problem(
  p,
  paste0("CPA-", percentile),
  file.path(sol_dir, paste0("s-CPA-", percentile, "-", scenario, "-Refugia.rds"))
)

#===============================================================================
# QUICK FULL-REGION SOLUTION PLOT — CLIMATE REFUGIA
#===============================================================================

# Load archived solution used for manuscript figures.
# This avoids re-solving the prioritisation problem when only regenerating HD plots.
archived_solution_path_Refugia <- file.path("Data", "Solutions", "s-CPA-10-245-Refugia.rds")
solution_R <- readRDS(archived_solution_path_Refugia)

# Define tighter full-region bounding box in lon/lat
bbox_full <- st_as_sfc(
  st_bbox(c(xmin = 98, xmax = 166, ymin = -12, ymax = 19), crs = 4326)
)

# Transform bounding box to the CRS of the solution
bbox_full <- st_transform(bbox_full, st_crs(solution_R))
bbox_full <- st_bbox(bbox_full)

# Add breathing room around the map
xpad <- (bbox_full["xmax"] - bbox_full["xmin"]) * 0.05
ypad <- (bbox_full["ymax"] - bbox_full["ymin"]) * 0.05

Refugia_Full_Region_245_10 <- ggplot() +
  geom_sf(
    data = solution_R,
    aes(fill = as.factor(solution_1)),
    colour = NA,
    size = 0.1
  ) +
  geom_sf(
    data = landmass,
    colour = "grey50",
    fill = "grey90",
    alpha = 0.9,
    size = 0.1
  ) +
  coord_sf(
    xlim = c(bbox_full["xmin"] - xpad, bbox_full["xmax"] + xpad),
    ylim = c(bbox_full["ymin"] - ypad, bbox_full["ymax"] + ypad),
    expand = FALSE
  ) +
  scale_fill_manual(
    values = c("1" = "#295E11", "0" = "#FED674"),
    breaks = c("0", "1"),
    labels = c("Not Selected", "Selected"),
    name   = "Planning Units"
  ) +
  theme_bw() +
  theme(
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text    = element_text(color = "black", size = 10),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 16)
  ) +
  labs(
    title = "Climate refugia",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(
  filename = file.path(fig_dir, "Refugia_Full_Region_245_10.png"),
  plot     = Refugia_Full_Region_245_10,
  width    = 9,
  height   = 4.5,
  dpi      = 300
)

#===============================================================================
# REGIONAL PLOTS
#===============================================================================

# Helper function to crop planning units using bounding box
crop_pus <- function(bbox_values) {
  PUs_region <- readRDS(file = file.path(rds_dir1, "PUs.rds"))
  
  bbox <- st_as_sfc(st_bbox(bbox_values, crs = 4326))
  bbox <- st_transform(bbox, st_crs(PUs_region))
  
  PUs_region %>%
    filter(st_intersects(geometry, bbox, sparse = FALSE))
}


# Helper function to plot regional solution
plot_region <- function(PUs_region, y_breaks, x_breaks = waiver()) {
  ggplot() +
    geom_sf(data = solution_R, aes(fill = as.factor(solution_1)), colour = NA, size = 0.1, show.legend = FALSE) +
    geom_sf(data = landmass, colour = "grey30", fill = "grey90", alpha = 0.9, size = 0.1, show.legend = FALSE) +
    coord_sf(
      xlim = c(st_bbox(PUs_region)$xmin, st_bbox(PUs_region)$xmax),
      ylim = c(st_bbox(PUs_region)$ymin, st_bbox(PUs_region)$ymax)
    ) +
    scale_fill_manual(values = c("1" = "#295E11", "0" = "#FED674")) +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_continuous(breaks = x_breaks) +
    theme_bw() +
    theme(
      axis.ticks   = element_line(color = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.text    = element_text(color = "black", size = 14),
      axis.title   = element_blank(),
      plot.title   = element_blank()
    )
}


# Crop planning units
PUs_ph   <- crop_pus(c(xmin = 113, xmax = 126, ymin = 5, ymax = 14))
PUs_west <- crop_pus(c(xmin = 119, xmax = 132, ymin = -8, ymax = 1))
PUs_east <- crop_pus(c(xmin = 147, xmax = 160, ymin = -11.5, ymax = -2.5))


# Create regional plots
plot_ph <- plot_region(
  PUs_region = PUs_ph,
  y_breaks   = c(5, 10, 15, 20)
)

plot_west <- plot_region(
  PUs_region = PUs_west,
  y_breaks   = seq(-10, 5, by = 5),
  x_breaks   = seq(95, 130, by = 10)
)

plot_east <- plot_region(
  PUs_region = PUs_east,
  y_breaks   = seq(-10, 5, by = 5),
  x_breaks   = seq(130, 165, by = 10)
)


# Combine and save plots
combined_plot_R <- (
  plot_ph +
    theme(plot.margin = margin(0, 0, 10, 0))
) / (
  plot_west +
    theme(plot.margin = margin(10, 0, 10, 0))
) / (
  plot_east +
    theme(plot.margin = margin(10, 0, 0, 0))
)

ggsave(
  filename = file.path(fig_dir, "Refugia_Plan_245_10perc_PH.png"),
  plot     = combined_plot_R,
  width    = 6.8,
  height   = 14,
  dpi      = 300
)


#===============================================================================
# C) CLIMATE-UNINFORMED
# Run SpatPrio for SSP 2-4.5 at 10% ONLY
#===============================================================================

scenario <- "245"
percentile <- 10

# No climate metric is used here – Climate-uninformed approach
# Directly use the species distribution data
# Skip the CPA step because no filtering based on climate

# Prepare spatial data and targets
out_sf <- prepare_out_sf(out_spp_C_bin, scenario)

targets <- C_target_df

target_check  <- clean_targets(out_sf, targets)
targets_fixed <- target_check$targets_fixed

# Build and solve spatial prioritisation problem
p <- prepare_problem(out_sf, targets_fixed)

set.seed(123)

solution_uninformed <- solve_problem(
  p,
  "Climate-uninformed",
  file.path(sol_dir, paste0("s-CPA-", percentile, "-", scenario, "-Uninformed.rds"))
)

#===============================================================================
# QUICK FULL-REGION SOLUTION PLOT — CLIMATE-UNINFORMED
#===============================================================================

# Load archived solution used for manuscript figures.
# This avoids re-solving the prioritisation problem when only regenerating HD plots.
archived_solution_path_uninformed <- file.path(
  "Data", "Solutions", "s-CPA-10-245-Uninformed.rds"
)

solution_uninformed <- readRDS(archived_solution_path_uninformed)

# Define tighter full-region bounding box in lon/lat
bbox_full <- st_as_sfc(
  st_bbox(c(xmin = 98, xmax = 166, ymin = -12, ymax = 19), crs = 4326)
)

# Transform bounding box to the CRS of the solution
bbox_full <- st_transform(bbox_full, st_crs(solution_uninformed))
bbox_full <- st_bbox(bbox_full)

# Add breathing room around the map
xpad <- (bbox_full["xmax"] - bbox_full["xmin"]) * 0.05
ypad <- (bbox_full["ymax"] - bbox_full["ymin"]) * 0.05

Uninformed_Full_Region_245_10 <- ggplot() +
  geom_sf(
    data = solution_uninformed,
    aes(fill = as.factor(solution_1)),
    colour = NA,
    size = 0.1
  ) +
  geom_sf(
    data = landmass,
    colour = "grey50",
    fill = "grey90",
    alpha = 0.9,
    size = 0.1
  ) +
  coord_sf(
    xlim = c(bbox_full["xmin"] - xpad, bbox_full["xmax"] + xpad),
    ylim = c(bbox_full["ymin"] - ypad, bbox_full["ymax"] + ypad),
    expand = FALSE
  ) +
  scale_fill_manual(
    values = c("1" = "#295E11", "0" = "#FED674"),
    breaks = c("0", "1"),
    labels = c("Not Selected", "Selected"),
    name   = "Planning Units"
  ) +
  theme_bw() +
  theme(
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text    = element_text(color = "black", size = 10),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 16)
  ) +
  labs(
    title = "Climate-uninformed",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(
  filename = file.path(fig_dir, "Uninformed_Full_Region_245_10.png"),
  plot     = Uninformed_Full_Region_245_10,
  width    = 9,
  height   = 4.5,
  dpi      = 300
)
#===============================================================================
# REGIONAL PLOTS
#===============================================================================

# Helper function to crop planning units using bounding box
crop_pus <- function(bbox_values) {
  PUs_region <- readRDS(file = file.path(rds_dir1, "PUs.rds"))
  
  bbox <- st_as_sfc(st_bbox(bbox_values, crs = 4326))
  bbox <- st_transform(bbox, st_crs(PUs_region))
  
  PUs_region %>%
    filter(st_intersects(geometry, bbox, sparse = FALSE))
}


# Helper function to plot regional solution
plot_region <- function(PUs_region, y_breaks, x_breaks = waiver()) {
  ggplot() +
    geom_sf(
      data = solution_uninformed,
      aes(fill = as.factor(solution_1)),
      colour = NA,
      size = 0.1,
      show.legend = FALSE
    ) +
    geom_sf(
      data = landmass,
      colour = "grey30",
      fill = "grey90",
      alpha = 0.9,
      size = 0.1,
      show.legend = FALSE
    ) +
    coord_sf(
      xlim = c(st_bbox(PUs_region)$xmin, st_bbox(PUs_region)$xmax),
      ylim = c(st_bbox(PUs_region)$ymin, st_bbox(PUs_region)$ymax)
    ) +
    scale_fill_manual(values = c("1" = "#295E11", "0" = "#FED674")) +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_continuous(breaks = x_breaks) +
    theme_bw() +
    theme(
      axis.ticks   = element_line(color = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.text    = element_text(color = "black", size = 14),
      axis.title   = element_blank(),
      plot.title   = element_blank()
    )
}


# Crop planning units
PUs_ph   <- crop_pus(c(xmin = 113, xmax = 126, ymin = 5, ymax = 14))
PUs_west <- crop_pus(c(xmin = 119, xmax = 132, ymin = -8, ymax = 1))
PUs_east <- crop_pus(c(xmin = 147, xmax = 160, ymin = -11.5, ymax = -2.5))


# Create regional plots
plot_ph <- plot_region(
  PUs_region = PUs_ph,
  y_breaks   = c(5, 10, 15, 20)
)

plot_west <- plot_region(
  PUs_region = PUs_west,
  y_breaks   = seq(-10, 5, by = 5),
  x_breaks   = seq(95, 130, by = 10)
)

plot_east <- plot_region(
  PUs_region = PUs_east,
  y_breaks   = seq(-10, 5, by = 5),
  x_breaks   = seq(130, 165, by = 10)
)


# Combine and save plots
combined_plot_uninformed <- (
  plot_ph +
    theme(plot.margin = margin(0, 0, 10, 0))
) / (
  plot_west +
    theme(plot.margin = margin(10, 0, 10, 0))
) / (
  plot_east +
    theme(plot.margin = margin(10, 0, 0, 0))
)

ggsave(
  filename = file.path(fig_dir, "Climate_Uninformed_Plan_245_10perc_PH.png"),
  plot     = combined_plot_uninformed,
  width    = 6.8,
  height   = 14,
  dpi      = 300
)