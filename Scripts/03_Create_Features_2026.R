# 03_Create_Features: This R script interpolates the features onto the planning region

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "03_Create_Features"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)

# Set filepaths for RDS
rds_dir1 <- file.path(out_dir, "RDS", "01_Setup_PUs")  # to retrieve PUs
rds_dir2 <- file.path(out_dir, "RDS", "02_Load_Shp")  # to retrieve shapefiles/biodiversity data

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Load libraries
library(spatialplanr) # devtools::install_github("MathMarEcol/spatialplanr")
library(oceandatr)
library(tidyverse)
library(sf)
library(spatialgridr)

# Load Planning Region
PUs <- readRDS(file = file.path(rds_dir1, "PUs.rds"))

# Load Extent Datasets
Coral <- readRDS(file = file.path(rds_dir2, "Coral_Ext.rds")) %>%
  sf::st_sf() %>%
  dplyr::select(geometry)

# Load IUCN Datasets
IUCN_Coral <- readRDS(file = file.path(rds_dir2, "IUCN_Coral.rds")) %>%
  sf::st_sf() %>%
  dplyr::mutate(sci_name = stringr::str_replace_all(sci_name, " ", "_"))

# Interpolate Extent Data onto PUs
out <- Coral %>%
  sf::st_intersection(PUs) %>%
  dplyr::mutate(Coral_area_km2 = as.numeric(units::set_units(sf::st_area(.), "km2"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(cellID) %>%
  dplyr::summarise(Coral_area_km2 = sum(Coral_area_km2),
                   .groups = "drop_last") %>%
  dplyr::left_join(PUs, ., by = "cellID")
out_Coral <- out[!is.na(out$Coral_area_km2), ]
saveRDS(out_Coral, file = file.path(rds_dir, "out_Coral.rds"))

# Interpolate IUCN Datasets onto PUs
out <- IUCN_Coral %>%
  sf::st_make_valid() %>%
  sf::st_crop(sf::st_bbox(out)) %>%
  sf::st_cast("POLYGON") %>%
  sf::st_intersection(out %>% filter(Coral_area_km2 > 0)) %>%
  dplyr::mutate(IUCNCoral_area_km2 = as.numeric(units::set_units(sf::st_area(.), "km2"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(cellID) %>%
  dplyr::summarise(IUCNCoral_area_km2 = sum(IUCNCoral_area_km2),
                   .groups = "drop_last") %>%
  dplyr::left_join(out, ., by = "cellID")

out <- out %>%
  replace_na(list(Coral_area_km2 = 0,
                  IUCNCoral_area_km2 = 0))

saveRDS(out, file = file.path(rds_dir, "PUs.fin.rds"))
out <- readRDS(file = file.path(rds_dir, "PUs.fin.rds"))

# Get coral species distribution from IUCN as separate binary features
# Use only planning units with coral area
PUs_coral <- out %>%
  filter(Coral_area_km2 > 0) %>%
  dplyr::select(cellID, geometry)

# Clean IUCN coral polygons
IUCN_Coral_clean <- IUCN_Coral %>%
  sf::st_make_valid() %>%
  dplyr::select(sci_name, geometry)

# Identify species polygons that intersect each coral planning unit
# This is faster than st_intersection because it does not cut geometries
IUCN_Coral_PU_overlap <- PUs_coral %>%
  sf::st_join(
    IUCN_Coral_clean,
    join = sf::st_intersects,
    left = FALSE
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(cellID, sci_name) %>%
  dplyr::mutate(value = 1)

# Convert long overlap table to wide species-by-planning-unit matrix
out_spp_C_table <- IUCN_Coral_PU_overlap %>%
  tidyr::pivot_wider(
    names_from = sci_name,
    values_from = value,
    values_fill = 0
  )

# Join back to coral planning units
out_spp_C <- PUs_coral %>%
  dplyr::left_join(out_spp_C_table, by = "cellID") %>%
  dplyr::mutate(
    dplyr::across(
      -c(cellID, geometry),
      ~ tidyr::replace_na(., 0)
    )
  )

# Add C_ prefix to species columns
species_columns <- names(out_spp_C)[!names(out_spp_C) %in% c("geometry", "cellID")]
new_column_names <- paste0("C_", species_columns)
names(out_spp_C)[!names(out_spp_C) %in% c("geometry", "cellID")] <- new_column_names

# Save
saveRDS(out_spp_C, file = file.path(rds_dir, "PUs.fin_spp_C.rds"))

# Quick check
spp_cols <- names(out_spp_C)[grepl("^C_", names(out_spp_C))]

out_spp_C %>%
  sf::st_drop_geometry() %>%
  dplyr::select(all_of(spp_cols)) %>%
  dplyr::mutate(richness = rowSums(.)) %>%
  dplyr::summarise(
    min_richness = min(richness, na.rm = TRUE),
    mean_richness = mean(richness, na.rm = TRUE),
    max_richness = max(richness, na.rm = TRUE),
    unique_richness_values = length(unique(richness))
  )