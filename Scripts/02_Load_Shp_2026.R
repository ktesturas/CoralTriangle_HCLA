# 02_Load_Shp: This R script loads all shapefiles necessary for spatial planning in the Coral Triangle Region onto R

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "02_Load_Shp"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Load libraries
library(sf)
library(tidyverse)
library(spatialplanr)

# Read Extent Spatial Datasets
Coraltxt <- file.path(
  "Data",
  "Features",
  "Coral Reefs",
  "Allen Coral Atlas Files",
  "Reef Extent",
  "CT_Reef_fin.shp"
) # Allen Coral Atlas Reef extent

# Read IUCN Spatial Datasets
IUCNCoraltxt <- file.path(
  "Data",
  "Features",
  "Coral Reefs",
  "IUCN",
  "IUCN+CT",
  "IUCNfin_cropped.shp"
) # IUCN Coral distribution

####### CORALS #######
# Read Allen Extent shapefile
Coral <- Coraltxt %>%
  st_read() %>%
  st_make_valid() %>%
  st_transform(crs = "ESRI:54009")
saveRDS(Coral, file = file.path(rds_dir, "Coral_Ext.rds"))

# Read Coral IUCN shapefile
IUCN_Coral <- IUCNCoraltxt %>%
  st_read() %>%
  st_make_valid() %>%
  st_transform(crs = "ESRI:54009")
saveRDS(IUCN_Coral, file = file.path(rds_dir, "IUCN_Coral.rds"))