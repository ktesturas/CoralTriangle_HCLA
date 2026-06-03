# 01_Setup_PUs: This R script prepares the Coral Triangle Planning Region

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "01_Setup_PUs"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Load Libraries
#install.packages("remotes")
#remotes::install_github("emlab-ucsb/spatialgridr")
library(spatialgridr) # https://emlab-ucsb.github.io/spatialgridr/
# gitcreds::gitcreds_set() # for reference: https://carpentries.github.io/sandpaper-docs/github-pat.html#creating-a-new-github-personal-access-token
# remotes::install_github("emlab-ucsb/oceandatr")
library(oceandatr) # https://emlab-ucsb.github.io/oceandatr/
library(sf)
library(tidyverse)
#remotes::install_github("https://github.com/MathMarEcol/utils4mme")
library(utils4mme)

# Get EEZs for each country
countries <- c("Philippines", "Malaysia", "Indonesia", "Solomon Islands", "Papua New Guinea", "East Timor")
eez <- purrr::map_df(countries, spatialgridr::get_boundary) %>%
  sf::st_union() %>% # Remove unneeded boundaries
  sf::st_sf() # Needed due to some errors in oceandatr at the moment. Will be fixed shortly

# Get PUs within the EEZ region
PUs <- spatialgridr::get_grid(boundary = eez, crs = "ESRI:54009", resolution = 25000, output = "sf_hex") %>%
  dplyr::mutate(cellID = row_number(),PU_area_km2 = as.numeric(units::set_units(sf::st_area(.), "km2")))
PUs
saveRDS(PUs, file = file.path(rds_dir, "PUs.rds")) # Save as RDS as file is heavy

# Plot to check
gg <- ggplot() + geom_sf(data = PUs, linewidth = 0.01) + geom_sf(data = eez, colour = "red", fill = NA)
ggsave(file.path(fig_dir, "PUs.png"), gg) # Its pretty big so better to save it to a file, rather than display on screen.

# Add PU area and Ecoregion
PUs <- readRDS(file = file.path(rds_dir, "PUs.rds")) %>%
  dplyr::mutate(cellID = row_number(),PU_area_km2 = as.numeric(units::set_units(sf::st_area(.), "km2"))) %>%
  match_MEOW()
saveRDS(PUs, file = file.path(rds_dir, "PUs.rds"))