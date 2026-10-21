# 00_Dependencies: This R Script sets up all shared libraries, paths, and project setup

# Load required packages
library(here)
library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)

# Set project directory
project_dir <- here::here()

# Optional: set working directory to project root
setwd(project_dir)

# Define common folders
data_dir    <- file.path(project_dir, "Data")
scripts_dir <- file.path(project_dir, "Scripts")
out_dir     <- file.path(project_dir, "Outputs")
rds_dir     <- file.path(out_dir, "RDS")
figures_dir <- file.path(out_dir, "Figures")

# Specific output folders
rds_dir3 <- file.path(rds_dir, "03_Create_Features")

# Quick check
message("Project directory set to: ", project_dir)
message("Current working directory: ", getwd())