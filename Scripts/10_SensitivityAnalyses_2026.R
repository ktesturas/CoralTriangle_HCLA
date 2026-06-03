# 10_RunCPAThresholdSensitivity_ssp245.R
# Run five CPA-based prioritization approaches across increasing percentile thresholds for all emission scenarios

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Load packages
library(spatialplanr)
library(prioritizr)
library(dplyr)
library(sf)
library(here)

script_name <- "10_SensitivityAnalyeses"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)
sol_dir <- file.path(out_dir, "Solutions", script_name)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)

# Set file paths for RDS
rds_dir3 <- file.path(out_dir, "RDS", "03_Create_Features")     # to retrieve out_Coral and PUs.fin_spp_C.rds

# Load PUs and calculate median coral area
PUs <- readRDS(file.path(rds_dir3, "out_Coral.rds")) %>%
  dplyr::select(-any_of(c("EcoRegion", "Province", "Realm", "LatZone", "Area_MEOW")))

median_coral_area <- median(PUs$Coral_area_km2, na.rm = TRUE)

# Keep only PUs with coral areas >= median coral area
PUs <- PUs %>%
  filter(Coral_area_km2 >= median_coral_area)

# Load coral species occurrence data
out_spp_C <- readRDS(file.path(rds_dir3, "PUs.fin_spp_C.rds"))

# Load functions
source("Utils/fClimatePriorityArea_CSapproach.R")
source("Utils/fClimateFeatureOverlap_CSapproach.R")
source("Utils/fFeature_CSapproach.R")
source("Utils/fUtils_Plotting.R")

# Prepare species data
out_spp_C_bin <- out_spp_C %>%
  as_tibble() %>%
  mutate(across(-any_of(c("cellID", "geometry")), ~ if_else(. >= 0.5, 1, 0))) %>%
  select(-which(colSums(select(., -geometry)) == 0)) %>%
  sf::st_as_sf()

out_spp_C_bin <- out_spp_C_bin %>%
  filter(cellID %in% PUs$cellID) %>%
  left_join(st_drop_geometry(PUs), by = "cellID") %>%
  select(-PU_area_km2, -Coral_area_km2) %>%
  st_as_sf()

# Define available scenarios
scenarios <- c("126", "245", "585")

# Set directories
metrics_dir <- file.path("Data", "Metrics")
combi_dir   <- file.path("Outputs", "RDS", "05_Calc_Combi_Metrics")

# Loop through each scenario
for (scenario in scenarios) {
  
  # Load climate exposure metrics
  for (var in c("tos", "mhw")) {
    
    df <- readRDS(
      file = file.path(metrics_dir, paste0(var, "_ssp", scenario, "_roc.rds"))
    )
    
    colnames(df)[2] <- names(df)[2]  # retain original name
    
    assign(
      paste0(var, "_ssp", scenario),
      df,
      envir = .GlobalEnv
    )
  }
  
  # Load combined climate-smart metrics: Adaptation hotspots
  assign(
    paste0("Combined_ssp", scenario, "_AH"),
    readRDS(
      file.path(combi_dir, paste0("combined_metric_ssp", scenario, ".rds"))
    ) %>%
      rename(metric = names(.)[2])
  )
  
  # Load combined climate-smart metrics: Climate refugia
  assign(
    paste0("Combined_ssp", scenario, "_R"),
    readRDS(
      file.path(combi_dir, paste0("combined_metric_ssp", scenario, "_Refugia.rds"))
    ) %>%
      rename(metric = names(.)[2])
  )
  
  # Load normalized AH–Refugia difference metric
  assign(
    paste0("Combined_ssp", scenario, "_norm"),
    readRDS(
      file.path(combi_dir, paste0("Combined_ssp", scenario, "_AH_R.rds"))
    ) %>%
      select(cellID, normalized_diff, geometry) %>%
      rename(metric = normalized_diff)
  )
}

# Targets and cost
C_target_df <- out_spp_C_bin %>%
  as_tibble() %>%
  select(-geometry, -cellID) %>%
  names() %>%
  tibble::enframe(name = NULL, value = "feature") %>%
  mutate(target = 0.3)

UniformCost <- PUs %>% select(-Coral_area_km2) %>% rename(cost = PU_area_km2)

# Prepare spatial data
prepare_out_sf <- function(CPA, scenario) {
  tos_var <- get(paste0("tos_ssp", scenario))
  mhw_var <- get(paste0("mhw_ssp", scenario))
  combined_ah_var <- get(paste0("Combined_ssp", scenario, "_AH"))
  combined_r_var <- get(paste0("Combined_ssp", scenario, "_R"))
  combined_norm_var <- get(paste0("Combined_ssp", scenario, "_norm"))
  
  out_df <- CPA %>%
    st_drop_geometry() %>%
    left_join(st_drop_geometry(UniformCost), by = "cellID") %>%
    left_join(st_drop_geometry(tos_var), by = "cellID") %>%
    left_join(st_drop_geometry(mhw_var), by = "cellID") %>%
    left_join(st_drop_geometry(combined_ah_var), by = "cellID") %>%
    rename(metric_AH = metric) %>%
    left_join(st_drop_geometry(combined_r_var), by = "cellID") %>%
    rename(metric_R = metric) %>%
    left_join(st_drop_geometry(combined_norm_var), by = "cellID") %>%
    rename(metric_norm = metric)
  
  out_sf <- st_sf(
    out_df,
    geometry = st_geometry(CPA),
    crs = st_crs(CPA)
  )
  
  return(out_sf)
}

# Prepare and solve prioritization
prepare_problem <- function(out_sf, targets) {
  problem(out_sf, features = targets$feature, cost_column = "cost") %>%
    add_min_set_objective() %>%
    add_relative_targets(targets$target) %>%
    add_binary_decisions() %>%
    add_rsymphony_solver(gap = 0.1, time_limit = 600, verbose = FALSE)
}

solve_problem <- function(problem, approach, file_name) {
  solution <- solve(problem, force = TRUE) %>%
    mutate(approach = approach)
  saveRDS(solution, file = file_name)
  return(solution)
}



###############################################
##### LOOP OVER SCENARIOS #####################
###############################################

# Define thresholds and scenarios
thresholds <- seq(0, 30, by = 1.5)
scenarios <- c("126", "245", "585")  # Add/remove scenarios as needed

# Initialize results summary dataframe and solution list
results_summary <- data.frame(
  threshold = numeric(),
  approach = character(),
  scenario = character(),
  area_km2 = numeric(),
  median_metric_R = numeric(),
  mean_metric_R = numeric(),
  median_metric_AH = numeric(),
  mean_metric_AH = numeric(),
  median_tos_transformed = numeric(),
  mean_tos_transformed = numeric(),
  median_mhw = numeric(),
  mean_mhw = numeric(),
  stringsAsFactors = FALSE
)

solutions_list <- list()  # Store each full solution object

# MAIN LOOP
for (scenario in scenarios) {
  message(paste("Starting scenario SSP", scenario))
  
  for (percentile in thresholds) {
    message(paste("  Processing threshold:", percentile))
    
    # Define approaches and their corresponding metric dataframes
    approaches <- list(
      "Climate-uninformed" = NULL,
      "Adaptation hotspots" = get(paste0("Combined_ssp", scenario, "_AH")),
      "Refugia" = get(paste0("Combined_ssp", scenario, "_R")),
      "Split" = list(
        AH = get(paste0("Combined_ssp", scenario, "_AH")),
        R = get(paste0("Combined_ssp", scenario, "_R"))
      ),
      "Combined" = get(paste0("Combined_ssp", scenario, "_norm"))
    )
    
    for (approach in names(approaches)) {
      message(paste("    →", approach))
      
      # Create CPA based on approach
      if (approach == "Climate-uninformed") {
        CPA <- out_spp_C_bin
        targets <- C_target_df
      } else if (approach == "Split") {
        CPA <- fClimatePriorityArea_CSapproach_Combined_bin(
          featuresDF = out_spp_C_bin,
          percentile = percentile,
          metricDF_R = approaches[[approach]]$R,
          metricDF_AH = approaches[[approach]]$AH,
          direction = 1
        )
        targets <- fAssignTargets_CPA_Combined(CPA, C_target_df, csTarget = 1)
      } else if (approach == "Combined") {
        CPA <- fClimateFeatureOverlap_CSapproach_Combined2_bin(
          featuresDF = out_spp_C_bin,
          percentile = percentile,
          metricDF = approaches[[approach]],
          direction = 1,
          upper_percentile = TRUE
        )
        targets <- fAssignTargets_CPA(CPA, C_target_df, refugiaTarget = 1)
      } else {
        CPA <- fClimatePriorityArea_CSapproach_bin(
          featuresDF = out_spp_C_bin,
          percentile = percentile,
          metricDF = approaches[[approach]],
          direction = 1
        )
        targets <- fAssignTargets_CPA(CPA, C_target_df, refugiaTarget = 1)
      }
      
      # Solve spatial prioritization
      out_sf <- prepare_out_sf(CPA, scenario)
      prob <- prepare_problem(out_sf, targets)
      s <- solve_problem(prob, approach,
                         file.path(sol_dir, paste0("soln_", approach, "_", percentile, "_ssp", scenario, ".rds")))
      
      # Store solution to list
      solution_name <- paste0("ssp", scenario, "_", approach, "_", percentile)
      solutions_list[[solution_name]] <- s
      
      # Compute summary statistics for selected units
      selected <- s %>% filter(solution_1 == 1)
      area_selected <- sum(selected$cost, na.rm = TRUE)
      
      median_R  <- if ("metric_R" %in% names(selected)) median(selected$metric_R, na.rm = TRUE) else NA
      mean_R    <- if ("metric_R" %in% names(selected)) mean(selected$metric_R, na.rm = TRUE) else NA
      median_AH <- if ("metric_AH" %in% names(selected)) median(selected$metric_AH, na.rm = TRUE) else NA
      mean_AH   <- if ("metric_AH" %in% names(selected)) mean(selected$metric_AH, na.rm = TRUE) else NA
      
      median_tos <- if ("tos_transformed" %in% names(selected)) median(selected$tos_transformed, na.rm = TRUE) else NA
      mean_tos   <- if ("tos_transformed" %in% names(selected)) mean(selected$tos_transformed, na.rm = TRUE) else NA
      median_mhw <- if ("mhw" %in% names(selected)) median(selected$mhw, na.rm = TRUE) else NA
      mean_mhw   <- if ("mhw" %in% names(selected)) mean(selected$mhw, na.rm = TRUE) else NA
      
      # Append to results summary
      results_summary <- rbind(results_summary, data.frame(
        threshold = percentile,
        approach = approach,
        scenario = paste0("ssp", scenario),
        area_km2 = area_selected,
        median_metric_R = median_R,
        mean_metric_R = mean_R,
        median_metric_AH = median_AH,
        mean_metric_AH = mean_AH,
        median_tos_transformed = median_tos,
        mean_tos_transformed = mean_tos,
        median_mhw = median_mhw,
        mean_mhw = mean_mhw,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save summary table per scenario
  write.csv(results_summary %>% filter(scenario == paste0("ssp", scenario)),
            file.path(rds_dir, paste0("summary_ssp", scenario, ".csv")),
            row.names = FALSE)
  saveRDS(results_summary %>% filter(scenario == paste0("ssp", scenario)),
          file.path(rds_dir, paste0("summary_ssp", scenario, ".rds")))
}

# Save full results summary and all solution objects
write.csv(results_summary, file.path(rds_dir, "summary_all.csv"), row.names = FALSE)
saveRDS(results_summary, file.path(rds_dir, "summary_all.rds"))
saveRDS(solutions_list, file.path(sol_dir, "all_solutions_list.rds"))
