# 11_Plot_Sensitivity_Analyses

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Load packages
library(dplyr)
library(sf)
library(here)
library(ggplot2)
library(patchwork)
library(stringr)
library(purrr)

script_name <- "11_Plot_SensitivityAnalyeses"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)
sol_dir <- file.path(out_dir, "Solutions", script_name)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)

####################################################################
# Recreate CPA threshold sensitivity summary from saved solution RDS
####################################################################

# Folder containing the saved soln_*.rds files
cpa_loop_dir <- here::here("Data", "CPA_Loop")

# Load planning units so percent_area can be recalculated
rds_dir3 <- file.path(out_dir, "RDS", "03_Create_Features")

PUs <- readRDS(file.path(rds_dir3, "out_Coral.rds")) %>%
  dplyr::select(-any_of(c("EcoRegion", "Province", "Realm", "LatZone", "Area_MEOW")))

median_coral_area <- median(PUs$Coral_area_km2, na.rm = TRUE)

PUs <- PUs %>%
  filter(Coral_area_km2 >= median_coral_area)

n_total_units <- nrow(PUs)

####################################################################
# Helper function to summarise one solution file
####################################################################

summarise_solution_file <- function(file_path) {
  
  file_name <- basename(file_path)
  
  # Extract approach, threshold, and scenario from file name
  file_info <- stringr::str_match(
    file_name,
    "^soln_(.*)_([0-9.]+)_ssp([0-9]+)\\.rds$"
  )
  
  approach_raw <- file_info[, 2]
  threshold_value <- as.numeric(file_info[, 3])
  scenario_value <- paste0("ssp", file_info[, 4])
  
  s <- readRDS(file_path)
  
  selected <- s %>%
    filter(solution_1 == 1)
  
  area_selected <- sum(selected$cost, na.rm = TRUE)
  
  median_R <- if ("metric_R" %in% names(selected)) {
    median(selected$metric_R, na.rm = TRUE)
  } else {
    NA
  }
  
  mean_R <- if ("metric_R" %in% names(selected)) {
    mean(selected$metric_R, na.rm = TRUE)
  } else {
    NA
  }
  
  median_AH <- if ("metric_AH" %in% names(selected)) {
    median(selected$metric_AH, na.rm = TRUE)
  } else {
    NA
  }
  
  mean_AH <- if ("metric_AH" %in% names(selected)) {
    mean(selected$metric_AH, na.rm = TRUE)
  } else {
    NA
  }
  
  median_tos <- if ("tos_transformed" %in% names(selected)) {
    median(selected$tos_transformed, na.rm = TRUE)
  } else {
    NA
  }
  
  mean_tos <- if ("tos_transformed" %in% names(selected)) {
    mean(selected$tos_transformed, na.rm = TRUE)
  } else {
    NA
  }
  
  median_mhw <- if ("mhw" %in% names(selected)) {
    median(selected$mhw, na.rm = TRUE)
  } else {
    NA
  }
  
  mean_mhw <- if ("mhw" %in% names(selected)) {
    mean(selected$mhw, na.rm = TRUE)
  } else {
    NA
  }
  
  data.frame(
    threshold = threshold_value,
    approach = approach_raw,
    scenario = scenario_value,
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
  )
}

####################################################################
# Rebuild summary table from all solution RDS files
####################################################################

solution_files <- list.files(
  cpa_loop_dir,
  pattern = "^soln_.*_ssp[0-9]+\\.rds$",
  full.names = TRUE
)

summary_df <- purrr::map_dfr(solution_files, summarise_solution_file)

# Rename approaches
summary_df <- summary_df %>%
  mutate(
    approach = recode(
      approach,
      "Adaptation hotspots" = "Climate adaptation hotspots",
      "Refugia" = "Climate refugia"
    )
  )

# Recalculate percent_area using equal-area assumption
summary_df <- summary_df %>%
  mutate(
    percent_area = (area_km2 / median(PUs$PU_area_km2, na.rm = TRUE)) / n_total_units * 100
  )

# Arrange for cleaner plotting and saving
summary_df <- summary_df %>%
  arrange(scenario, approach, threshold)

####################################################################
# Save rebuilt summary tables
####################################################################

saveRDS(
  summary_df,
  file.path(rds_dir, "summary_all_rebuilt.rds")
)

write.csv(
  summary_df,
  file.path(rds_dir, "summary_all_rebuilt.csv"),
  row.names = FALSE
)

for (scenario in unique(summary_df$scenario)) {
  
  saveRDS(
    summary_df %>% filter(scenario == !!scenario),
    file.path(rds_dir, paste0("summary_", scenario, "_rebuilt.rds"))
  )
  
  write.csv(
    summary_df %>% filter(scenario == !!scenario),
    file.path(rds_dir, paste0("summary_", scenario, "_rebuilt.csv")),
    row.names = FALSE
  )
}

####################################################################
# Plot CPA threshold sensitivity analyses
####################################################################

# Define color palette for approaches
approach_colors <- c(
  "Climate-uninformed" = "#8c510a",
  "Climate refugia" = "#d8b365",
  "Climate adaptation hotspots" = "#5ab4ac",
  "Split" = "#01665e",
  "Combined" = "#c7eae5"
)

# Define theme
custom_theme <- theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )
library(cowplot)

# Function to generate selected plots for a given scenario
make_plots_per_scenario <- function(df, scenario_name) {
  
  df <- df %>%
    filter(scenario == scenario_name)
  
  p1 <- ggplot(df, aes(x = threshold, y = percent_area, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Area selected (%)",
      title = "Area selected (%)"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p2 <- ggplot(df, aes(x = threshold, y = mean_metric_R, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Mean climate refugia score",
      title = "Mean climate refugia score"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p3 <- ggplot(df, aes(x = threshold, y = median_metric_R, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Median climate refugia score",
      title = "Median climate refugia score"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p4 <- ggplot(df, aes(x = threshold, y = mean_metric_AH, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Mean adaptation hotspot score",
      title = "Mean adaptation hotspot score"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p5 <- ggplot(df, aes(x = threshold, y = median_metric_AH, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Median adaptation hotspot score",
      title = "Median adaptation hotspot score"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p6 <- ggplot(df, aes(x = threshold, y = mean_tos_transformed, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = expression("Mean warming score (" * Delta * degree * "C " * yr^{-1} * ")"),
      title = bquote(bold("Mean warming score (" * Delta * degree * "C " * yr^{-1} * ")"))
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p7 <- ggplot(df, aes(x = threshold, y = median_tos_transformed, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = expression("Median warming score (" * Delta * degree * "C " * yr^{-1} * ")"),
      title = expression("Median warming score (" * Delta * degree * "C " * yr^{-1} * ")")
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p8 <- ggplot(df, aes(x = threshold, y = mean_mhw, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Mean MHW score (Total DHD)",
      title = "Mean MHW score (Total DHD)"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  p9 <- ggplot(df, aes(x = threshold, y = median_mhw, color = approach)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = approach_colors) +
    labs(
      x = "Percentile threshold",
      y = "Median MHW score",
      title = "Median MHW score"
    ) +
    custom_theme +
    theme(plot.margin = margin(12, 12, 28, 12))
  
  legend_plot <- p1 +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.2, "lines")
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        override.aes = list(linewidth = 1.5, size = 3)
      )
    )
  
  legend_only <- cowplot::get_legend(legend_plot)
  
  # Remove legends from individual plots
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p4 <- p4 + theme(legend.position = "none")
  p6 <- p6 + theme(legend.position = "none")
  p8 <- p8 + theme(legend.position = "none")
  
  # Arrange plots with vertical legend in the empty second-row space
  final_plot <- (p1 | p2 | p4) /
    (p6 | p8 | wrap_elements(legend_only)) +
    plot_layout(
      heights = c(1, 1),
      widths = c(1, 1, 0.8)
    )
  
  return(final_plot)
}

# Scenarios to plot
scenarios_to_plot <- c("ssp126", "ssp245", "ssp585")

# Generate, show, and save plots per scenario
for (scenario in scenarios_to_plot) {
  
  plot <- make_plots_per_scenario(summary_df, scenario)
  
  print(plot)
  
  ggsave(
    filename = file.path(fig_dir, paste0("figure_climate_smart_metrics_", scenario, ".png")),
    plot = plot,
    width = 16,
    height = 7,
    dpi = 300
  )
}
