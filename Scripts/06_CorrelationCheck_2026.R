# 06_Correlation_Tests: This R script tests and plots the correlation between adaptation hotspot and climate refugia scores

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))

# Set output directories
script_name <- "06_Correlation_Tests"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)
csv_dir <- file.path(out_dir, "CSV", script_name)

# Set filepaths for RDS
rds_dir5 <- file.path(out_dir, "RDS", "05_Calc_Combi_Metrics")  # to retrieve metric layers

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# Load libraries
library(dplyr)
library(sf)
library(ggplot2)
library(readr)

# Load combined metric layers
combined_metric_ssp126 <- readRDS(file.path(rds_dir5, "combined_metric_ssp126.rds"))
combined_metric_ssp245 <- readRDS(file.path(rds_dir5, "combined_metric_ssp245.rds"))
combined_metric_ssp585 <- readRDS(file.path(rds_dir5, "combined_metric_ssp585.rds"))

combined_metric_ssp126_refugia <- readRDS(file.path(rds_dir5, "combined_metric_ssp126_Refugia.rds"))
combined_metric_ssp245_refugia <- readRDS(file.path(rds_dir5, "combined_metric_ssp245_Refugia.rds"))
combined_metric_ssp585_refugia <- readRDS(file.path(rds_dir5, "combined_metric_ssp585_Refugia.rds"))

# Function to prepare joined data
prep_cor_data <- function(ah_sf, r_sf, scenario_label) {
  
  ah_name <- names(ah_sf)[2]
  r_name  <- names(r_sf)[2]
  
  ah_tbl <- ah_sf %>%
    st_drop_geometry() %>%
    dplyr::select(cellID, AH = all_of(ah_name))
  
  r_tbl <- r_sf %>%
    st_drop_geometry() %>%
    dplyr::select(cellID, R = all_of(r_name))
  
  joined <- ah_tbl %>%
    left_join(r_tbl, by = "cellID") %>%
    filter(!is.na(AH), !is.na(R)) %>%
    mutate(scenario = scenario_label)
  
  return(joined)
}

# Build dataframe for all scenarios
cor_plot_df <- bind_rows(
  prep_cor_data(combined_metric_ssp126, combined_metric_ssp126_refugia, "SSP1-2.6"),
  prep_cor_data(combined_metric_ssp245, combined_metric_ssp245_refugia, "SSP2-4.5"),
  prep_cor_data(combined_metric_ssp585, combined_metric_ssp585_refugia, "SSP5-8.5")
)

# Function to run correlation tests
run_cor_test <- function(df) {
  pearson_test <- cor.test(df$AH, df$R, method = "pearson")
  spearman_test <- cor.test(df$AH, df$R, method = "spearman", exact = FALSE)
  
  tibble(
    n = nrow(df),
    pearson_r = unname(pearson_test$estimate),
    pearson_p = pearson_test$p.value,
    pearson_ci_low = pearson_test$conf.int[1],
    pearson_ci_high = pearson_test$conf.int[2],
    spearman_rho = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value
  )
}

# Run tests for each scenario
cor_results <- cor_plot_df %>%
  group_by(scenario) %>%
  group_modify(~ run_cor_test(.x)) %>%
  ungroup()
print(cor_results)

# Save test results
# write_csv(cor_results, file.path(csv_dir,"correlation_tests_AH_vs_Refugia.csv"))
# saveRDS(cor_results, file.path(rds_dir, "correlation_tests_AH_vs_Refugia.rds"))

# Prepare labels for scatter plot
cor_labels <- cor_plot_df %>%
  group_by(scenario) %>%
  summarise(
    pearson_r = cor(AH, R, method = "pearson", use = "complete.obs"),
    spearman_rho = cor(AH, R, method = "spearman", use = "complete.obs"),
    x = min(R, na.rm = TRUE) + 0.05 * diff(range(R, na.rm = TRUE)),
    y = max(AH, na.rm = TRUE) - 0.05 * diff(range(AH, na.rm = TRUE)),
    label = paste0(
      "Pearson r = ", round(pearson_r, 2),
      "\nSpearman rho = ", round(spearman_rho, 2)
    ),
    .groups = "drop"
  )

# Scatter plot
gg_cor <- ggplot(cor_plot_df, aes(x = R, y = AH)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_text(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  facet_wrap(~ scenario, ncol = 3) +
  labs(
    x = "Climate refugia score",
    y = "Climate adaptation hotspot score",
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    strip.text = element_text(size = 13)
  )

# Save plot
ggsave(
  filename = file.path(fig_dir, "correlation_AH_vs_Refugia.png"),
  plot = gg_cor,
  width = 15,
  height = 5,
  dpi = 300
)