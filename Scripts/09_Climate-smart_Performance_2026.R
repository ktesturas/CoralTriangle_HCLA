# 09_Climate-smart_Performance: This R script evaluates the performance of climate-smart spatial prioritisation approaches.

#===============================================================================
# Set output directories
#===============================================================================

# Load shared dependencies and project paths
library(here)
source(here::here("Scripts", "00_Dependencies.R"))
source("Utils/fUtils_Plotting.R")

script_name <- "09_Climate-smart_Performance"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Set filepath for archived spatial prioritisation solutions
archived_sol_dir <- file.path("Data", "Solutions")

#===============================================================================
# Load required packages
#===============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sf)
library(corrplot)
library(irr)
library(RColorBrewer)
library(ggridges)


#===============================================================================
# Set scenario and percentile
#===============================================================================

scenario <- "245"
percentile <- 10

#===============================================================================
# Load archived spatial prioritisation solutions
#===============================================================================
# Note:
# These archived RDS files replace the need to rerun spatial prioritisation.
# They are expected to contain the following columns:
# cellID, solution_1, tos_transformed, mhw, metric_R, metric_AH, metric_norm, cost

solution_uninformed <- readRDS(
  file.path(archived_sol_dir, "s-CPA-10-245-Uninformed.rds")
)

solution_AH <- readRDS(
  file.path(archived_sol_dir, "s-CPA-10-245.rds")
)

solution_R <- readRDS(
  file.path(archived_sol_dir, "s-CPA-10-245-Refugia.rds")
)

solution_Split <- readRDS(
  file.path(archived_sol_dir, "s-CPA-10-245-Combination1.rds")
)

solution_Combined <- readRDS(
  file.path(archived_sol_dir, "s-CPA-10-245-Combination2.rds")
)

#===============================================================================
# Prepare combined dataframe for performance plots
#===============================================================================

solution_list <- list(
  solution_uninformed,
  solution_AH,
  solution_R,
  solution_Split,
  solution_Combined
)

approach_names <- c(
  "Climate-uninformed",
  "Climate adaptation hotspots",
  "Climate refugia",
  "Split",
  "Combined"
)

# Assign approach names to each solution
for (i in seq_along(solution_list)) {
  solution_list[[i]]$approach <- approach_names[i]
}

# Combine all solutions into one dataframe
df <- solution_list %>%
  map_dfr(~ .x %>%
            dplyr::select(
              cellID,
              solution_1,
              tos_transformed,
              mhw,
              metric_R,
              metric_AH,
              metric_norm,
              approach,
              cost
            )) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    approach = factor(approach, levels = approach_names)
  )

# Save combined dataframe
saveRDS(
  df,
  file = file.path(rds_dir, paste0("Combined_df_ssp", scenario, ".rds"))
)


#===============================================================================
# Kernel density plots
#===============================================================================
# Note:
# This section uses your existing create_kernel_plot2() function.
# Make sure this function is loaded before running this script.

legend_text <- expression('Δ'^"o"*'C yr'^"-1"*'')
x_label <- expression('Rate of climate warming (Δ'^"o"*'C yr'^"-1"*')')

gg_01 <- create_kernel_plot2(
  df,
  "tos_transformed",
  x_label,
  "Climate Warming",
  file.path(fig_dir, paste0("CW_Kernel_", scenario, ".png")),
  "B",
  legend_text
)

ggsave(
  filename = file.path(fig_dir, paste("CW_Kernel_", scenario, "_", percentile, ".png")),
  plot = gg_01,
  width = 18.5,
  height = 11,
  dpi = 300
)


legend_text <- "total DHD"
x_label <- expression('Sum of cumulative MHW intensities (total DHD)')

gg_02 <- create_kernel_plot2(
  df,
  "mhw",
  x_label,
  "DHD",
  file.path(fig_dir, paste0("MHW_Kernel_", scenario, ".png")),
  "C",
  legend_text
)

ggsave(
  filename = file.path(fig_dir, paste("DHD_Kernel_", scenario, "_", percentile, ".png")),
  plot = gg_02,
  width = 18.5,
  height = 11,
  dpi = 300
)

#===============================================================================
# Cohen's kappa / correlation matrix
#===============================================================================

# Function to plot correlation matrix
fPlot_CorrPlot <- function(matrix) {
  
  rownames(matrix) <- matrix[, 1]
  matrix_f <- matrix[, -1]
  
  matrix_f <- as.matrix(matrix_f)
  storage.mode(matrix_f) <- "numeric"
  
  corrplot(
    matrix_f,
    method = "shade",
    is.corr = TRUE,
    cl.lim = c(-1, 1),
    tl.col = "black",
    tl.cex = 1.0,          # text size of row/column labels
    addCoef.col = "black",
    number.cex = 0.9,      # text size of correlation values
    col = COL2("BrBG", 200),
    tl.srt = 45
  )
}

# Function to select solution column and rename
select_solution <- function(obj, col_name) {
  
  obj %>%
    sf::st_drop_geometry() %>%
    as_tibble() %>%
    dplyr::select(solution_1) %>%
    dplyr::rename(!!sym(col_name) := solution_1)
}

# Prepare shorter approach names for the matrix
approach_names_matrix <- c(
  "Climate-uninformed",
  "Climate adaptation hotspots",
  "Climate refugia",
  "Split",
  "Combined"
)

# Build list of solution selections for matrix
object_list <- lapply(seq_along(approach_names_matrix), function(i) {
  select_solution(solution_list[[i]], approach_names_matrix[i])
})

# Calculate and plot the correlation matrix
matrix <- fGetCorrMatrix(object_list)

png(
  height = 1200,
  width = 1200,
  res = 200,
  file = file.path(fig_dir, paste0("CorrMatrix_Approaches_ssp", scenario, ".png")),
  type = "cairo"
)

Corrplot <- fPlot_CorrPlot(matrix)

dev.off()


#===============================================================================
# Diverging plot: climate refugia and adaptation hotspot scores
#===============================================================================

summary_df <- df %>%
  filter(solution_1 == 1) %>%
  group_by(approach) %>%
  summarise(
    metric_R_median = median(metric_R, na.rm = TRUE),
    metric_R_mean = mean(metric_R, na.rm = TRUE),
    metric_AH_median = median(metric_AH, na.rm = TRUE),
    metric_AH_mean = mean(metric_AH, na.rm = TRUE),
    .groups = "drop"
  )

plot_df_median <- summary_df %>%
  select(approach, metric_R_median, metric_AH_median) %>%
  pivot_longer(
    cols = c(metric_R_median, metric_AH_median),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value = ifelse(metric == "metric_R_median", -value, value),
    metric = ifelse(
      metric == "metric_R_median",
      "Climate refugia",
      "Climate adaptation hotspots"
    )
  )

plot_df_mean <- summary_df %>%
  select(approach, metric_R_mean, metric_AH_mean) %>%
  pivot_longer(
    cols = c(metric_R_mean, metric_AH_mean),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value = ifelse(metric == "metric_R_mean", -value, value),
    metric = ifelse(
      metric == "metric_R_mean",
      "Climate refugia",
      "Climate adaptation hotspots"
    )
  )

Diverging_Median_245 <- ggplot(
  plot_df_median,
  aes(x = value, y = approach, fill = metric)
) +
  geom_col(position = "identity", width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(
    values = c(
      "Climate refugia" = "#4F6DB8",
      "Climate adaptation hotspots" = "#CA3C25"
    )
  ) +
  scale_x_continuous(
    limits = c(-50, 50),
    breaks = seq(-50, 50, by = 10)
  ) +
  labs(
    x = "Metric score",
    y = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "gray60", linetype = "dotted"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 13)
  )

ggsave(
  filename = file.path(fig_dir, paste0("Diverging_Median_ssp", scenario, ".png")),
  plot = Diverging_Median_245,
  width = 7,
  height = 5,
  dpi = 300
)

Diverging_Mean_245 <- ggplot(
  plot_df_mean,
  aes(x = value, y = approach, fill = metric)
) +
  geom_col(position = "identity", width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(
    values = c(
      "Climate refugia" = "#4F6DB8",
      "Climate adaptation hotspots" = "#CA3C25"
    )
  ) +
  scale_x_continuous(
    limits = c(-50, 50),
    breaks = seq(-50, 50, by = 10)
  ) +
  labs(
    x = "Metric score",
    y = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "gray60", linetype = "dotted"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 13)
  )

ggsave(
  filename = file.path(fig_dir, paste0("Diverging_Mean_ssp", scenario, ".png")),
  plot = Diverging_Mean_245,
  width = 7,
  height = 5,
  dpi = 300
)

#===============================================================================
# Cost plot
#===============================================================================


df_selected <- df %>%
  filter(solution_1 == 1)

total_cost <- df %>%
  distinct(cellID) %>%
  nrow()

df_summary_245 <- df_selected %>%
  mutate(
    approach = case_when(
      approach == "Adaptation hotspots" ~ "Climate adaptation hotspots",
      approach == "Adaptation hotspot" ~ "Climate adaptation hotspots",
      approach == "Climate adaptation hotspot" ~ "Climate adaptation hotspots",
      TRUE ~ approach
    )
  ) %>%
  group_by(approach) %>%
  summarise(
    total_cost_approach = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percent_cost = (total_cost_approach / total_cost) * 100,
    approach = factor(
      approach,
      levels = c(
        "Climate-uninformed",
        "Climate adaptation hotspots",
        "Climate refugia",
        "Split",
        "Combined"
      )
    )
  )

fill_colors <- c(
  "Climate-uninformed" = "#1A936F",
  "Climate adaptation hotspots" = "#1A936F",
  "Climate refugia" = "#1A936F",
  "Split" = "#1A936F",
  "Combined" = "#1A936F"
)

Cost_group_245 <- ggplot(
  df_summary_245,
  aes(x = approach, y = percent_cost, fill = approach)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(),
    alpha = 1,
    colour = "black",
    size = 0.5
  ) +
  labs(
    x = "",
    y = "% Area"
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 40)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray60", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 20),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none",
    legend.title = element_blank(),
    panel.border = element_rect(size = 0.7),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 18, color = "black"),
    plot.margin = margin(t = 15, r = 15, b = 15, l = 45)
  )

Cost_group_245

ggsave(
  filename = file.path(fig_dir, paste0("Cost_Percent_Area_ssp", scenario, ".png")),
  plot = Cost_group_245,
  width = 7.5,
  height = 8,
  dpi = 300
)

#===============================================================================
# Save performance objects
#===============================================================================

saveRDS(
  summary_df,
  file = file.path(rds_dir, paste0("Metric_summary_ssp", scenario, ".rds"))
)

saveRDS(
  df_summary_245,
  file = file.path(rds_dir, paste0("Cost_summary_ssp", scenario, ".rds"))
)

saveRDS(
  matrix,
  file = file.path(rds_dir, paste0("CorrMatrix_Approaches_ssp", scenario, ".rds"))
)
