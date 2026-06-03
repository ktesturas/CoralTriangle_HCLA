# Planning Layers

library(here)
source(here::here("Scripts", "00_Dependencies.R"))

script_name <- "Planning Layers"

fig_dir <- file.path(out_dir, "Figures", script_name)
rds_dir <- file.path(out_dir, "RDS", script_name)
sol_dir <- file.path(out_dir, "Solutions", script_name)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)


# Coral Diagnostics
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(here)
library(scales)
library("rnaturalearthdata")

# Filter out rows that fall below the median
# Load data sets
rds_dir3 <- file.path(out_dir, "RDS", "03_Create_Features")

PUs <- readRDS(file.path(rds_dir3, "out_Coral.rds"))

median_coral_area <- median(PUs$Coral_area_km2, na.rm = TRUE)
PUs<- PUs %>%
  filter(Coral_area_km2 >= median_coral_area)
out_spp_C <- readRDS(file.path(rds_dir3, "PUs.fin_spp_C.rds"))

# Load functions
source("Utils/fClimatePriorityArea_CSapproach.R") # CPA approach
source("Utils/fClimateFeatureOverlap_CSapproach.R") # CF overlap approach
source("Utils/fFeature_CSapproach.R") # Feature approach
source("Utils/fUtils_Plotting.R") # For plotting

# Convert feature probability values to binary data using median values
out_spp_C_bin <- out_spp_C %>%
  as_tibble() %>%
  mutate(across(-any_of(c("cellID", "geometry")), ~ case_when(. >= 0.5 ~ 1, TRUE ~ 0))) %>%
  select(-which(colSums(select(., -geometry)) == 0)) %>%
  sf::st_as_sf()

out_spp_C_bin <- out_spp_C_bin %>%
  filter(cellID %in% PUs$cellID) %>%
  left_join(st_drop_geometry(PUs), by = "cellID") %>%
  st_as_sf()


# A. Planning region and coral features 
# 1. Load planning region and data

rds_dir1 <- file.path(out_dir, "RDS", "01_Setup_PUs") # to retrieve PUs

PUs <- readRDS(file.path(rds_dir1, "PUs.rds"))
landmass <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(crs = "ESRI:54009")

# 2. Prepare the data
PUs <- PUs %>%
  mutate(has_coral = ifelse(cellID %in% out_spp_C_bin$cellID, TRUE, FALSE))


# 3. Plot planning region with coral features highlighted
plan_dom<-ggplot() +
  geom_sf(data = PUs, aes(fill = has_coral), size = 0.1) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "white"), name = "Coral Features") +
  geom_sf(data = landmass, fill = "grey30", color = "black") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 27),
    plot.title = element_text(size = 27, hjust = 0.5), # Center and enlarge the plot title
    legend.position = "none",  # Move the legend to the right
    legend.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 10)),  # Bold legend title and add margin
    legend.text = element_text(size = 22, hjust = 0.5),   # Center legend text
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(2.00, "cm"),
    legend.box.margin = margin(l = 10),  # Move the legend box to the right
    strip.text = element_text(size = 27),  # Increase facet label text size
    plot.margin = margin(t = 20, r = 40, b = 20, l = 20)  # Add margin around the plot to avoid cutting off axis text
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
plan_dom

ggsave(
  filename = file.path(fig_dir, "Planning_Region.png"),
  plot     = plan_dom,
  width    = 9,
  height   = 6,
  dpi      = 300
)





# B. Plot species richnes
Coral_spp_rich <- out_spp_C_bin %>%
  dplyr::select(
    -any_of(c(
      "EcoRegion",
      "Province",
      "Realm",
      "LatZone",
      "Area_MEOW",
      "PU_area_km2",
      "Coral_area_km2"
    ))
  ) %>% 
  sf::st_drop_geometry() %>%
  mutate(across(where(is.numeric) & !matches("cellID"), ~ as.numeric(. >= 0.5))) %>% # Convert to binary if needed
  mutate(richness = rowSums(select(., where(is.numeric), -cellID))) %>%
  select(cellID, richness) %>%
  sf::st_as_sf(., geometry = sf::st_geometry(out_spp_C_bin)) %>%
  filter(richness > 0)

# Filter out values below 300
Coral_spp_rich <- Coral_spp_rich %>%
  filter(richness >= 298)


breaks <- c(300, 400, 500, 600)

spp_rich <- ggplot(Coral_spp_rich) +
  geom_sf(aes(fill = richness), color = NA) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Species Richness",
    breaks = breaks,
    limits = c(300, 600),
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 20,
      barheight = 1
    )
  ) +
  geom_sf(data = landmass, fill = "grey30", color = "black") +
  coord_sf(xlim = st_bbox(PUs)$xlim, ylim = st_bbox(PUs)$ylim, expand = TRUE) +
  scale_y_continuous(
    breaks = seq(-90, 90, by = 10),
    labels = function(y) {
      ifelse(y == 0, "0°",
             ifelse(y > 0, paste0(y, "°N"), paste0(abs(y), "°S")))
    }
  ) +
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 27),
    plot.title = element_text(size = 27, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 22, hjust = 0.5, margin = margin(b = 5)),
    legend.text = element_text(size = 22, hjust = 0.5),
    legend.key.width = unit(3.50, "cm"),
    legend.key.height = unit(0.50, "cm"),
    legend.box.margin = margin(l = -10),
    strip.text = element_text(size = 27),
    plot.margin = margin(t = 20, r = 40, b = 20, l = 20)
  )

spp_rich

ggsave(
  filename = file.path(fig_dir, "Species_Richness.png"),
  plot     = spp_rich,
  width    = 11,
  height   = 7,
  dpi      = 300
)



# C) Plot frequency distribution of coral area in a planning unit

# Load coral planning units
PUs_all <- readRDS(file.path(rds_dir3, "out_Coral.rds"))

# Calculate rounded-up median coral area cut-off
median_coral_area <- ceiling(median(PUs_all$Coral_area_km2, na.rm = TRUE))

# Create retained planning units separately
PUs_retained <- PUs_all %>%
  filter(Coral_area_km2 >= median_coral_area)

# Check how many planning units are retained
n_retained <- nrow(PUs_retained)

# Plot full frequency distribution with cut-off line
right <- ggplot(PUs_all, aes(x = Coral_area_km2)) +
  geom_histogram(
    binwidth = 4,
    fill = "skyblue",
    color = "black"
  ) +
  geom_vline(
    xintercept = median_coral_area,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    x = "Reef area (km²)",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "top",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 28),
    plot.title = element_text(size = 26),
    legend.text = element_text(size = 20),
    axis.line.x.top = element_line(color = "black", linewidth = 1),
    axis.line.y.right = element_line(color = "black", linewidth = 1)
  )

right

ggsave(
  filename = file.path(fig_dir, "Histogram_right_skew.png"),
  plot = right,
  width = 10,
  height = 7,
  dpi = 300
)

