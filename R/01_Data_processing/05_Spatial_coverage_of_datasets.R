#----------------------------------------------------------#
# Fossil pollen data can predict robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#----------------------------------------------------------#

#----------------------------------------------------------#
#
#               Spatial coverage of datasets ----
#          
#----------------------------------------------------------#

#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# Load datasets ----
#-----------------------------------------------#
# Surface_samples 
surface_pollen_filtered_clim_zone <-
  readr::read_rds("Inputs/Data/surface_pollen_filtered_030124.rds") #2735 samples/datasets

# Fossil pollen records
top_fossil_pollen_500yr_filtered_clim_zone <-
  readr::read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds") %>% #363 samples
  dplyr::group_by(dataset_id, lat, long) %>% 
  tidyr::nest() %>%   # 71 datasets
  dplyr::ungroup()

# Make a base map with climate zones ----
# Read the raster points from the geo-tiff file published in Beck et al. 2018
beck_translation_table <- 
  read_csv("Inputs/Data/Biomes_spatial/koppen_link.csv") %>% 
  dplyr::rename(
    ecozone_koppen_30 = genzone,
    ecozone_koppen_15 = genzone_cluster,
    ecozone_koppen_5 = broadbiome
  ) %>% 
  # Aggregate the climate zones as suggested by John
  dplyr::mutate(
    climate_zone_revised = dplyr::case_when(
      ecozone_koppen_15 == "Arid_Desert" ~ "Arid",
      ecozone_koppen_15 == "Arid_Steppe" ~ "Arid",
      ecozone_koppen_15 == "Cold_Dry_Summer" ~ "Cold-seasonally dry",
      ecozone_koppen_15 == "Cold_Dry_Winter" ~ "Cold-seasonally dry",
      ecozone_koppen_15 == "Temperate_Dry_Summer" ~ "Temperate",
      ecozone_koppen_15 == "Temperate_Dry_Winter" ~ "Temperate",
      ecozone_koppen_15 == "Temperate_Without_dry_season" ~ "Temperate",
      ecozone_koppen_15 == "Polar_Tundra" ~ "Polar",
      ecozone_koppen_15 == "Polar_Frost" ~ "Polar frost",
      ecozone_koppen_15 == "Cold_Without_dry_season" ~ "Cold without dry season",
      ecozone_koppen_15 == "Tropical_Rainforest" ~ "Tropical rainforest",
      ecozone_koppen_15 == "Tropical_Monsoon" ~ "Tropical monsoon",
      ecozone_koppen_15 == "Tropical_Savannah" ~ "Tropical savannah",
      TRUE ~ ecozone_koppen_15
    )
  ) 

beck_raster_file <-
  raster::raster(
    here::here(
      "Inputs/Data/Biomes_spatial/Beck_KG_V1_present_0p083.tif"
    )
  )

# Extract the required raster points
raster_file <- aggregate(beck_raster_file, fact = 10) 
# Reduce the raster dimension (by factor of 10), otherwise would be too big file

raster_df <-
  # Convert raster points into a dataframe
  as.data.frame(beck_raster_file, xy = TRUE) %>%
  tibble::as_tibble() %>%
  # Extract the rater points of only required area
  dplyr::filter(x > 20 & y > 0) %>% # for required lat, long
  dplyr::rename(raster_values = Beck_KG_V1_present_0p083) %>%
  dplyr::mutate(raster_values = round(raster_values, digits = 0)) %>%
  # Assign the names of climate zone to the raster values
  dplyr::left_join(beck_translation_table, by = c("raster_values")) %>%
  dplyr::filter(!raster_values == 0) 

base_map <-
  surface_pollen_filtered_clim_zone %>%
  ggplot2::ggplot(aes(x = long, y = lat)) +
  ggplot2::coord_fixed(
    ylim = c(20.00, 80.00),
    xlim = c(65.00, 135.00)
    ) +
  ggplot2::geom_tile(
    data = raster_df,
    aes(x = x, y = y, fill = climate_zone_revised),
    inherit.aes = FALSE,
    alpha = 0.75
    ) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::labs(x = expression(paste('Longitude ', (degree ~ E))),
                y = expression(paste('Latitude ', (degree ~ N))),
                fill = "Climate-zones"
                ) +
  ggplot2::theme_classic() +
  ggplot2::borders(colour = color_common,
          linewidth = 0.2) +
  ggplot2::theme(
    axis.title = element_text(color = color_common, size = 24),
    axis.text = element_text(colour = color_common, size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 18), 
    legend.key.size = unit(0.75, 'cm')
    ) +
  ggplot2::guides(fill = guide_legend(title.position = "top"))

surface_samples_distribution <- 
  base_map +
  ggplot2::geom_point(color = "blue", size = 1.5) 
  
filtered_data_surface_samples <-
  surface_pollen_filtered_clim_zone %>%
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) 

insetrect <- 
  data.frame(xmin = min(filtered_data_surface_samples$long) - 0.75, 
             xmax = max(filtered_data_surface_samples$long) + 0.75,
             ymin = min(filtered_data_surface_samples$lat) - 0.75,  
             ymax = max(filtered_data_surface_samples$lat) + 0.75) 

final_fig_surface_sample_distribution <-
  surface_samples_distribution +
  ggplot2::geom_rect(
    data = insetrect,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax),
    alpha = 0,
    colour = color_common,
    linewidth = 1,
    inherit.aes = FALSE) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 
                           unit = "cm")
        )

#ggsave(final_fig_surface_sample_distribution,
#       file = "Outputs/Figure/surface_samples_distribution_260124.tiff", 
#       height = 15, 
#       width = 15, 
#       units = "cm", 
#       dpi = 400,
#       compression = "lzw")


# Fossil datasets
fossil_samples_distribution <- 
  base_map +
  ggplot2::geom_point(color = "blue", 
             size = 1.5, 
             data = top_fossil_pollen_500yr_filtered_clim_zone) 

filtered_data_fossil_samples <-
  surface_pollen_filtered_clim_zone %>%
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) 

insetrect <- 
  data.frame(xmin = min(filtered_data_fossil_samples$long) - 0.75, 
             xmax = max(filtered_data_fossil_samples$long) + 0.75,
             ymin = min(filtered_data_fossil_samples$lat) - 0.75,  
             ymax = max(filtered_data_fossil_samples$lat) + 0.75) 

final_fig_fossil_sample_distribution <-
  fossil_samples_distribution +
  ggplot2::geom_rect(
    data = insetrect,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax),
    alpha = 0,
    colour = color_common,
    linewidth = 1,
    inherit.aes = FALSE) 

#ggsave(final_fig_fossil_sample_distribution,
#       file = "Outputs/Figure/fossil_samples_distribution_260124.tiff", 
#       height = 15, 
#       width = 15, 
#       units = "cm", 
#       dpi = 400,
#       compression = "lzw")

final_map <- 
  ggpubr::ggarrange(
    final_fig_surface_sample_distribution,
    final_fig_fossil_sample_distribution,
    ncol = 2,
    nrow = 1,
    labels = c("(a)", "(b)"),
    label.x = 0.20,
    label.y = 0.99,
    font.label = list(size = 24, color = color_common),
    common.legend = TRUE,
    legend = "bottom"
    )

#ggsave(final_map,
#       file = "Outputs/Figure/samples_distribution_210524.tiff", 
#       height = 15,
#       width = 30,
#       units = "cm", 
#       dpi = 400,
#       compression = "lzw"
#       )
