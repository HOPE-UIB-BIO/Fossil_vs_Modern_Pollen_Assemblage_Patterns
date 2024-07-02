#----------------------------------------------------------#
# Fossil pollen data can predict robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#----------------------------------------------------------#

#----------------------------------------------------------#
#        Prepare data of phylodiversity estimation ----
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# A. Surface_samples ----
#-----------------------------------------------#

dat_harmonised <- 
  read_rds("Inputs/Data/data_harmonised_121223.rds")

dat_phylodiv <- 
  dat_harmonised %>% 
  dplyr::select(
    lat, 
    long,
    sample_id,
    harmonised_percentages
    ) %>% 
  dplyr::filter(!sample_id == "201") %>% #erroneous data
  tidyr::unnest(harmonised_percentages)  #2789 datasets

which(duplicated(dat_phylodiv$sample_id)) # 178
  
data_for_phylodiv_estimation <- 
dat_phylodiv %>% 
  dplyr::slice(-178,) %>% 
  dplyr::arrange(lat) %>% 
  tibble::column_to_rownames("sample_id") %>% 
  dplyr::mutate_if(is.numeric, 
                   ~replace(., is.na(.), 
                            0)
                   ) %>% 
  dplyr::select_if(colSums(.) != 0) 
  
  
data_for_phylodiv_estimation <- 
  data_for_phylodiv_estimation[,colSums(data_for_phylodiv_estimation > 0) > 2] # Filter out taxa with < 3 occurrences #107 taxa

surface_pollen_filtered <- 
  data_for_phylodiv_estimation[rowSums(data_for_phylodiv_estimation > 0) > 2,] %>% # filter out samples with < 3 taxa
  tibble::rownames_to_column("sample_id")  #2783 samples


# Add climate zones (Beck et al. 2018) ----

beck_translation_table <- 
  readr::read_csv("Inputs/Data/Biomes_spatial/koppen_link.csv") %>% 
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

surface_pollen_filtered_geo_tif <- 
  RUtilpol::geo_assign_tif(
    data_source = surface_pollen_filtered, 
    tif_file_name = "Inputs/Data/Biomes_spatial/Beck_KG_V1_present_0p083.tif",  
    fill_na = FALSE,
    na_as_value = NULL,
    distance_step = 500,
    n_max_step = 10)

surface_pollen_filtered_clim_zone <- 
  surface_pollen_filtered_geo_tif %>% 
  dplyr::left_join(
    beck_translation_table,
    by = "raster_values"
    ) %>% 
  dplyr::select(sample_id,
                lat,
                long,
                ecozone_koppen_5, 
                ecozone_koppen_15, 
                ecozone_koppen_30,
                climate_zone_revised,
                everything(),
                -raster_values) %>% 
  dplyr::filter(!is.na(climate_zone_revised)) # 2735 samples

#-----------------------------------------------#
# B. Top samples from fossil pollen data ----
#-----------------------------------------------#

# Note: This data was used directly from  "Latitudinal gradients in the 
#  phylogenetic assembly of angiosperms in Asia during the Holocene"
#  Therefore, initial data filtering and harmonisation was already done there.

top_fossil_pollen <-
  readr::read_rds("Inputs/Data/data_for_main_analysis_191223.rds") %>%
  dplyr::select(dataset_id,
                long,
                lat,
                levels_filtered,
                harmonised_fam_angiosperms_percentages) %>%
  dplyr::mutate(
    levels_filtered_500yr = purrr::map(levels_filtered,
                                       ~.x %>% 
                                         dplyr::filter(age >= 0 & age <= 500)
                                       ),
    nsamples_500yr = purrr::map_dbl(levels_filtered_500yr,
                                    ~ nrow(.x)
                                    )
    )

top_fossil_pollen_500yr <- 
  top_fossil_pollen %>% 
  dplyr::filter(nsamples_500yr > 0) %>%  
  dplyr::mutate(harmonised_fam_angiosperms_percentages_500yr = 
                  purrr::map2(
                    .x = harmonised_fam_angiosperms_percentages,
                    .y = levels_filtered_500yr,
                    ~ .x %>% 
                      dplyr::inner_join(.y %>% 
                                          dplyr::select(sample_id),
                                        by = "sample_id"
                      ) %>% 
                      dplyr::select(sample_id, everything())
                  )
                ) %>% 
  dplyr::select(
    dataset_id,
    lat,
    long,
    harmonised_fam_angiosperms_percentages_500yr
    ) %>% 
  tidyr::unnest(harmonised_fam_angiosperms_percentages_500yr) %>% 
  dplyr::mutate(sample_id = paste(sample_id, "_", row.names(.), sep = "")) # To avoid duplicated row names

lat_sample_id_top_fossil_pollen_500yr <- 
  top_fossil_pollen_500yr %>% 
  dplyr::select(dataset_id, lat, long, sample_id)

top_fossil_pollen_500yr_filtered_1 <- 
  top_fossil_pollen_500yr %>% 
  dplyr::select(-c(dataset_id, lat, long)) %>% 
  column_to_rownames("sample_id") %>% 
  dplyr::mutate_if(is.numeric, 
                   ~replace(., is.na(.), 
                            0)
                   ) %>% 
  dplyr::select_if(colSums(.) != 0)

top_fossil_pollen_500yr_filtered <- 
  top_fossil_pollen_500yr_filtered_1[rowSums(top_fossil_pollen_500yr_filtered_1 > 0) > 2,] %>% # filter out samples with < 3 taxa
  rownames_to_column("sample_id") %>% 
  dplyr::inner_join(lat_sample_id_top_fossil_pollen_500yr,
                    by = "sample_id") %>% 
  dplyr::select(dataset_id, sample_id, lat, long, everything()) # 79 taxa, 363 samples

# Add climate zone (Beck et al. 2018)
top_fossil_pollen_500yr_filtered_geo_tif <- 
  RUtilpol::geo_assign_tif(
    data_source = top_fossil_pollen_500yr_filtered, 
    tif_file_name = "Inputs/Data/Biomes_spatial/Beck_KG_V1_present_0p083.tif",  
    fill_na = FALSE,
    na_as_value = NULL,
    distance_step = 500,
    n_max_step = 10)

top_fossil_pollen_500yr_filtered_clim_zone <- 
  top_fossil_pollen_500yr_filtered_geo_tif %>% 
  dplyr::left_join(
    beck_translation_table,
    by = "raster_values"
    ) %>% 
  dplyr::select(
    dataset_id,
    sample_id,
    lat,
    long,
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised,
    everything(),
    -raster_values
    )

#-----------------------------------------------#
# Save filtered data ----
#-----------------------------------------------#

readr::write_rds(surface_pollen_filtered_clim_zone,
                 file = "Inputs/Data/surface_pollen_filtered_030124.rds",
                 compress = "gz")

readr::write_rds(top_fossil_pollen_500yr_filtered_clim_zone,
          file = "Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds",
          compress = "gz")
