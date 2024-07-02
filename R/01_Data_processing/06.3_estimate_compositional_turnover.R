#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#
#               Surface sample pollen data 
#          
#             Compositional turnover (DCCA axis 1) ----
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

# Load processed datasets
surface_pollen_filtered_clim_zone <- 
  read_rds("Inputs/Data/surface_pollen_filtered_030124.rds")
surface_samples_sp <- 
  surface_pollen_filtered_clim_zone %>% 
  dplyr::select(
    -c(
      lat, 
      long,
      ecozone_koppen_5,
      ecozone_koppen_15,
      ecozone_koppen_30,
      climate_zone_revised
      )
    )
surface_samples_lat <- 
  surface_pollen_filtered_clim_zone %>% 
  dplyr::select(
    sample_id, 
    lat, 
    long,
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised
    )


top_fossil_pollen_500yr_filtered_clim_zone <- 
  read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds")  
top_fossil_pollen_500yr_sp <- 
  top_fossil_pollen_500yr_filtered_clim_zone %>% 
  dplyr::select(
    -c(
      dataset_id,
      lat, 
      long,
      ecozone_koppen_5,
      ecozone_koppen_15,
      ecozone_koppen_30,
      climate_zone_revised
      )
    )
top_fossil_pollen_500yr_lat <- 
  top_fossil_pollen_500yr_filtered_clim_zone %>% 
  dplyr::select(
    dataset_id, 
    sample_id, 
    lat, 
    long,
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised
    )

top_fossil_pollen_1000yr_filtered_clim_zone <- 
  read_rds("Inputs/Data/top_fossil_pollen_1000yr_filtered_080124.rds")
top_fossil_pollen_1000yr_sp <- 
  top_fossil_pollen_1000yr_filtered_clim_zone %>% 
  dplyr::select(
    -c(
      dataset_id,
      lat, 
      long,
      ecozone_koppen_5,
      ecozone_koppen_15,
      ecozone_koppen_30,
      climate_zone_revised
      )
    )
top_fossil_pollen_1000yr_lat <- 
  top_fossil_pollen_1000yr_filtered_clim_zone %>% 
  dplyr::select(
    dataset_id, 
    sample_id, 
    lat, 
    long,
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised
    )

#modern_veg_data <- 
#  read_rds("Inputs/Data/data_gbif_for_phylodiversity_160124.rds") %>% 
#             dplyr::mutate(
#               sample_id = paste(
 #                "sample_", row.names(.), 
 #                sep = ""
 #              )
 #            ) %>% 
 #            dplyr::select(sample_id, lat, everything())
#modern_veg_data_sp <- 
 # modern_veg_data %>% 
 # dplyr::select(-lat)
#modern_veg_data_lat <- 
#  modern_veg_data %>% 
#  dplyr::select(sample_id, lat)

# DCCA
dcca_surface_samples <- 
  REcopol::fit_ordination(
    data_source_community = surface_samples_sp,
    data_source_predictors = surface_samples_lat,
    var_name_pred = "lat",
    sel_method = "constrained",
    sel_complexity = "poly_2",
    transform_to_percentage = FALSE,
    tranformation = "none"
    )
turnover_surface_samples <- 
  dcca_surface_samples$case_r %>% 
  #dplyr::select(sample_id, axis_1, axis_2) %>% 
  dplyr::inner_join(surface_samples_lat, by = "sample_id") %>% 
  dplyr::mutate(data_type = paste("surface_samples")) %>% 
  dplyr::group_by(lat) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(dataset_id = paste("dataset_", row.names(.), sep = "")) %>% 
  tidyr::unnest(data) %>% 
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    everything()
    ) 
  
#dplyr::mutate(DCCA_axis_1 = ifelse(sitename == "Reiersdalvann", 
#                                   (max(DCCA_axis_1) - DCCA_axis_1) + 
#                                      min(DCCA_axis_1), DCCA_axis_1)) # To invert the curve

dcca_top_fossil_pollen_500yr <- 
  REcopol::fit_ordination(
    data_source_community = top_fossil_pollen_500yr_sp,
    data_source_predictors = top_fossil_pollen_500yr_lat,
    var_name_pred = "lat",
    sel_method = "constrained",
    sel_complexity = "poly_2",
    transform_to_percentage = FALSE,
    tranformation = "none"
    )

turnover_top_fossil_pollen_500yr <- 
  dcca_top_fossil_pollen_500yr$case_r %>% 
  #dplyr::select(sample_id, axis_1, axis_2) %>% 
  dplyr::inner_join(top_fossil_pollen_500yr_lat, by = "sample_id") %>% 
  dplyr::mutate(axis_1 = max(axis_1) - axis_1 + min(axis_1)) %>% 
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    everything()
    )  %>% 
  dplyr::mutate(data_type = paste("fossil_pollen_top_500yr"))


dcca_top_fossil_pollen_1000yr <- 
  REcopol::fit_ordination(
    data_source_community = top_fossil_pollen_1000yr_sp,
    data_source_predictors = top_fossil_pollen_1000yr_lat,
    var_name_pred = "lat",
    sel_method = "constrained",
    sel_complexity = "poly_2",
    transform_to_percentage = FALSE,
    tranformation = "none"
    )
turnover_top_fossil_pollen_1000yr <- 
  dcca_top_fossil_pollen_1000yr$case_r %>% 
  #dplyr::select(sample_id, axis_1, axis_2) %>% 
  dplyr::inner_join(top_fossil_pollen_1000yr_lat, by = "sample_id") %>% 
  dplyr::mutate(axis_1 = max(axis_1) - axis_1 + min(axis_1)) %>% 
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    everything()
    ) %>%   
  dplyr::mutate(data_type = paste("fossil_pollen_top_1000yr"))  


#dcca_modern_samples <- 
#  REcopol::fit_ordination(
#    data_source_community = modern_veg_data_sp,
#    data_source_predictors = modern_veg_data_lat,
 #   var_name_pred = "lat",
 #   sel_method = "constrained",
 #   sel_complexity = "poly_2",
 #   transform_to_percentage = FALSE,
 #   tranformation = "none"
 # )
#turnover_modern_samples <- 
 # dcca_modern_samples$case_r %>% 
 # dplyr::select(sample_id, axis_1, axis_2) %>% 
 # dplyr::inner_join(modern_veg_data_lat, by = "sample_id") %>% 
 # dplyr::select(sample_id, lat, everything()) %>% 
#  dplyr::mutate(data_type = paste("modern_vegetation"))


turnover_combined <- 
  turnover_surface_samples %>% 
  bind_rows(turnover_top_fossil_pollen_500yr) %>% 
  bind_rows(turnover_top_fossil_pollen_1000yr) %>% 
  #bind_rows(turnover_modern_samples) %>% 
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    data_type,
    everything()
  ) 

write_rds(turnover_combined,
          file = "Inputs/Data/turnover_combined_050224.rds",
          compress = "gz")

#-------------------------------#
# DCA
#-------------------------------#
dca_surface_samples <- 
  REcopol::fit_ordination(
    data_source_community = surface_samples_sp,
    data_source_predictors = surface_samples_lat,
    var_name_pred = "lat",
    sel_method = "unconstrained",
    sel_complexity = "poly_2",
    transform_to_percentage = FALSE,
    tranformation = "none"
    )
composition_variation_surface_samples <- 
  dca_surface_samples$case_r %>% 
  dplyr::select(sample_id, axis_1, axis_2) %>% 
  dplyr::inner_join(surface_samples_lat, by = "sample_id") %>% 
  dplyr::mutate(axis_1 = max(axis_1) - axis_1 + min(axis_1)) %>% 
  dplyr::mutate(data_type = paste("surface_samples")) %>% 
  dplyr::group_by(lat) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(dataset_id = paste("dataset_", row.names(.), sep = "")) %>% 
  tidyr::unnest(data) %>% 
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    everything()
  ) 


dca_top_fossil_pollen_500yr <- 
  REcopol::fit_ordination(
    data_source_community = top_fossil_pollen_500yr_sp,
    data_source_predictors = top_fossil_pollen_500yr_lat,
    var_name_pred = "lat",
    sel_method = "unconstrained",
    sel_complexity = "poly_2",
    transform_to_percentage = FALSE,
    tranformation = "none"
  )

composition_variation_top_fossil_pollen_500yr <- 
  dca_top_fossil_pollen_500yr$case_r %>% 
  dplyr::select(sample_id, axis_1, axis_2) %>% 
  dplyr::inner_join(top_fossil_pollen_500yr_lat, by = "sample_id") %>% 
  dplyr::mutate(axis_1 = max(axis_1) - axis_1 + min(axis_1)) %>% 
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    everything()
  )  %>% 
  dplyr::mutate(data_type = paste("fossil_pollen_top_500yr"))


dca_combined <- 
  composition_variation_surface_samples %>% 
  dplyr::bind_rows(composition_variation_top_fossil_pollen_500yr) %>%  
  dplyr::select(
    dataset_id,
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5, 
    ecozone_koppen_15, 
    ecozone_koppen_30, 
    climate_zone_revised, 
    data_type,
    everything()
  ) 

write_rds(dca_combined,
          file = "Inputs/Data/dca_combined_300424.rds",
          compress = "gz")


