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
#                 Data output summary ----
#
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source configuration ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# Extract_summary
n_taxa_surface_samples <- 
  readr::read_rds("Inputs/Data/surface_pollen_filtered_030124.rds") %>% 
  dplyr::select(
    -c(sample_id, 
       lat, 
       long, 
       ecozone_koppen_5,
       ecozone_koppen_15,
       ecozone_koppen_30,
       climate_zone_revised)
    ) %>% 
  colnames(.) %>%
  tibble::enframe(name = NULL,
          value = "taxa") %>%
  dplyr::distinct() %>% 
  nrow(.)

n_taxa_top_fossil_pollen_500yr <- 
  readr::read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds") %>% 
  dplyr::select(
    -c(dataset_id, 
       sample_id, 
       lat,
       long,
       ecozone_koppen_5,
       ecozone_koppen_15,
       ecozone_koppen_30,
       climate_zone_revised)
    ) %>% 
  colnames(.) %>%
  tibble::enframe(name = NULL,
          value = "taxa") %>%
  dplyr::distinct() %>% 
  nrow(.)

phylo_div_full <- readr::read_rds("Inputs/Data/phylo_div_full_090123.rds")
turnover_combined <- readr::read_rds("Inputs/Data/turnover_combined_050224.rds")

surface_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "surface_samples") %>% 
  dplyr::summarise(datasets = nrow(.),
            samples = nrow(.),
            mean_taxon_richness = mean(ntaxa),
            min_mpd = min(ses_mpd),
            max_mpd = max(ses_mpd),
            mean_mpd = mean(ses_mpd),
            min_mntd = min(ses_mntd),
            max_mntd = max(ses_mntd),
            mean_mntd = mean(ses_mntd),
            data_type = unique(.$data_type)
            ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(n_taxa = paste(n_taxa_surface_samples))

top_500_yr_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "top_500_yr") %>% 
  dplyr::mutate(datasets = phylo_div_full %>% 
                  dplyr::filter(data_type == "top_500_yr") %>% 
                  dplyr::select(dataset_id) %>% 
                  dplyr::distinct() %>% 
                  nrow(.)
                ) %>% 
  dplyr::summarise(
    datasets = mean(datasets),
    samples = nrow(.),
    mean_taxon_richness = mean(ntaxa),
    min_mpd = min(ses_mpd),
    max_mpd = max(ses_mpd),
    mean_mpd = mean(ses_mpd),
    min_mntd = min(ses_mntd),
    max_mntd = max(ses_mntd),
    mean_mntd = mean(ses_mntd),
    data_type = unique(.$data_type)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_taxa = paste(n_taxa_top_fossil_pollen_500yr))

phylodiversity <- 
  dplyr::bind_rows(surface_samples,
            top_500_yr_samples) %>% 
  tibble::as_tibble()

compositional_turnover <- 
  turnover_combined %>% 
  dplyr::mutate(data_type = ifelse(data_type == "fossil_pollen_top_500yr",
                                   "top_500_yr", data_type)
                ) %>% 
  dplyr::select(data_type, axis_1) %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(total_gardient_length = 
                  purrr::map_dbl(data,
                                 ~ max(.x$axis_1) - min(.x$axis_1)
                                 ),
                min_turnover = purrr::map_dbl(data,
                                              ~ min(.x$axis_1)
                                              ),
                mean_turnover = purrr::map_dbl(data,
                                              ~ mean(.x$axis_1)
                                              )
                ) %>% 
  dplyr::select(-data) %>% 
  tibble::as_tibble()


output_summary <- 
  phylodiversity %>% 
  dplyr::left_join(compositional_turnover, by = "data_type")

readr::write_csv(output_summary,
                 file = "Outputs/Table/Data_outputs_summary_050224.csv")
