#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
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
  read_rds("Inputs/Data/surface_pollen_filtered_030124.rds") %>% 
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
  enframe(name = NULL,
          value = "taxa") %>%
  distinct() %>% 
  nrow(.)

n_taxa_top_fossil_pollen_500yr <- 
  read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds") %>% 
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
  enframe(name = NULL,
          value = "taxa") %>%
  distinct() %>% 
  nrow(.)

n_taxa_top_fossil_pollen_1000yr <- 
  read_rds("Inputs/Data/top_fossil_pollen_1000yr_filtered_080124.rds") %>% 
  dplyr::select(
    -c(
      dataset_id,
      sample_id,
      lat,
      long,
      ecozone_koppen_5,
      ecozone_koppen_15,
      ecozone_koppen_30,
      climate_zone_revised
    )
    ) %>% 
  colnames(.) %>%
  enframe(name = NULL,
          value = "taxa") %>%
  distinct() %>% 
  nrow(.)


phylo_div_full <- read_rds("Inputs/Data/phylo_div_full_090123.rds")
turnover_combined <- read_rds("Inputs/Data/turnover_combined_050224.rds")

surface_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "surface_samples") %>% 
  summarise(datasets = nrow(.),
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
                  distinct() %>% 
                  nrow(.)
  ) %>% 
  summarise(
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

top_1000_yr_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "top_1000_yr") %>% 
  dplyr::mutate(datasets = phylo_div_full %>% 
                  dplyr::filter(data_type == "top_1000_yr") %>% 
                  dplyr::select(dataset_id) %>% 
                  distinct() %>% 
                  nrow(.)
  ) %>% 
  summarise(
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
  dplyr::mutate(n_taxa = paste(n_taxa_top_fossil_pollen_1000yr))

phylodiversity <- 
  bind_rows(surface_samples,
            top_500_yr_samples,
            top_1000_yr_samples) %>% 
  as_tibble()
compositional_turnover <- 
  turnover_combined %>% 
  dplyr::mutate(data_type = ifelse(data_type == "fossil_pollen_top_500yr",
                                   "top_500_yr", data_type),
                data_type = ifelse(data_type == "fossil_pollen_top_1000yr",
                                   "top_1000_yr", data_type)
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
  as_tibble()


output_summary <- 
  phylodiversity %>% 
  dplyr::left_join(compositional_turnover, by = "data_type")

write_csv(output_summary,
          file = "Outputs/Table/Data_outputs_summary_050224.csv")
