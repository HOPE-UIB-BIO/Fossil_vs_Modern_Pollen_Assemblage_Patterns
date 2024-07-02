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
#        Estimate phylogenetic relatedness (MPD, MNTD) ----
#
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# Load harmonised filtered data of surface samples ----
#-----------------------------------------------#
surface_pollen_filtered_clim_zone <- 
  read_rds("Inputs/Data/surface_pollen_filtered_030124.rds")

#-----------------------------------------------#
# Estimate MPD/MNTD of surface samples ----
#-----------------------------------------------#
set.seed(2012)

mpd_surface_samples <- 
  get_phylogenetic_diversity(
    counts = surface_pollen_filtered_clim_zone %>% 
      dplyr::select(
        -c(
          lat,
          long,
          ecozone_koppen_5, 
          ecozone_koppen_15, 
          ecozone_koppen_30,
          climate_zone_revised
          )
        ),
    backbone_tree = ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_surface_samples.tre",
        sep = ""
        )
      ),
    type = "mpd",
    null.model = "phylogeny.pool", 
    abundance.weighted = TRUE,
    runs = 999
    )

set.seed(2012)

mntd_surface_samples <- 
  get_phylogenetic_diversity(
    counts = surface_pollen_filtered_clim_zone %>% 
      dplyr::select(
        -c(
          lat,
          long,
          ecozone_koppen_5, 
          ecozone_koppen_15, 
          ecozone_koppen_30,
          climate_zone_revised
        )
      ),
    backbone_tree = ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_surface_samples.tre",
        sep = ""
      )
    ),
    type = "mntd",
    null.model = "phylogeny.pool", 
    abundance.weighted = TRUE,
    runs = 999
  )

phylodiversity_estimated_surface_samples <- 
  surface_pollen_filtered_clim_zone %>% 
  dplyr::select(
    sample_id,
    lat,
    long,
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised
    ) %>% 
  dplyr::inner_join(mpd_surface_samples %>% 
                      dplyr::select(-runs),
                    by = "sample_id") %>% 
  dplyr::inner_join(mntd_surface_samples %>% 
                     dplyr::select(-ntaxa),
                   by = "sample_id") %>% 
  drop_na()

write_rds(
  phylodiversity_estimated_surface_samples,
  file = "Inputs/Data/phylodiversity_estimated_surface_samples_050124.rds",
  compress = "gz"
  )

#-----------------------------------------------#
# Load harmonised filtered data of top fossil pollen samples ----
#-----------------------------------------------#
top_fossil_pollen_500yr_filtered_clim_zone <- 
  read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds")

#-----------------------------------------------#
# Estimate MPD/MNTD of surface samples ----
#-----------------------------------------------#
set.seed(2012)

mpd_top_500yr_samples <- 
  get_phylogenetic_diversity(
    counts = top_fossil_pollen_500yr_filtered_clim_zone %>% 
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
      ),
    backbone_tree = ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_top_fossil_samples.tre",
        sep = ""
      )
    ),
    type = "mpd",
    null.model = "phylogeny.pool", 
    abundance.weighted = TRUE,
    runs = 999
  )

set.seed(2012)

mntd_top_500yr_samples <- 
  get_phylogenetic_diversity(
    counts = top_fossil_pollen_500yr_filtered_clim_zone %>% 
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
      ),
    backbone_tree = ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_top_fossil_samples.tre",
        sep = ""
      )
    ),
    type = "mntd",
    null.model = "phylogeny.pool", 
    abundance.weighted = TRUE,
    runs = 999
  )

phylodiversity_estimated_top_500yr_samples <- 
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
    ) %>% 
  dplyr::inner_join(mpd_top_500yr_samples %>% 
                      dplyr::select(-runs),
                    by = "sample_id") %>% 
  dplyr::inner_join(mntd_top_500yr_samples %>% 
                      dplyr::select(-ntaxa),
                    by = "sample_id") %>% 
  drop_na()


readr::write_rds(
  phylodiversity_estimated_top_500yr_samples,
  file = "Inputs/Data/phylodiversity_estimated_top_500yr_samples_080124.rds",
  compress = "gz"
  )

# Combine all estimates
phylodiversity_estimated_surface_samples <- 
  readr::read_rds("Inputs/Data/phylodiversity_estimated_surface_samples_050124.rds") 
phylodiversity_estimated_top_500yr_samples <- 
  readr::read_rds("Inputs/Data/phylodiversity_estimated_top_500yr_samples_080124.rds")

names(phylodiversity_estimated_surface_samples)
names(phylodiversity_estimated_top_500yr_samples)

phylo_div_full <- 
  phylodiversity_estimated_surface_samples %>% 
  dplyr::mutate(data_type = paste("surface_samples")) %>% 
  dplyr::bind_rows(phylodiversity_estimated_top_500yr_samples %>% 
              dplyr::mutate(data_type = paste("top_500_yr"))
            ) %>% 
  dplyr::select(dataset_id, sample_id, lat, everything())

readr::write_rds(phylo_div_full,
          file = "Inputs/Data/phylo_div_full_090123.rds",
          compress = "gz")

