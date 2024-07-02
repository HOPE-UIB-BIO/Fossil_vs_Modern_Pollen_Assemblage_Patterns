#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#
#               Surface sample pollen data 
#          
#        Estimate phylogenetic relatedness (MPD, MNTD) ----
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

top_fossil_pollen_1000yr_filtered_clim_zone <- 
  read_rds("Inputs/Data/top_fossil_pollen_1000yr_filtered_080124.rds")

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


write_rds(
  phylodiversity_estimated_top_500yr_samples,
  file = "Inputs/Data/phylodiversity_estimated_top_500yr_samples_080124.rds",
  compress = "gz"
  )

# Top 1000-yrs samples
set.seed(2012)

mpd_top_1000yr_samples <- 
  get_phylogenetic_diversity(
    counts = top_fossil_pollen_1000yr_filtered_clim_zone %>% 
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

mntd_top_1000yr_samples <- 
  get_phylogenetic_diversity(
    counts = top_fossil_pollen_1000yr_filtered_clim_zone %>% 
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

phylodiversity_estimated_top_1000yr_samples <- 
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
    ) %>%  
  dplyr::inner_join(mpd_top_1000yr_samples %>% 
                      dplyr::select(-runs),
                    by = "sample_id") %>% 
  dplyr::inner_join(mntd_top_1000yr_samples %>% 
                      dplyr::select(-ntaxa),
                    by = "sample_id") %>% 
  drop_na()


write_rds(
  phylodiversity_estimated_top_1000yr_samples,
  file = "Inputs/Data/phylodiversity_estimated_top_1000yr_samples_080124.rds",
  compress = "gz"
  )

#-----------------------------------------------#
# Load modern vegetation data of taxa occurrence ----
#-----------------------------------------------#
taxon_data_transform <- 
  read_rds("Inputs/Data/data_gbif_for_phylodiversity_160124.rds") %>% 
  dplyr::mutate(
    sample_id = paste(
      "sample_", row.names(.), 
      sep = ""
      )
    ) %>% 
  dplyr::select(sample_id, lat, everything())
modern_data_for_analysis <- 
  taxon_data_transform %>% 
  column_to_rownames("sample_id") %>% 
  dplyr::select(-lat)

modern_data_for_analysis <- 
  modern_data_for_analysis[rowSums(modern_data_for_analysis > 0) > 2,] %>% # filter out samples with < 3 taxa
  rownames_to_column("sample_id")

#-----------------------------------------------#
# Estimate MPD/MNTD of modern vegetation samples ----
#-----------------------------------------------#
set.seed(2012)

mpd_modern_vegetation <- 
  get_phylogenetic_diversity(
    counts = modern_data_for_analysis,
    backbone_tree = ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_modern_vegetation.tre",
        sep = ""
      )
    ),
    type = "mpd",
    null.model = "phylogeny.pool", 
    abundance.weighted = FALSE,
    runs = 999
  )

set.seed(2012)

mntd_modern_vegetation <- 
  get_phylogenetic_diversity(
    counts = modern_data_for_analysis,
    backbone_tree = ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_modern_vegetation.tre",
        sep = ""
      )
    ),
    type = "mntd",
    null.model = "phylogeny.pool", 
    abundance.weighted = FALSE,
    runs = 999
  )

phylodiversity_estimated_modern_vegetation <- 
  taxon_data_transform %>% 
  dplyr::select(sample_id, lat) %>% 
  dplyr::inner_join(mpd_modern_vegetation %>% 
                      dplyr::select(-runs),
                    by = "sample_id") %>% 
  dplyr::inner_join(mntd_modern_vegetation %>% 
                      dplyr::select(-ntaxa),
                    by = "sample_id") %>% 
  drop_na()


write_rds(
  phylodiversity_estimated_modern_vegetation,
  file = "Inputs/Data/phylodiversity_estimated_modern_vegetation_170124.rds",
  compress = "gz"
  )

# Combine all estimates
phylodiversity_estimated_surface_samples <- 
  read_rds("Inputs/Data/phylodiversity_estimated_surface_samples_050124.rds") 
phylodiversity_estimated_top_500yr_samples <- 
  read_rds("Inputs/Data/phylodiversity_estimated_top_500yr_samples_080124.rds")
phylodiversity_estimated_top_1000yr_samples <- 
  read_rds("Inputs/Data/phylodiversity_estimated_top_1000yr_samples_080124.rds") 
phylodiversity_estimated_modern_vegetation <- 
  read_rds("Inputs/Data/phylodiversity_estimated_modern_vegetation_170124.rds")


names(phylodiversity_estimated_surface_samples)
names(phylodiversity_estimated_top_500yr_samples)
names(phylodiversity_estimated_top_1000yr_samples)
names(phylodiversity_estimated_modern_vegetation)

phylo_div_full <- 
  phylodiversity_estimated_surface_samples %>% 
  dplyr::mutate(data_type = paste("surface_samples")) %>% 
  bind_rows(phylodiversity_estimated_top_500yr_samples %>% 
              dplyr::mutate(data_type = paste("top_500_yr"))
            ) %>% 
  bind_rows(phylodiversity_estimated_top_1000yr_samples %>% 
              dplyr::mutate(data_type = paste("top_1000_yr"))
            ) %>% 
  bind_rows(phylodiversity_estimated_modern_vegetation %>% 
              dplyr::mutate(data_type = paste("modern_vegetaton"))
            ) %>% 
  dplyr::select(dataset_id, sample_id, lat, everything())

write_rds(phylo_div_full,
          file = "Inputs/Data/phylo_div_full_090123.rds",
          compress = "gz")

