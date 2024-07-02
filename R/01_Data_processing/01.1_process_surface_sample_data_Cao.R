#----------------------------------------------------------#
# Fossil pollen data can predict robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#                     
#----------------------------------------------------------#

#----------------------------------------------------------#
#
#         Surface sample pollen data from Cao ----
#                          
#----------------------------------------------------------#

#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

data_cao_samples <- 
  readr::read_csv(
    "Inputs/Data/surface_samples_cao/mp.csv"
  ) %>% 
  dplyr::rename(sample_id = No.) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(raw_counts = -group_cols()) %>% 
  dplyr::ungroup()
  
data_cao_lat_long <- 
  read_csv(
    "Inputs/Data/surface_samples_cao/mc.csv"
  ) %>% 
  dplyr::select(sample_id = No.,
                sitename = `Site name`,
                long,
                lat,
                Sedtype) %>% 
  dplyr::mutate(percentage = "TRUE",
                depositional_env = 
                  ifelse(Sedtype == "1", "Lake sediment surface", Sedtype),
                depositional_env = 
                  ifelse(Sedtype == "2", "Moss polster", depositional_env),
                depositional_env = 
                  ifelse(Sedtype == "3", "Snow, ice and glacier", depositional_env),
                depositional_env = 
                  ifelse(Sedtype == "4", "Soil surface", depositional_env),
                depositional_env = 
                  ifelse(Sedtype == "5", "Marine sediment surface", depositional_env),
                depositional_env = 
                  ifelse(Sedtype == "6", "Dust flux", depositional_env)
                ) %>% 
  dplyr::select(-Sedtype) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(metadata = -group_cols()) %>% 
  dplyr::ungroup()

data_cao <- 
  inner_join(data_cao_lat_long,
             data_cao_samples,
             by = "sample_id")
write_rds(data_cao,
          file = "Inputs/Data/surface_samples_cao/surface_samples_cao_221123.rds",
          compress = "gz")

taxa_cao <-
  data_cao %>%
  dplyr::mutate(taxa =
                  purrr::map(raw_counts,
                             ~ colnames(.x) %>%
                               enframe(name = NULL,
                                       value = "taxon_name")
                             )
                ) %>%
  dplyr::select(taxa) %>%
  tidyr::unnest(taxa) %>% 
  distinct()
write_csv(taxa_cao,
          file = "Inputs/Tables/taxa_names_cao_data_221123.csv")
