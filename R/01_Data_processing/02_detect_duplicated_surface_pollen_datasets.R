#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#
#               Surface sample pollen data 
#          
#               Detect duplicated datasets ----
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# Load raw data ----
#-----------------------------------------------#
epd <- 
  read_rds(
    "Inputs/Data/surface_samples_epd/data_empd_211123.rds"
  ) %>% 
  tidyr::unnest(metadata) %>% 
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>% 
  dplyr::mutate_at("sample_id", as.character) %>% 
  dplyr::mutate_at("percentage", as.character)

cao <- 
  read_rds(
    "Inputs/Data/surface_samples_cao/surface_samples_cao_221123.rds"
  ) %>% 
  tidyr::unnest(metadata) %>% 
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>% 
  dplyr::mutate_at("sample_id", as.character) %>% 
  dplyr::mutate_at("percentage", as.character)

pg <- 
  read_rds(
    "Inputs/Data/surface_samples_cao/surface_samples_pangaea_220923.rds"
  ) %>% 
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>% 
  dplyr::mutate_at("sample_id", as.character) %>% 
  dplyr::mutate_at("percentage", as.character)


dup_epd_cao_lat <- 
  epd %>%
  dplyr::filter(lat %in% cao$lat)

dup_epd_cao_long <- 
  epd %>%
  dplyr::filter(long %in% cao$long)

# None are duplicated!

dup_epd_pg_lat <- 
  pg %>%
  dplyr::filter(lat %in% epd$lat)

dup_epd_pg_long <- 
  pg %>%
  dplyr::filter(long %in% epd$long)


dup_cao_pg_lat <- 
  cao %>%
  dplyr::filter(lat %in% pg$lat)

dup_cao_pg_long <- 
  cao %>%
  dplyr::filter(long %in% pg$long)
# None are duplicated!

exclude <- 
  c(
    "TMo-CF_001", "TMo-CF_002", "TMo-CF_003", "TMo-CF_004", "TMo-CF_005", 
    "TMo-CF_006", "TMo-CF_007", "TMo-CF_008", "TMo-CF_009", "TMo-CF_010", 
    "TMo-CF_011", "TMo-CF_012", "TMo-CF_013",  "TMo-CF_014", "TMo-CF_015",
    "TMo-CF_016", "MRUT1M04", "07-NE-12", "07-NE-13", "07-NE-14", "07-NE-15",
    "07-WI-07", "07-WI-08", "07-WI-09", "07-WI-10", "07-WI-11", "07-WI-12"
    )

pg_1 <- 
  pg %>% 
  dplyr::filter(!sample_id  %in% exclude) %>% 
  dplyr::select(
    -c(
      publication,
      counts
      )
    ) %>% 
  dplyr::mutate(
    raw_counts = purrr::map(raw_counts,
                            ~ .x %>% 
                              mutate_if(is.double, as.numeric)
                            )
    ) %>% 
  tidyr::unnest(raw_counts)

# 131 dataframe contains Lamiaceae as character, this data has too low parecentages

pg_1 <- 
  pg %>% 
  dplyr::filter(!sample_id  %in% exclude) %>% 
  dplyr::select(
    -c(
      publication,
      counts
    )
  ) %>% 
  dplyr::mutate(
    raw_counts = purrr::map(raw_counts,
                            .f = ~ {
                              col_names <- colnames(.x)
                              suppressWarnings(
                                dat <-
                                  .x %>%
                                  dplyr::mutate_at(
                                    col_names,
                                    as.numeric
                                    )
                                )
                              }
                            )
                 ) %>% 
  dplyr::select(-dataset_id)


names(pg_1)
names(epd)
names(cao)


full_data <- 
  bind_rows(pg_1, epd, cao)

write_rds(full_data,
          file = "Inputs/Data/data_combined_raw_121223.rds",
          compress = "gz")
