#----------------------------------------------------------#
# Fossil pollen data can predict robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#----------------------------------------------------------#

#----------------------------------------------------------#
#                  Harmonise data ----
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#

source("R/00_Config_file.R")

#-----------------------------------------------#

full_data <- 
  readr::read_rds(
    "Inputs/Data/data_combined_raw_121223.rds"
  )

dat_harmonised <-
  full_data %>%
  dplyr::mutate(
    harmonised_counts =
      purrr::map2(
        .x = sample_id, 
        .y = raw_counts,
        .f = ~ {
          harm_table <-
            readr::read_csv(
              "Inputs/Tables/surface_samples_harmonisation_table_221123.csv",
              show_col_types = FALSE
              )
          message(
            msg = paste0(.x)
          )
          
          # transform all taxa to numeric
          taxa_names <- names(.y)
           
          suppressWarnings(dat <-
                             .y %>%
                             dplyr::mutate_at(
                               taxa_names,
                               as.numeric
                             )
                           )
          res <- 
            dat %>%
            as.data.frame() %>%
            dplyr::mutate(
              dplyr::across(
                where(is.numeric),
                ~ tidyr::replace_na(., 0)
              )
            ) %>%
            tidyr::gather(
              key = "taxon_name",
              value = "counts"
              ) %>%
            dplyr::left_join(
              harm_table, 
              by = "taxon_name"
              ) %>%
            dplyr::filter(!family_level == "delete") %>%
            dplyr::group_by(
              family_level
            ) %>%
            dplyr::summarise(
              .groups = "keep",
              counts = sum(counts)
            ) %>%
            dplyr::rename(harmonised_counts = `family_level`) %>%
            tidyr::spread(
              harmonised_counts, 
              counts
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select_if(colSums(.) != 0) %>%
            tibble::as_tibble()
          return(res)
        }
      ),
    harmonised_percentages = 
      purrr::map2(
        .x = percentage,
        .y = harmonised_counts,
        .f = ~ {
          if(.x == "TRUE") {
            dat <- round(.y, digits = 3)
          } else {
            dat <- 
              round((.y/rowSums(.y))*100, 
                    digits = 3)
          }
        }
      )
  )

readr::write_rds(
  dat_harmonised,
  file = "Inputs/Data/data_harmonised_121223.rds",
  compress = "gz"
  )

