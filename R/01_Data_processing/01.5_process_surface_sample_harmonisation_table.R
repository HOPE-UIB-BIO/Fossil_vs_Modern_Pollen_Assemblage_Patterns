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
# Prepapre surface pollen data harmonisation table ----
#                          
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

pangaea_taxa <- 
  readr::read_csv("Inputs/Tables/taxa_names_surface_samples_pangaea_150923.csv") %>% 
  dplyr::distinct(taxon_name)
empd_taxa <- 
  readr::read_csv("Inputs/Tables/taxa_names_empd_211123.csv") %>% 
  dplyr::distinct()
cao_taxa <- 
  readr::read_csv("Inputs/Tables/taxa_names_cao_data_221123.csv")


full_list <- 
  empd_taxa %>% 
  dplyr::full_join(cao_taxa,
                   by = "taxon_name") %>% 
  dplyr::full_join(pangaea_taxa,
                   by = "taxon_name") %>% 
  dplyr::mutate(family_level = NA) %>% 
  dplyr::arrange(taxon_name)

# readr::write_csv(full_list,
#          file = "Inputs/Tables/surface_samples_harmonisation_table_221123.csv")

# Note: this table is updated to created full harmonisation table!
