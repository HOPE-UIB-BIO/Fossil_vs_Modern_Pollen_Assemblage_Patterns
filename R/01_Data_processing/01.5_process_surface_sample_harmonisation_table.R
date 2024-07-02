#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#
#         Surface sample pollen data harmonisation table ----
#                          
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

pangaea_taxa <- 
  read_csv("Inputs/Tables/taxa_names_surface_samples_pangaea_150923.csv") %>% 
  distinct(taxon_name)
empd_taxa <- 
  read_csv("Inputs/Tables/taxa_names_empd_211123.csv") %>% 
  distinct()
cao_taxa <- 
  read_csv("Inputs/Tables/taxa_names_cao_data_221123.csv")


full_list <- 
  empd_taxa %>% 
  dplyr::full_join(cao_taxa,
                   by = "taxon_name") %>% 
  dplyr::full_join(pangaea_taxa,
                   by = "taxon_name") %>% 
  dplyr::mutate(family_level = NA) %>% 
  dplyr::arrange(taxon_name)

write_csv(full_list,
          file = "Inputs/Tables/surface_samples_harmonisation_table_221123.csv")
