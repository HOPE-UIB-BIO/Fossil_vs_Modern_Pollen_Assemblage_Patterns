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
#         Surface sample pollen data from EPD ----
#                          
#----------------------------------------------------------#

#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")


# Surface samples of Eurasian region are available at "https://doi.pangaea.de/10.1594/PANGAEA.909130"

read_excel_allsheets <- 
    function(file_path, 
             tibble = TRUE) {
    sheets <- readxl::excel_sheets(file_path)
    sheet_list <- 
      purrr::map(sheets,
                 ~ readxl::read_excel(file_path, 
                                      sheet = .x)
                 )
    names(sheet_list) <- sheets
    sheet_list
    }
  
mysheets <- 
  read_excel_allsheets("Inputs/Data/surface_samples_epd/EMPD2_count_v1.xlsx")
# Ignore warnings!  

metadata <- 
  mysheets$metadata %>% 
  dplyr::select(sample_id = SampleName,
                sitename = SiteName,
                lat = Latitude,
                long = Longitude,
                percentage = ispercent,
                depositional_env = SampleType) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>% 
  tidyr::nest(.key = "metadata")

counts <- 
  mysheets$counts %>% 
  dplyr::select(
    sample_id = SampleName,
    taxon = original_varname,
    raw_counts = count
  ) %>%
 tidyr::pivot_wider(
    names_from = taxon,
    values_from = raw_counts,
    values_fill = 0 # empty values as '0', otherwise turns NA
  ) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(.key = "raw_counts")


data_empd <- 
  dplyr::inner_join(metadata,
             counts,
             by = "sample_id")

readr::write_rds(data_empd,
                 file = "Inputs/Data/surface_samples_epd/data_empd_211123.rds",
                 compress = "gz")

empd_taxa <- 
  data_empd %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(taxa = 
                  purrr::map(
                    raw_counts,
                    ~  colnames(.x) %>% 
                      tibble::enframe(name = NULL,
                                      value = "taxa")
                    )
                ) %>% 
  dplyr::select(taxa) %>% 
  tidyr::unnest(taxa) %>% 
  dplyr::distinct()
  
  

harmonisation_table <- 
  mysheets$p_vars %>% 
  dplyr::select(taxon_name = original_varname,
                eco_group = groupid) %>% 
  dplyr::filter(taxon_name %in% empd_taxa$taxa) 

write_csv(harmonisation_table,
          file = "Inputs/Tables/taxa_names_empd_211123.csv")
