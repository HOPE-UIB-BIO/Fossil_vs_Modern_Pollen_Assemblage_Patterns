#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#
#          Modern vegetation data download from GBIF
#          
#----------------------------------------------------------#
library(tidyverse)
library(rgbif) # To download gbif data
library(taxize) # Useful for finding GBIF taxon IDs, etc.

#1. First go to GBIF webpage, enter the data filter criteria (e.g, lat-long ranges, elev,
# presence/absence, basis of record, target taxon or group, etc) manually. Aim here is
# to know the 'Taxon Key', and 'format of polygon' to be used to download data in R.

#2. After setting manually the filters, taxon key is displayed in the link of the web
# browser (e.g., https://www.gbif.org/occurrence/search?taxon_key=196&taxon_key=220&occurrence_status=present).
#3. Then click download tab, and the format od polygon to be specified is found there.

# Get these information and download the data using 'rgbif' package in r and exit the manual
# download or just agree to the terms and conditions of GBIF and manually download the data.

# Alternatively, Taxon Key can also be found out using 'taxize' package in R.
# e.g.,
get_gbifid(sci = 'Liliopsida')

# Download data of modern species occurrence from GBIF using R
taxon_keys <- c(220, 196) # 120 = Magnoliopsida, 196 = Liliopsida
download_data <- 
  rgbif::occ_download(
  rgbif::pred_within("POLYGON((75 25,125 25,125 66,75 66,75 25))"),
  pred_in("taxonKey", taxon_keys),
  pred_in("basisOfRecord", 
         c('PRESERVED_SPECIMEN',
           'HUMAN_OBSERVATION',
           'OBSERVATION',
           'MACHINE_OBSERVATION'
           )
         ),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  user = "kuber.bhatta1980",
 pwd = "kuber20bhatta37519",
 email = "kuber.bhatta@gmail.com",
 format = "SIMPLE_CSV"
 ) 

download_data 
# check 'Download key' and follow instructions given there.
# Citation: GBIF Occurrence Download https://doi.org/10.15468/dl.fvp4at Accessed from 
# R via rgbif (https://github.com/ropensci/rgbif) on 2024-01-16

# Check data processing status with
rgbif::occ_download_wait('0057625-231120084113126')
# Total records: 4486100

# After it finishes, use
raw_data_gbif <- 
  rgbif::occ_download_get(
    key = "0057625-231120084113126",
    path = "Inputs/Data/gbif_download/",
    overwrite = FALSE
    ) %>%
  rgbif::occ_download_import(.) # 4486100 records

write_rds(
  raw_data_gbif,
  file = "Inputs/Data/gbif_download/raw_data_gbif_160124.rds",
  compress = "gz"
  )

test_1 <- 
  read_rds("Inputs/Data/gbif_download/raw_data_gbif_160124.rds") %>% 
  as_tibble()

test_2 <- 
  test_1 %>% 
  dplyr::select(
    gbifID, 
    family, 
    class,
    scientificName,
    occurrenceStatus,
    decimalLatitude,
    decimalLongitude
    ) %>% 
  dplyr::mutate(
    lat = round(decimalLatitude, 2),
    long = round(decimalLongitude, 2)
    ) %>% 
  dplyr::mutate_at("gbifID", as.character) %>% 
  dplyr::filter(class == "Magnoliopsida" | class == "Liliopsida") %>% 
  dplyr::filter(!lat == "NA") %>% 
  dplyr::filter(!family == "NA") %>% 
  dplyr::select(
    -c(class, decimalLatitude, decimalLongitude)
    ) %>% 
  column_to_rownames("gbifID") %>% 
  distinct() %>% 
  rownames_to_column("gbifID") %>% 
  arrange(scientificName) %>% 
  dplyr::filter(!family == "Nymphaeaceae") %>% 
  dplyr::filter(!family == "Pontederiaceae") %>% 
  dplyr::filter(!family == "Najadaceae") %>% 
  dplyr::filter(!family == "Hydrocharitaceae") %>% 
  dplyr::filter(!family == "Lemnaceae") %>% 
  dplyr::filter(!family == "Alismataceae") #aquatic families removed
 
 
# Backbone tree
tree <- 
  ape::read.tree(
    paste(
      "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
      "Ramirez_Barahona_etal_2020_raxml_processed.tre",
      sep = ""
      )
    ) # 442 tip labels

taxa <- 
  test_2 %>% 
  dplyr::filter(!family == "") %>%  
  dplyr::select(family) %>% 
  distinct()        # 296 unique families

lacking_families <- 
  taxa %>% 
  dplyr::filter(!family %in% tree$tip.label) %>% 
  dplyr::filter(!is.na(family)) #24

for_correction <- 
  test_2 %>% 
  dplyr::filter(family %in% lacking_families$family) %>% 
  arrange(family) %>% 
  dplyr::mutate_at("gbifID", as.character) %>%
  dplyr::filter(!gbifID == "1917020115") # This ID lacks family (only order)

write_csv(for_correction,
          file = "Inputs/Tables/taxa_for_correction_gbif.csv")


# Families that were lacking in the backbone tree were corrected manually
  
# Corrected taxa
corrected_taxa <- 
  readr::read_csv("Inputs/Tables/taxa_gbif_corrected_160124.csv",
                  show_col_types = FALSE) %>% 
  dplyr::select(-family) %>% 
  dplyr::rename(family = family_corrected) %>% 
  dplyr::mutate_at("gbifID", as.character)

test_3 <- 
  test_2 %>% 
  as_tibble() %>% 
  dplyr::filter(!family == "") %>% 
  dplyr::filter(!gbifID %in% corrected_taxa$gbifID) %>% 
  bind_rows(corrected_taxa) %>% 
  dplyr::mutate(occurrence = ifelse(occurrenceStatus == "PRESENT", 1, 0)) %>% 
  dplyr::select(
    -c(scientificName,
       occurrenceStatus)
    ) 

lacking_families_1 <- 
  test_3 %>% 
  dplyr::filter(!family %in% tree$tip.label) # 0, All OK!

taxon_data <- 
  test_3 %>% 
  arrange(lat) %>% 
  dplyr::select(-long)
  
lat_data <- 
  taxon_data %>% 
  dplyr::select(gbifID, lat)

taxon_data_transform <- 
  taxon_data %>% 
  dplyr::group_by(lat, family) %>% 
  dplyr::summarise(
    .groups = "keep",
    occurrence = mean(occurrence)
    ) %>% 
  tidyr::spread(
    family, 
    occurrence
    ) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ tidyr::replace_na(., 0)
      )
    ) %>%
  dplyr::ungroup() %>% 
  dplyr::select_if(colSums(.) != 0) %>%
  tibble::as_tibble()

write_rds(
  taxon_data_transform,
  file = "Inputs/Data/data_gbif_for_phylodiversity_160124.rds",
  compress = "gz"
  )  

# TEST-----
# Randomly sample two samples from every 0.1 degree latitude
gbif <- 
  read_rds("Inputs/Data/data_gbif_for_phylodiversity_160124.rds") %>% 
  arrange(lat) %>% 
  dplyr::mutate(lat_grouped = ceiling(lat/0.1)*0.1) %>% 
  group_by(lat_grouped) %>% 
  slice_sample(n = 2) %>% 
  ungroup() %>% 
  dplyr::select(-lat_grouped) %>% 
  dplyr::select_if(colSums(.) != 0) %>%
  tibble::as_tibble()

write_rds(gbif,
          file = "Inputs/Data/gbif_test_data_310124.rds",
          compress = "gz")
