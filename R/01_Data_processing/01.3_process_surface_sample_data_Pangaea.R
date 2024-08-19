#----------------------------------------------------------#
# Fossil pollen data can reconstruct robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#----------------------------------------------------------#

#----------------------------------------------------------#
#
#         Surface sample pollen data from Pangaea ----
#                          
#----------------------------------------------------------#

#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

# Download surface sample data (14.09.2023)
surface_sample1 <- 
  pangaear::pg_search(query = 
                        "surface pollen", 
                      bbox = c(75, 25, 125, 66)
                      ) %>%   #min_long, min_lat, max_long, max_lat
  dplyr::select(doi, citation)

surface_sample2 <- 
  pangaear::pg_search(query = 
                        "surface sample", 
                      bbox = c(75, 25, 125, 66)
                      ) %>% 
  dplyr::select(doi, citation)

surface_sample3 <- 
  pangaear::pg_search(query = 
                        "modern pollen", 
                      bbox = c(75, 25, 125, 66)
                      ) %>% 
  dplyr::select(doi, citation)

doi_tibble <- 
  bind_rows(surface_sample1, surface_sample2, surface_sample3) %>% 
  distinct()

# Checked all DOIs one-by-one in pangaea.de, following are the potential datasets
req_doi <- 
  c(
    "10.1594/PANGAEA.849666",
    "10.1594/PANGAEA.808958",
    "10.1594/PANGAEA.808953",
    "10.1594/PANGAEA.829753",
    "10.1594/PANGAEA.871524",
    "10.1594/PANGAEA.933664",
    "10.1594/PANGAEA.717143",
    "10.1594/PANGAEA.940115",
    "10.1594/PANGAEA.909130",
    "10.1594/PANGAEA.910726"
  )

raw_data <- 
  purrr::map(
    .x = req_doi, 
    .f = function(.x) {
      raw_dat <- pangaear::pg_data(.x)
      proc_dat <- 
        purrr::map(
          .x = raw_dat,
          .f = function(.x) {
            doi <- .x$doi
            citation <- .x$metadata$citation 
            data <- list(.x$data)
            metadata <- list(.x$metadata)
            result <- 
              tibble(
                doi,
                citation,
                data,
                metadata)
            return(result)
            }
          ) 
      }
    ) %>% 
  unlist(., recursive = FALSE) %>% 
  dplyr::bind_rows()

# Individually checked DOIs
sel_doi <- c(
  "10.1594/PANGAEA.849666",
  "10.1594/PANGAEA.808952",
  "10.1594/PANGAEA.808953",
  "10.1594/PANGAEA.808953",
  "10.1594/PANGAEA.829749",
  "10.1594/PANGAEA.829750",
  "10.1594/PANGAEA.829751",
  "10.1594/PANGAEA.871523",
  "10.1594/PANGAEA.871521",
  "10.1594/PANGAEA.933659",
  "10.1594/PANGAEA.717143",
  "10.1594/PANGAEA.910726"
 )
raw_data_1 <- 
  raw_data %>% 
  dplyr::filter(
    doi %in% sel_doi
  )
readr::write_rds(
  raw_data_1,
  file = "Inputs/Data/data_raw_surface_samples_140923.rds",
  compress = "gz"
  )  


# Extract taxa names
raw_taxa <- 
  purrr::map(raw_data_1$metadata, ~ .x$parameters) %>% 
  unlist(., recursive = FALSE) %>% 
  purrr::map(., ~ as.data.frame(.x , sep = "\t") %>% 
               slice(1,)
             ) %>% 
  dplyr::bind_rows() %>% 
  dplyr::rename(taxa = .x) %>% 
  dplyr::distinct() 

clean_names <- 
  purrr::map(raw_taxa$taxa, 
             ~ str_replace(as.character(.x), "[%]", "#")) %>% 
  do.call(rbind.data.frame, .) %>% 
  dplyr::distinct() %>% 
  purrr::map(., ~ strsplit(as.character(.x), " [#] (", fixed = TRUE)) %>% 
  unlist(., recursive = FALSE) %>%
  do.call(rbind.data.frame, .)
names(clean_names) <- as.character(c('taxon.name', 'abbrev'))

clean_names <- 
  clean_names %>% 
  dplyr::mutate(
    abbrev = gsub(".*\\(","", abbrev),
    abbrev = gsub("\\).*","", abbrev)
    ) %>% 
  dplyr::distinct(taxon.name, abbrev, .keep_all = TRUE)

#readr::write_csv(
#  clean_names,
#  file = "Inputs/Tables/taxa_names_surface_samples_pangaea_150923.csv"
#  )  


# Process individual dataset
raw_dat <- 
  readr::read_rds(
    "Inputs/Data/data_raw_surface_samples_140923.rds"
    ) %>% 
  dplyr::distinct()

#10.1594/PANGAEA.849666 
dt1 <- 
  raw_dat %>% 
  dplyr::filter(
    doi == "10.1594/PANGAEA.849666"
   )
  
dataset_id_849666 <- "PANGAEA.849666"
publication_849666 <- dt1$citation
depositional_env_849666 <- "Moss_soil_and_lake"
counts_849666_raw <- 
  dt1$data[[1]][,-64] %>% # bot = botryococcus is duplicated with bot = botrychium
  dplyr::select(
    sample_id = Event,
    lat = Latitude,
    long = Longitude,
    everything(),
    -c(
      "Elevation [m]",
      "Depth [m]"
      )  
    ) 
lat_long_849666 <- 
  counts_849666_raw %>% 
  dplyr::select(sample_id,
                lat, 
                long) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(lat_long = -group_cols()) %>% 
  dplyr::ungroup()

counts_849666 <- 
  counts_849666_raw %>% 
  dplyr::select(
    -c(lat, long)
    ) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(counts = -group_cols()) %>% 
  dplyr::ungroup()

lat_long_sample_id_849666 <- 
  dplyr::inner_join(
    counts_849666, 
    lat_long_849666,
    by = "sample_id")


dat1_849666 <- 
  tibble(
    dataset_id = dataset_id_849666, 
    depositional_env = depositional_env_849666,
    lat_long_sample_id_849666,
    publication = publication_849666
    ) %>% 
  dplyr::mutate(percentage = FALSE) %>% 
  tidyr::unnest(lat_long) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)


#10.1594/PANGAEA.808952
dt2 <- 
  raw_dat %>% 
  dplyr::filter(
    doi == "10.1594/PANGAEA.808952"
  )
dt2$citation
as <- dt2$data[[1]] #fossil pollen data: EXCLUDE


# 10.1594/PANGAEA.808953
dt3 <- 
  raw_dat %>% 
  dplyr::filter(
    doi == "10.1594/PANGAEA.808953"
  )
dataset_id_808953 <- "PANGAEA.808953"
publication_808953 <- dt3$citation[1]
depositional_env_808953 <- "lake"

#"Urtae [#]" at locations 90 and 97."Poac [#]" at locations 111 and 112 are duplicated names
#"Urtae [#]" is urticaceae at both locations and "Poac [#]" is Poaceae, verified from Parameters

Urticaceae <- dt3$data[[1]][,96] + dt3$data[[1]][,103]
Poaceae <-  dt3$data[[1]][,117] + dt3$data[[1]][,118]

counts_808953 <- 
  dt3$data[[1]][,-c(96, 103, 117, 118)] %>%  
  dplyr::select(
    sample_id = Event,
    everything(),
    -c(
      "Sample label",
      "Analyst",
      "No",
      "Ratio",
      "AP [%]",
      "Pollen indet [#]" 
      )
    ) %>% 
  dplyr::bind_cols(
    Urticaceae, 
    Poaceae
    ) %>% 
  dplyr::slice(n = 1:8) %>%  # bottom 3 rows are NA
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(
    counts = -group_cols()
    ) %>% 
  dplyr::ungroup()

event <- 
  dt3[1,]$metadata[[1]]$events %>% 
  stringr::str_split(., "; ") 

lat_long_808953 <-
  tibble(
    sample_id = c(
      "ML_CF_023",
      "ML_CF_018",
      "Tso_Moriri",
      "ML_CF_011",
      "ML_CF_007",
      "ML_CF_005a",
      "ML_CF_003a",
      "ML_CF_001"
    ),
    lat = c(
      33.766000,
      33.251000,
      32.929440,
      32.762000,
      32.557000,
      32.364000,
      32.346000,
      32.317000
    ),
    long = c(
      77.764000,
      77.841000,
      78.323330,
      77.400000,
      76.995000,
      77.249000,
      77.222000,
      77.187000
    )
  ) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(lat_long = -group_cols()) %>% 
  dplyr::ungroup()

lat_count_808953 <- 
  dplyr::inner_join(
    lat_long_808953,
    counts_808953,
    by = "sample_id"
    )

dat_808953 <- 
  tibble(
    dataset_id = dataset_id_808953, 
    depositional_env = depositional_env_808953,
    lat_count_808953,
    publication = publication_808953
  ) %>% 
  dplyr::mutate(percentage = FALSE) %>% 
  tidyr::unnest(lat_long) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)
  

# 10.1594/PANGAEA.829749
dt4 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.829749")
as <- dt4$data[[1]]  # NPP: Not quitable


# 10.1594/PANGAEA.829750 
dt5 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.829750")
dt5$data # Fossil pollen data: NOT suitable


# 10.1594/PANGAEA.829751
dt6 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.829751")

dataset_id_829751 <- "PANGAEA.829751"
depositional_env_829751 <- "lake"
publication_829751 <- dt6$citation
lat_long_829751 <-
  dt6$data[[1]] %>%
  dplyr::select(sample_id = Event,
                lat = Latitude,
                long = Longitude) %>%
  dplyr::group_by(sample_id) %>%
  tidyr::nest(
    lat_long = -group_cols()
    ) %>%
  dplyr::ungroup()

counts_829751 <- 
  dt6$data[[1]] %>%
  dplyr::select(
    -c(
      Latitude,
      Longitude,
      `Depth [m]`,
      `Mark found [#]`
      ),
    sample_id = Event,
    everything()
    ) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(
    counts = -group_cols()
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::inner_join(lat_long_829751,
             by = "sample_id")
dat_829751 <- 
  tibble(
    dataset_id = dataset_id_829751, 
    depositional_env = depositional_env_829751,
    counts_829751,
    publication = publication_829751
  ) %>% 
  dplyr::mutate(percentage = FALSE) %>% 
  tidyr::unnest(lat_long) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)


#10.1594/PANGAEA.871523
dt7 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.871523")

dataset_id_871523 <- "PANGAEA.871523"
depositional_env_871523 <- 
  dt7$data[[1]] %>% 
  dplyr::mutate(depositional_env = 
                  ifelse(
                    stringr::str_detect(Event, "River") == TRUE,
                    "River",
                    "Marine"
                    )
                ) %>% 
  dplyr::select(sample_id = Event,
                depositional_env)
  
publication_871523 <- dt7$citation

lat_long_871523 <- 
  data.frame(
    strsplit(
      as.character(dt7$metadata[[1]]$events), 
      '; ',  
      fixed = TRUE
      )
    ) %>% 
  rlang::set_names("V") %>% 
  dplyr::mutate(
    V1 = purrr::map(V, 
                    function(.x){
                      dt <- data.frame(
                        strsplit(
                          as.character(.x), 
                          ' * ',  
                          fixed = TRUE
                          )
                        ) %>% 
                        rlang::set_names("V1")
        dt1 <- 
          data.frame(
            strsplit(
              as.character(dt$V1), 
              ': ', 
              fixed = TRUE
              )
            ) %>% 
          `colnames<-`(.[1, ]) %>%  # Mame first row as column name
          .[-1, ] 
        
        colnames(dt1)[1] <- "sample_id"
        return(dt1)
        }
      )
    ) %>% 
  dplyr::select(V1) %>% 
  tidyr::unnest(V1) %>% 
  dplyr::select(
    sample_id,
    lat = LATITUDE,
    long = LONGITUDE
  ) 
counts_871523 <- 
  dt7$data[[1]] %>%
  dplyr::select(
    -c(
      `Depth [m]`,
      `Depth top [m]`,
      `Depth bot [m]`,
      `Lyc-T [#]`,
      `Samp m [g]`,
      `Sample ID`,
      Texture
    ),
    sample_id = Event,
    everything() 
  )

which(duplicated(counts_871523$sample_id)) #85, "Liaodong_Bay-SP73"
test <- counts_871523 %>% 
  dplyr::filter(sample_id == "Liaodong_Bay-SP73") %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise_if(is.numeric, sum)

counts_871523a <- 
  counts_871523 %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise_if(is.numeric, sum) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(
    counts = -group_cols()
    ) %>% 
  dplyr::ungroup()

comb_dat_871523 <- 
  dplyr::inner_join(
    depositional_env_871523,
    counts_871523a,
    by = "sample_id"
    ) %>% 
  dplyr::distinct() %>% 
  dplyr::inner_join(
    lat_long_871523,
    by = "sample_id"
    )

dat_871523 <- 
  tibble(
    dataset_id = dataset_id_871523, 
    comb_dat_871523,
    publication = publication_871523
    ) %>% 
  dplyr::mutate(percentage = FALSE) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)


#10.1594/PANGAEA.871523
dt8 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.871521") # duplicated, exclude


#10.1594/PANGAEA.933659
dt9 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.933659")

dataset_id_933659 <- "PANGAEA.933659"
publication_933659 <- dt9$citation
lat_long_933659 <- 
  dt9$data[[1]] %>% 
  dplyr::select(
    sample_id = Event,
    lat = Latitude,
    long = Longitude
                )
depo_env_933659 <- 
  data.frame(
    strsplit(dt9$metadata[[1]]$events,
             '; ',  
             fixed = TRUE
             )
    ) %>% 
  rlang::set_names("V") %>% 
  dplyr::mutate(
    V1 = purrr::map(V, 
                    function(.x){
                      dt <- data.frame(
                        strsplit(
                          as.character(.x), 
                          ' * ',  
                          fixed = TRUE
                        )
                      ) %>% 
                        rlang::set_names("V1")
                      dt1 <- 
                        data.frame(
                          strsplit(
                            as.character(dt$V1), 
                            ': ', 
                            fixed = TRUE
                          )
                        ) %>% 
                        `colnames<-`(.[1, ]) %>%  # Mame first row as column name
                        .[-1, ] 
                      
                      colnames(dt1)[1] <- "sample_id"
                      return(dt1)
                    }
    )
  ) %>% 
  dplyr::select(V1) %>% 
  tidyr::unnest(V1)

depositional_env_933659 <- "soil_moss"

which(duplicated(dt9$data[[1]]$Event))
counts_933659 <- 
  dt9$data[[1]] %>% 
  dplyr::select(
    -c(
      Latitude,
      Longitude,
      `Elevation [m]`,
      `Date/Time`
      ),
    sample_id = Event,
    everything()
  ) %>% 
  dplyr::slice(-49,) %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(counts = -group_cols()) %>% 
  dplyr::ungroup()

counts_latlong_933659 <- 
  dplyr::inner_join(
  lat_long_933659, 
  counts_933659,
  by = "sample_id"
  )

dat_933659 <- 
  tibble(
    dataset_id = dataset_id_933659,
    counts_latlong_933659,
    depositional_env = depositional_env_933659,
    publication = publication_933659
  ) %>% 
  dplyr::mutate(percentage = TRUE) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)


#10.1594/PANGAEA.717143
dt10 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.717143")

dataset_id_717143 <- "PANGAEA.717143"
publication_717143 <- dt10$citation
depositional_env_717143 <- "Lacustrine"
lat_long_717143 <- 
  tibble(
    lat = dt10$metadata[[1]]$events$LATITUDE,
    long = dt10$metadata[[1]]$events$LONGITUDE
  )
counts_717143 <- 
  dt10$data[[1]] %>% 
  dplyr::select(
    -c(
      "Sample label",
      "Depth [m]",
      "Age [ka BP]",
      "Slide count [#]",
      "Samp vol [cm**3]",
      "Lyc.added [#]",
      "Mark add [#]",
      "Mark found [#]",
      "Pollen conc [#/cm**3]"  
    )
  ) %>% 
  dplyr::mutate(sample_id = "TK_220") %>% 
  dplyr::select(-"Sec label") %>% 
  dplyr::group_by(sample_id) %>% 
  tidyr::nest(counts = -group_cols()) %>% 
  dplyr::mutate(percentage = FALSE) %>% 
  dplyr::ungroup()

dat_717143 <- 
  tibble(
    dataset_id = dataset_id_717143,
    publication = publication_717143,
    depositional_env = depositional_env_717143,
    lat_long_717143,
    counts_717143
    ) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)


#10.1594/PANGAEA.910726
dt11 <- 
  raw_dat %>% 
  dplyr::filter(doi == "10.1594/PANGAEA.910726")
dataset_id_910726 <- "PANGAEA.910726"
depositional_env_910726 <- "Moss"
publication_910726 <- dt11$citation

lat_long_910726 <-
  data.frame(
    strsplit(
      dt11$metadata[[1]]$events,
      '; ',
      fixed = TRUE
      )
    ) %>%
  rlang::set_names("V") %>%
  dplyr::mutate(
    V1 = purrr::map(
      V,
      function(.x) {
        dt <- 
          data.frame(
          strsplit(
            as.character(.x),
            ' * ',
            fixed = TRUE
            )
          ) %>%
          rlang::set_names("V1") 
        }
    )
  )
        
lat_long_proc <-  
  purrr::map(lat_long_910726$V1,
                     ~ .x %>%
                       dplyr::mutate(
                         V2 = 
                           ifelse(
                             stringr::str_detect(
                               V1, "river|LATITUDE|LONGITUDE"
                               ) == TRUE,
                             "Include",
                             "Exclude"
                             )
                         )
             ) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(V2 == "Include") %>%
  dplyr::select(-c(V2)) %>%
  dplyr::mutate(
    event = ifelse(
      stringr::str_detect(
        V1, "river"
        ) == TRUE, 
      "event", 
      "other"),
    lat = ifelse(
      stringr::str_detect(
        V1, "LATITUDE"
        ) == TRUE, 
      "lat", 
      "other"),
    long = ifelse(
      stringr::str_detect(
        V1, "LONGITUDE"
        ) == TRUE, 
      "long", 
      "other")
    )

  sample_id <- 
    lat_long_proc %>% 
    dplyr::select(V1,
                  event) %>% 
    dplyr::filter(event == "event") %>% 
    dplyr::select(sample_id = V1) %>% 
    dplyr::mutate(sample_label 
                  = c(
                    "B15",
                    "B16",
                    "B17", 
                    "B18-1",
                    "B18-2", 
                    "B18-3", 
                    "B29-1", 
                    "B29-2", 
                    "B29-3", 
                    "B30",
                    "B31"
                  )
                  ) %>% 
    dplyr::select(sample_label)
  
  lat_comp <- 
    lat_long_proc %>% 
    dplyr::select(V1,
                  lat) %>% 
    dplyr::filter(lat == "lat") %>% 
    dplyr::select(-lat) %>% 
    dplyr::rename(lat = V1)
    
  lat <- 
    data.frame(
      strsplit(lat_comp$lat, ': ',
             fixed = TRUE)) %>% 
    dplyr::slice(-1,) %>% 
    tidyr::gather() %>% 
  dplyr::select(lat = value) %>% 
  dplyr::mutate_at("lat", as.numeric)
      
  long_comp <- 
    lat_long_proc %>% 
    dplyr::select(V1,
                  long) %>% 
    dplyr::filter(long == "long") %>% 
    dplyr::select(-long) %>% 
    dplyr::rename(long = V1)
  
  long <- 
    data.frame(
      strsplit(long_comp$long, ': ',
               fixed = TRUE)) %>% 
    dplyr::slice(-1,) %>% 
    tidyr::gather() %>% 
    dplyr::select(long = value) %>% 
    dplyr::mutate_at("long", as.numeric)
  
  
lat_long <- 
  dplyr::bind_cols(
    sample_id,
    lat,
    long
    ) 

counts_910726 <- 
  dt11$data[[1]] %>% 
  dplyr::select(
    sample_id = Event,
    sample_label = "Sample label",
    everything()
    ) %>% 
  dplyr::group_by(
    sample_id,
    sample_label
    ) %>% 
  tidyr::nest(counts = -group_cols())

counts_lat_long_910726 <- 
  inner_join(
    lat_long,
    counts_910726,
    by = "sample_label") %>% 
dplyr::select(-sample_label) %>% 
  tibble::as_tibble() %>% 
  dplyr::ungroup()
    
dat_910726 <- 
  tibble(
    dataset_id = dataset_id_910726,
    depositional_env = depositional_env_910726,
    publication = publication_910726,
    counts_lat_long_910726
    ) %>% 
  dplyr::mutate(percentage = TRUE) %>% 
  dplyr::mutate_at("lat", as.numeric) %>% 
  dplyr::mutate_at("long", as.numeric)
  

surface_samples_pangaea <- 
  dplyr::bind_rows(
    dat_717143,
    dat_808953,
    dat_829751,
    dat_871523,
    dat_910726,
    dat_933659,
    dat1_849666
  )
  
counts_taxa <- 
  purrr::map(
    surface_samples_pangaea$counts,
    ~ colnames(.x) %>% 
      tibble::enframe(
        name = NULL,
        value = "taxa"
        ) %>% 
      dplyr::bind_rows()
    ) %>% 
  dplyr::bind_rows() %>% 
  dplyr::distinct()

clean_names <-
  counts_taxa %>%
  dplyr::mutate(
    abbrev = strsplit(
      as.character(taxa),
      " [#]",
      fixed = TRUE
      ),
    abbrev = 
      strsplit(
        as.character(abbrev),
        " [%]",
        fixed = TRUE
      )
    ) %>% 
  tidyr::unnest(abbrev) %>% 
  dplyr::distinct()

table_taxa <- 
  readr::read_csv(
    "Inputs/Tables/taxa_names_surface_samples_pangaea_150923.csv"
  ) %>% 
  dplyr::bind_rows(
    tibble(taxon.name = 
             c(
               "Rosaceae_nonoperculate",
               "Rosaceae_operculate"
               ),
           abbrev = 
             c(
               "c(\"Rosae\", \" (non-operculate)\")", 
               "c(\"Rosae\", \" (operculate)\")" 
               )
           )
    )

  
unique_taxa_dataset <- 
  clean_names %>% 
  dplyr::filter(
    !abbrev %in% table_taxa$abbrev
  ) # All can be discarded
unique_taxa_table <- 
  table_taxa %>% 
  dplyr::filter(
    !abbrev %in% clean_names$abbrev
  ) # doesn't matter


# Assign original names to taxa
raw_dat <- 
  surface_samples_pangaea %>% 
  dplyr::mutate(
    raw_counts = 
      purrr::map(
        .x = counts,
        .f = ~ {
          dat <- .x
          res <- 
            tidyr::gather(
              dat,
              key = "abbrev",
              value = "counts"
              ) %>% 
            dplyr::mutate(
              abbrev = strsplit(
                as.character(abbrev),
                " [#]",
                fixed = TRUE
              ),
              abbrev = 
                strsplit(
                  as.character(abbrev),
                  " [%]",
                  fixed = TRUE
                )
              ) %>% 
            tidyr::unnest(abbrev) %>% 
            dplyr::distinct() %>% 
            dplyr::inner_join(
              table_taxa,
              by = "abbrev"
            ) %>% 
            dplyr::select(
              taxon.name,
              counts
              ) %>% 
          tidyr::spread(
            taxon.name, counts
          )
            
          return(res)
        }
        )
    )

readr::write_rds(
  raw_dat,
  file = "Inputs/Data/surface_samples_cao/surface_samples_pangaea_220923.rds",
  compress = "gz"
)




