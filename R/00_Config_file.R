#----------------------------------------------------------#
# Fossil pollen data can predict robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#                     
#                     Configuration file 
#----------------------------------------------------------#

#----------------------------------------------------------#
# 1. Load packages -----
#----------------------------------------------------------#

# define packages

package_list <- 
  c(
    "ape",
    "assertthat",
    "corrplot", 
    "cowplot",
    "devtools",
    "dplyr",
    "ggplot2",
    "ggpubr",
    "ggmap",
    "ggpmisc",
    "ggspatial",
    "ggtext",
    "gratia",
    "grid",
    "gridExtra",
    "here",   
    "leaflet",
    "lmodel2",
    "maps",
    "mgcv", 
    "pangaear",
    "performance",
    "picante",
    "raster",
    "RColorBrewer",
    "REcopol",
    "renv", 
    "RFossilpol",
    "rpart",
    "roxygen2", 
    "scales",
    "sf",
    "sjmisc",
    "sjstats", 
    "tidyverse",  
    "usethis",
    "vegan",
    "viridis",
    "xaringan"
    )

# load all packages
sapply(package_list, library, character.only = TRUE)

# Note: new required packages should be installed here!

if(!exists("update_repo_packages")){
  update_repo_packages <- TRUE
}

if(update_repo_packages == TRUE){
  
  if (!exists("already_synch")){
    already_synch <- FALSE
  }
  
  if(already_synch == FALSE){
    library(here)
    # synchronise the package versions
    renv::restore(lockfile = here::here( "renv/library_list.lock"))
    already_synch <- TRUE
    
    # save snapshot of package versions
    renv::snapshot(lockfile =  "renv/library_list.lock")  # do only for update
  }
}

#----------------------------------------------------------#
# 2. Define space -----
#----------------------------------------------------------#
# project directory is set up by 'here' package, Adjust if needed 
current_dir <- here::here()

#----------------------------------------------------------#
# 3. Source functions -----
#----------------------------------------------------------#
invisible(lapply(
    list.files(
      path = "R/functions/",
      pattern = "*.R",
      recursive = TRUE,
      full.names = TRUE
    ),
    source)
  )


#----------------------------------------------------------#
# 4. criteria for age limits of period of interest
#----------------------------------------------------------#
continents <- 
  c("North America", "Europe", "Asia", "Latin America", "Africa", "Oceania")

continental_age_limits <- 
  tibble(
    continent = continents,
    young_age = c(rep(1e3, 2), rep(1e3, 4)),
    old_age = c(rep(8e3, 2), rep(6e3, 4)),
    end_of_interest_period = c(rep(12e3, 6)))

#----------------------------------------------------------#
# 5. criteria to filter out levels (sample depths)
#----------------------------------------------------------#
min_n_levels <-  5  # at least 5 levels within period of interest 
# (young_age to end_of_interest_period)

min_n_grains <- 25 # each level at least 25 counted individuals 

target_n_grains <- 150 # ideal number of counts

percentage_samples <- 50 # threshold of number of samples with ideal counts

maximum_age_extrapolation <- 3000 # how much age can be extrapolated beyond
# the oldest chronology control point

#----------------------------------------------------------#
# 6. Colour scheme ----
#----------------------------------------------------------#
color_pd_curve <- "#0072B2"
color_richness_curve <- "#D55E00"
color_phylo_age_curve <- "#666666" 
color_common <- "#000000" 
  
color_high_lat <- "#0099FF"
color_low_lat <- "#FF0000"

color_high_age <- "#009E73"
color_low_age <- "#E69F00"

# Assign unique colour to each climate zone
my_palette <-
  c("#FFCC99",
    "#993300",
    "#9999CC",
    "#3399FF",
    "#999999",
    "#00CCCC",
    "#99CC00",
    "#006600",
    "#996600"
  ) %>% 
  rlang::set_names(
    nm = c(
      "Arid",
      "Cold-seasonally dry",
      "Cold without dry season",
      "Polar",
      "Polar frost",
      "Temperate",
      "Tropical monsoon",
      "Tropical rainforest",
      "Tropical savannah"
    )
  )

