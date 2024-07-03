#----------------------------------------------------------#
# Fossil pollen data can predict robust spatial patterns of biodiversity 
#                        in the past
#
#                       Project setup
#                 
#                         K. Bhatta 
#
#                           2024
#                         
#----------------------------------------------------------#

# Script to prepare all necessary components of environment to run the Project.
#  Needs to be run only once

#----------------------------------------------------------#
# Step 1: Install 'renv' package -----
#----------------------------------------------------------#

utils::install.packages("renv")


#----------------------------------------------------------#
# Step 2: Deactivate 'renv' package -----
#----------------------------------------------------------#

# deactivate to make sure that packages are updated on the machine
renv::deactivate()


#----------------------------------------------------------#
# Step 3: Create a list of packages
#----------------------------------------------------------#

package_list <- 
  c(
    "ape",
    "assertthat",
    "cowplot",
    "devtools",
    "ggpubr",
    "ggmap",
    "ggpmisc",
    "ggspatial",
    "ggtext",
    "glmm",
    "glmmTMB",
    "gratia",
    "grid",
    "gridExtra",
    "here",
    "httpgd",
    "janitor",
    "jsonlite",
    "languageserver",
    "leaflet",
    "lme4", 
    "lmodel2",
    "mgcv",
    "performance",
    "picante",
    "raster",
    "RColorBrewer",
    "remotes",
    "renv",       
    "roxygen2",   
    "tidyverse",  
    "usethis",
    "vegan",
    "viridis"
  )

# define helper function
install_packages <-
  function(pkgs_list) {

    # install all packages in the lst from CRAN
    sapply(pkgs_list, utils::install.packages, character.only = TRUE)

     #install 'RFossilpol' from GitHub
     remotes::install_github(
      "HOPE-UIB-BIO/R-Fossilpol-package"
      )
     
     # Install 'REcopol' from FitHub
     remotes::install_github(
       "HOPE-UIB-BIO/R-Ecopol-package")
  }

#----------------------------------------------------------#
# Step 4: Install packages to the machine
#----------------------------------------------------------#

install_packages(package_list)


#----------------------------------------------------------#
# Step 5: Activate 'renv' project
#----------------------------------------------------------#

renv::activate()


#----------------------------------------------------------#
# Step 6: Install packages to the project
#----------------------------------------------------------#

install_packages(package_list)


#----------------------------------------------------------#
# Step 7: Synchronize package versions with the project 
#----------------------------------------------------------#

library(here)

# if there is no lock file present make a new snapshot
if (
  isFALSE("library_list.lock" %in% list.files(here::here("renv")))
  ) {
  renv::snapshot(lockfile = here::here("renv/library_list.lock"))
} else {
  renv::restore(lockfile = here::here("renv/library_list.lock"))
}


#----------------------------------------------------------#
# Step 8: Run the project 
#----------------------------------------------------------#
