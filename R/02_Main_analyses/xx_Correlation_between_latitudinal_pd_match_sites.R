#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#   Correlation between latitudinal patterns of PD estimates ----
#
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source configuration ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
fossil_pollen <- 
  read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds")

surface_pollen <- 
  read_rds("Inputs/Data/surface_pollen_filtered_030124.rds")

# Make buffer envelop around each location of fossil pollen data to select nearby surface pollen data
get_modern_samples_buffer <- 
  function(location_x, 
           location_y, 
           crd_ref = CRS("+proj=longlat +datum=WGS84"), 
           projection = CRS("+proj=merc +datum=WGS84"), 
           width = 200*10^3) {
    
  # needs to be data.frame
  location_x = as.data.frame(location_x)
  location_y = as.data.frame(location_y)
  
  # make spatial point object
  pts_fossil_data <- 
    SpatialPoints(location_x %>%
                    dplyr::rename(Latitude = lat,
                                  Longitude = long) %>%
                    dplyr::select(Latitude, 
                                  Longitude),
                  proj4string = crd_ref)
  pts_surface_samples <- 
    SpatialPoints(location_y %>% 
                    dplyr::rename(Latitude = lat, 
                                  Longitude = long) %>% 
                    dplyr::select(Latitude, 
                                  Longitude), 
                  proj4string = crd_ref)
  
  pts_fossil_data_temp <- SpatialPointsDataFrame(pts_fossil_data, location_x)
  pts_surface_samples_temp <- SpatialPointsDataFrame(pts_surface_samples, location_y)
  
  # make buffers to find sites within certain distances from points
  pts_fossil_data_proj <- spTransform(pts_fossil_data_temp, projection)
  pts_surface_samples_proj <- spTransform(pts_surface_samples_temp, projection)
  site_buffer <- 
    gBuffer(pts_fossil_data_proj, 
            width = width, 
            byid = TRUE)
  #plot(site_Buffer)
  
  # overlay operation
  overlay_pts_fossil_data <- 
    over(pts_surface_samples_proj, 
         site_buffer) %>% 
    filter(!is.na(site))
  samples <- rownames(overlay_pts_fossil_data)
  samples
  }



# extract relevant data
data_input_fossil <- 
  fossil_pollen  %>%
  dplyr::select(sample_id, 
                lat, 
                long, 
                climate_zone_revised) %>% 
  tidyr::nest(sites = c(sample_id, 
                        lat, 
                        long)
              )
data_input_surface <- 
  surface_pollen %>% 
  dplyr::select(sample_id, 
                lat, 
                long,
                climate_zone_revised) %>% 
  tidyr::nest(sites = c(sample_id, 
                        lat, 
                        long)
              ) %>% 
  dplyr::filter(!is.na(climate_zone_revised))

# Identify surface samples within the given distance; Input data are geographical locations of surface pollen samples and fossil sites
buffer_sites <- purrr::map2(data_input_fossil$sites, data_input_fossil$climate_zone_revised, .f = get_modern_samples_buffer)
names.buffer <- map_chr(data.input$sites, "site")
buffer_sites_t <- tibble(site = names.buffer,
                         buffer_sites = buffer_sites)
# attach buffer site to the tibble
LIG.climate <- LIG.climate %>% left_join(buffer_sites_t, by = "site")