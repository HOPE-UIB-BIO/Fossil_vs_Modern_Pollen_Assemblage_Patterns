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
# Latitudinal trend of phylogenetic dispersion and compositional turnover ----
#
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source configuration ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
phylo_div_full <- 
  readr::read_rds("Inputs/Data/phylo_div_full_090123.rds")

turnover_combined <- 
  readr::read_rds("Inputs/Data/turnover_combined_050224.rds") %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest(.key = "turnover") %>% 
  dplyr::ungroup()

#--------------------------------------------------------#
# 3. Fit the GAM models ----
#--------------------------------------------------------#

#---------------------------------#
# A. Phylogenetic dispersion (sesMPD, sesMNTD) ----
#---------------------------------#

# Surface samples
surface_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "surface_samples") %>% 
  dplyr::select(-dataset_id) %>% 
  dplyr::group_by(lat) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(dataset_id = paste("dataset_", row.names(.), sep = "")) %>% 
  tidyr::unnest(data) %>% 
  tidyr::gather(
    c(ses_mpd, ses_mntd),
    key = "var",
    value = "estimate"
    ) %>% 
  dplyr::select(
    dataset_id, 
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised,
    var, 
    estimate
    ) %>% 
  dplyr::mutate_at("dataset_id", as.factor) %>% 
  dplyr::group_by(var) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup()

gam_mod_surface_samples <-
  surface_samples %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = data, 
        .f = ~ {
        data <- .x
        set.seed(2330)
        mod <-
          mgcv::gam(
            estimate ~
              s(lat, 
                k = 10, 
                bs = 'tp'
                ) + 
              s(dataset_id,
                bs= 're'),
            data = data,
            method = "REML",
            family = "gaussian",
            control = gam.control(
              trace = TRUE, 
              maxit = 200
              )
          )
        }
      ),
    
    predicted_gam =
      purrr::map2(
        .x = data, 
        .y = gam_model, 
        .f = ~ {
        data <- .x
        new_data_gam <-
          with(
            data,
            base::expand.grid(
              lat = seq(
                min(25.50), max(65), by = 0.1
                ),
              dataset_id = dataset_id[1]
              )
            )
         
        crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE)
        set.seed(2330)
        predicted_mod <-
          new_data_gam %>%
          dplyr::bind_cols(
            predict(
              .y,
              newdata = new_data_gam,
              type = "response",
              se.fit = TRUE,
              exclude = "dataset_id"
              )
            ) %>%
          
          dplyr::mutate(
            var = fit,
            lwr = fit - (crit * se.fit),
            upr = fit + (crit * se.fit)
            ) %>%
          dplyr::select(
            !dplyr::any_of(
              c("fit", "se.fit")
              )
            )
        return(predicted_mod)
        }
      )
    )


# Top 500-yr fossil pollen samples 
top_500_fossil_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "top_500_yr") %>% 
  tidyr::gather(
    c(ses_mpd, ses_mntd),
    key = "var",
    value = "estimate"
    ) %>% 
  dplyr::select(
    dataset_id, 
    sample_id, 
    lat, 
    long, 
    ecozone_koppen_5,
    ecozone_koppen_15,
    ecozone_koppen_30,
    climate_zone_revised,
    var, 
    estimate
    ) %>% 
  dplyr::mutate_at("dataset_id", as.factor) %>% 
  dplyr::mutate_at("sample_id", as.factor) %>% 
  dplyr::group_by(var) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup()


gam_mod_top_500_fossil_samples <-
  top_500_fossil_samples %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = data, 
        .f = ~ {
          data <- .x
          set.seed(2330)
          mod <-
            mgcv::gam(
              estimate ~
                s(lat, 
                  k = 10, 
                  bs = 'tp'
                ) + 
               s(dataset_id,
                  bs= 're'),
              data = data,
              method = "REML",
              family = "gaussian",
              control = gam.control(
                trace = TRUE, 
                maxit = 200
              )
            )
        }
      ),
    
    predicted_gam =
      purrr::map2(
        .x = data, 
        .y = gam_model, 
        .f = ~ {
          data <- .x
          new_data_gam <-
            with(
              data,
              base::expand.grid(
                lat = seq(
                  min(25.50), max(65), by = 0.1
                  ),
                dataset_id = dataset_id[1]
              )
            )
          
          crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE)
          set.seed(2330)
          predicted_mod <-
            new_data_gam %>%
            dplyr::bind_cols(
              predict(
                .y,
                newdata = new_data_gam,
                type = "response",
                se.fit = TRUE,
                exclude = "dataset_id"
              )
            ) %>%
            
            dplyr::mutate(
              var = fit,
              lwr = fit - (crit * se.fit),
              upr = fit + (crit * se.fit)
            ) %>%
            dplyr::select(
              !dplyr::any_of(
                c("fit", "se.fit")
              )
            )
          
          return(predicted_mod)
        }
      )
  )

#--------------------------------------------------------#
# Save outputs ----
#--------------------------------------------------------#

readr::write_rds(gam_mod_surface_samples,
                 file = "Outputs/Data/gam_mod_surface_samples_100124.rds",
                 compress = "gz")

readr::write_rds(gam_mod_top_500_fossil_samples,
          file = "Outputs/Data/gam_mod_top_500_fossil_samples_100124.rds",
          compress = "gz")

#---------------------------------#
# B. Compositional turnover ----
#---------------------------------#
gam_mod_turnover <-
  turnover_combined %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = turnover, 
        .f = ~ {
          data <- 
            .x %>% 
            dplyr::mutate_at("dataset_id", as.factor) %>% 
            dplyr::mutate_at("sample_id", as.factor)
          
          set.seed(2330)
          mod <-
            mgcv::gam(
              axis_1 ~
                s(lat, 
                  k = 5, 
                  bs = 'tp'
                ) + 
                s(dataset_id,
                  bs= 're'),
              data = data,
              method = "REML",
              family = "gaussian",
              control = gam.control(
                trace = TRUE, 
                maxit = 200
              )
            )
        }
      ),
    
    predicted_gam =
      purrr::map2(
        .x = turnover, 
        .y = gam_model, 
        .f = ~ {
          data <- .x
          new_data_gam <-
            with(
              data,
              base::expand.grid(
                lat = seq(
                  min(25.50), max(65), by = 0.1
                ),
                dataset_id = dataset_id[1]
              )
            )
          
          crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE)
          set.seed(2330)
          predicted_mod <-
            new_data_gam %>%
            dplyr::bind_cols(
              predict(
                .y,
                newdata = new_data_gam,
                type = "response",
                se.fit = TRUE,
                exclude = "dataset_id"
              )
            ) %>%
            
            dplyr::mutate(
              var = fit,
              lwr = fit - (crit * se.fit),
              upr = fit + (crit * se.fit)
            ) %>%
            dplyr::select(
              !dplyr::any_of(
                c("fit", "se.fit")
              )
            )
          
          return(predicted_mod)
        }
      )
  )

#--------------------------------------------------------#
# Save output ----
#--------------------------------------------------------#
readr::write_rds(gam_mod_turnover,
          file = "Outputs/Data/gam_mod_turnover_050224.rds",
          compress = "gz")
