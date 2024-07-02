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
#       Latitudinal analysis of phylogenetic dispersion
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
  read_rds("Inputs/Data/phylo_div_full_090123.rds")

turnover_combined <- 
  read_rds("Inputs/Data/turnover_combined_050224.rds") %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest(.key = "turnover") %>% 
  ungroup()

#--------------------------------------------------------#
# 3. Fit the GAM models ----
#--------------------------------------------------------#

#---------------------------------#
#A. Phylogenetic dispersion (sesMPD, sesMNTD) ----
#---------------------------------#
#Surface samples
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
write_rds(gam_mod_surface_samples,
          file = "Outputs/Data/gam_mod_surface_samples_100124.rds",
          compress = "gz")
#-----------------------------------------------#

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
write_rds(gam_mod_top_500_fossil_samples,
          file = "Outputs/Data/gam_mod_top_500_fossil_samples_100124.rds",
          compress = "gz")

#-----------------------------------------------#
# Top 1000-yr fossil pollen samples 
top_1000_fossil_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "top_1000_yr") %>% 
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


gam_mod_top_1000_fossil_samples <-
  top_1000_fossil_samples %>%
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
write_rds(gam_mod_top_1000_fossil_samples,
          file = "Outputs/Data/gam_mod_top_1000_fossil_samples_100124.rds",
          compress = "gz")

# Modern samples
modern_samples <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "modern_vegetaton") %>% 
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
  dplyr::select(dataset_id, sample_id, lat, var, estimate) %>% 
  dplyr::mutate_at("dataset_id", as.factor) %>% 
  dplyr::group_by(var) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup()

gam_mod_modern_samples <-
  modern_samples %>%
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
                  ),
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
                  min(lat), max(lat), by = 0.25
                )
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
                se.fit = TRUE#,
                #exclude = "dataset_id"
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
write_rds(gam_mod_modern_samples,
          file = "Outputs/Data/gam_mod_modern_samples_180124.rds",
          compress = "gz")


#---------------------------------#
# B. Compositional turnover
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
write_rds(gam_mod_turnover,
          file = "Outputs/Data/gam_mod_turnover_050224.rds",
          compress = "gz")

#---------------------------------#
# C. Compositional variation (DCA)
#---------------------------------#
gam_mod_dca <-
  dca_combined %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = dca, 
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
        .x = dca, 
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
write_rds(gam_mod_dca,
          file = "Outputs/Data/gam_mod_dca_060524.rds",
          compress = "gz")

#--------------------------------------------------------#
# 4. Plot GAM models ----
#--------------------------------------------------------#
# Surface samples
gam_mod_surface_samples_pd <- 
  read_rds("Outputs/Data/gam_mod_surface_samples_100124.rds")
gam_mod_turnover <- 
  read_rds("Outputs/Data/gam_mod_turnover_050224.rds")
gam_mod_dca <- 
  read_rds("Outputs/Data/gam_mod_dca_060524.rds")

# sesMPD and sesMNTD
plot_gam_mod_surface_samples_pd <-
  purrr::pmap(
    .l = list(
      gam_mod_surface_samples_pd$data, # ..1
      gam_mod_surface_samples_pd$predicted_gam, # ..2
      gam_mod_surface_samples_pd$var # ..3
    ),
    .f = ~ {
      ggplot2::ggplot(
        aes(
          x = lat, 
          y = estimate
          ),
        data = ..2
        ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = estimate,
            colour = climate_zone_revised
            ), 
          data = ..1,
          size = 1,
          alpha = 0.5
          ) + 
        ggplot2::scale_colour_manual(values = my_palette)+
        ggplot2::scale_fill_manual(values = my_palette) + 
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
            ),
          linewidth = 1.5,
          alpha = 1,
          data = ..2,
          colour = "#D55E00"
          ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
              'Latitude ', (degree ~ N)
            )
          ),
          y = paste(..3),
          colour = "Climate zone",
          fill = "Climate zone"
          ) + 
        
        theme(legend.position = "none",
              axis.title = element_text(
                size = 15,
                color = color_common
              ),
              axis.text = element_text(
                size = 13,
                color = color_common
              )
        )
    }
  )


# Compositional turnover
plot_gam_mod_turnover <-
  purrr::map2(
    .x = gam_mod_turnover$turnover, 
    .y = gam_mod_turnover$predicted_gam,
    .f = ~ {
      
      ggplot2::ggplot(
        aes(
          x = lat, 
          y = var
        ),
        data = .y
      ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = axis_1,
            colour = climate_zone_revised
          ), 
          data = .x,
          size = 0.9,
          alpha = 0.5
        ) + 
        ggplot2::scale_colour_manual(values = my_palette)+
        ggplot2::scale_fill_manual(values = my_palette) + 
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
          ),
          linewidth = 1.5,
          alpha = 1,
          data = .y,
          colour = "#D55E00"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
              'Latitude ', (degree ~ N)
            )
          ),
          y = paste("DCCA axis 1"),
          colour = "Climate zone",
          fill = "Climate zone"
        ) + 
        
        theme(legend.position = "none",
              axis.title = element_text(
                size = 15,
                color = color_common
              ),
              axis.text = element_text(
                size = 13,
                color = color_common
              )
        )
    }
  )
# Make plot separately to make large data points for fossil datasets, which have fewer and non-overlapping data points
plot_gam_mod_turnover_top_500 <-
  purrr::map2(
    .x = gam_mod_turnover[2,]$turnover, 
    .y = gam_mod_turnover[2,]$predicted_gam,
    .f = ~ {
      
      ggplot2::ggplot(
        aes(
          x = lat, 
          y = var
        ),
        data = .y
      ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = axis_1,
            colour = climate_zone_revised
          ), 
          data = .x,
          size = 2.5,
          alpha = 0.5
        ) + 
        ggplot2::scale_colour_manual(values = my_palette)+
        ggplot2::scale_fill_manual(values = my_palette) + 
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
          ),
          linewidth = 1.5,
          alpha = 1,
          data = .y,
          colour = "#D55E00"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
              'Latitude ', (degree ~ N)
            )
          ),
          y = paste("DCCA axis 1"),
          colour = "Climate zone",
          fill = "Climate zone"
        ) + 
        theme(legend.position = "none",
              axis.title = element_text(
                size = 15,
                color = color_common
              ),
              axis.text = element_text(
                size = 13,
                color = color_common
              )
        )
    }
  )


# DCA
plot_gam_mod_dca <-
  purrr::map2(
    .x = gam_mod_dca$dca, 
    .y = gam_mod_dca$predicted_gam,
    .f = ~ {
      
      ggplot2::ggplot(
        aes(
          x = lat, 
          y = var
        ),
        data = .y
      ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = axis_1,
            colour = climate_zone_revised
          ), 
          data = .x,
          size = 0.9,
          alpha = 0.5
        ) + 
        ggplot2::scale_colour_manual(values = my_palette)+
        ggplot2::scale_fill_manual(values = my_palette) + 
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
          ),
          linewidth = 1.5,
          alpha = 1,
          data = .y,
          colour = "#D55E00"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
              'Latitude ', (degree ~ N)
            )
          ),
          y = paste("DCA axis 1"),
          colour = "Climate zone",
          fill = "Climate zone"
        ) + 
        
        theme(legend.position = "none",
              axis.title = element_text(
                size = 15,
                color = color_common
              ),
              axis.text = element_text(
                size = 13,
                color = color_common
              )
        )
    }
  )

# Make plot separately to make large data points for fossil datasets, which have fewer and non-overlapping data points
plot_gam_mod_dca_top_500 <-
  purrr::map2(
    .x = gam_mod_dca[2,]$dca, 
    .y = gam_mod_dca[2,]$predicted_gam,
    .f = ~ {
      
      ggplot2::ggplot(
        aes(
          x = lat, 
          y = var
        ),
        data = .y
      ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = axis_1,
            colour = climate_zone_revised
          ), 
          data = .x,
          size = 2.5,
          alpha = 0.5
        ) + 
        ggplot2::scale_colour_manual(values = my_palette)+
        ggplot2::scale_fill_manual(values = my_palette) + 
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
          ),
          linewidth = 1.5,
          alpha = 1,
          data = .y,
          colour = "#D55E00"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
              'Latitude ', (degree ~ N)
            )
          ),
          y = paste("DCA axis 1"),
          colour = "Climate zone",
          fill = "Climate zone"
        ) + 
        theme(legend.position = "none",
              axis.title = element_text(
                size = 15,
                color = color_common
              ),
              axis.text = element_text(
                size = 13,
                color = color_common
              )
        )
    }
  )


fig_1 <- plot_gam_mod_surface_samples_pd[[1]] + 
  guides(fill = guide_legend(nrow = 3,
                             title.position = "top"),
         colour = guide_legend(nrow = 3,
                               title.position = "top")
         ) +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
        )
fig_2 <- plot_gam_mod_surface_samples_pd[[2]] 
fig_3 <- plot_gam_mod_turnover[[1]] 
fig_dca <- plot_gam_mod_dca[[1]]

composite_surface_samples <- 
  ggpubr::ggarrange(
    fig_1, fig_2, fig_3,
    ncol = 1,
    nrow = 3,
    labels = c("(a)", "(c)", "(e)"),
    font.label = list(size = 16, color = color_common),
    label.x = 0.15,
    common.legend = TRUE,
    legend = "bottom"
    ) 
final_fig_surface_samples <- 
  annotate_figure(
    composite_surface_samples,
    top = text_grob(
      "Surface pollen data",
      color = color_common, 
      size = 18
    )
  )

composite_surface_samples_dca <- 
  ggpubr::ggarrange(
    fig_1, fig_2, fig_dca,
    ncol = 1,
    nrow = 3,
    labels = c("(a)", "(c)", "(e)"),
    font.label = list(size = 16, color = color_common),
    label.x = 0.15,
    common.legend = TRUE,
    legend = "bottom"
  ) 
final_fig_surface_samples_dca <- 
  annotate_figure(
    composite_surface_samples_dca,
    top = text_grob(
      "Surface pollen data",
      color = color_common, 
      size = 18
    )
  )


# Top 500-yr fossil pollen samples
gam_mod_top_500_fossil_samples_pd <- 
  read_rds("Outputs/Data/gam_mod_top_500_fossil_samples_100124.rds")

plot_gam_500_pd <-
  purrr::pmap(
    .l = list(
      gam_mod_top_500_fossil_samples_pd$data, # ..1
      gam_mod_top_500_fossil_samples_pd$predicted_gam, # ..2
      gam_mod_top_500_fossil_samples_pd$var # ..3
    ),
    .f = ~ {
      
      ggplot2::ggplot(
        aes(
          x = lat, 
          y = estimate
          ),
        data = ..2
        ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = estimate,
            colour = climate_zone_revised
            ), 
          data = ..1,
          size = 2.5,
          alpha = 0.5
          ) + 
        ggplot2::scale_colour_manual(values = my_palette) +
        ggplot2::scale_fill_manual(values = my_palette) + 
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
            ),
          linewidth = 1.5,
          alpha = 1,
          data = ..2,
          colour = "#D55E00"
          ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
              'Latitude ', (degree ~ N)
            )
          ),
          y = paste(..3),
          colour = "Climate zone",
          fill = "Climate zone"
          ) + 
        
        theme(legend.position = "none",
              axis.title = element_text(
                size = 15,
                color = color_common
              ),
              axis.text = element_text(
                size = 13,
                color = color_common
              )
        )
    }
  )
fig_4 <- plot_gam_500_pd[[1]] +
  guides(fill = guide_legend(nrow = 3,
                             title.position = "top"),
         colour = guide_legend(nrow = 3,
                               title.position = "top")
         ) +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
        )
fig_5 <- plot_gam_500_pd[[2]]
fig_6 <- plot_gam_mod_turnover_top_500[[1]] 
fig_6_dca <- plot_gam_mod_dca_top_500[[1]] 


composite_top_500 <- 
  ggpubr::ggarrange(
    fig_4, fig_5, fig_6,
    ncol = 1,
    nrow = 3,
    common.legend = TRUE,
    legend = "bottom",
    labels = c("(b)", "(d)", "(f)"),
    font.label = list(size = 16, color = color_common),
    label.x = 0.15
    ) 

final_fig_top_500 <- 
  annotate_figure(
    composite_top_500,
    top = text_grob(
      "Fossil pollen data",
      color = color_common, 
      size = 18
    )
  )

fig_final_composite <- 
  ggpubr::ggarrange(
    final_fig_surface_samples,
    final_fig_top_500,
    ncol = 2,
    nrow = 1
    )


composite_top_500_dca <- 
  ggpubr::ggarrange(
    fig_4, fig_5, fig_6_dca,
    ncol = 1,
    nrow = 3,
    common.legend = TRUE,
    legend = "bottom",
    labels = c("(b)", "(d)", "(f)"),
    font.label = list(size = 16, color = color_common),
    label.x = 0.15
  ) 

final_fig_top_500 <- 
  annotate_figure(
    composite_top_500,
    top = text_grob(
      "Fossil pollen data",
      color = color_common, 
      size = 18
    )
  )

final_fig_top_500_dca <- 
  annotate_figure(
    composite_top_500_dca,
    top = text_grob(
      "Fossil pollen data",
      color = color_common, 
      size = 18
    )
  )

fig_final_composite <- 
  ggpubr::ggarrange(
    final_fig_surface_samples,
    final_fig_top_500,
    ncol = 2,
    nrow = 1
  )

fig_final_composite_dca <- 
  ggpubr::ggarrange(
    final_fig_surface_samples_dca,
    final_fig_top_500_dca,
    ncol = 2,
    nrow = 1
  )

  
ggplot2::ggsave(
  fig_final_composite,
  filename = paste(
    "Outputs/Figure/",
    "Phylogenetic_dispersion_turnover_latitudinal_composite_210524.tiff",
    sep = ""
    ),
  height = 25,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )

ggplot2::ggsave(
  fig_final_composite_dca,
  filename = paste(
    "Outputs/Figure/",
    "Phylogenetic_dispersion_dca_latitudinal_composite_210524.tiff",
    sep = ""
  ),
  height = 25,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )


#--------------------------------------------------------#
# Maps of PAPs ----
#--------------------------------------------------------#
# Base map
base_map <-
  phylo_div_full %>% 
  dplyr::filter(data_type == "surface_samples") %>%
  ggplot(aes(x = long, y = lat)) +
  coord_fixed(ylim = c(20.00, 80.00),
              xlim = c(65.00, 135.00)
              ) + 
  labs(x = expression(paste('Longitude ', (degree ~ E))),
       y = expression(paste('Latitude ', (degree ~ N))),
       fill = "Climate-zones"
       ) +
  theme_classic() +
  borders(colour = color_common,
          linewidth = 0.2) +
  theme(
    axis.title = element_text(color = color_common, size = 16),
    axis.text = element_text(colour = color_common, size = 12),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) 

# Surface samples
surface_samples_pd <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "surface_samples")
surface_samples_turnover <- 
  turnover_combined[1,]$turnover[[1]] 
surface_samples_dca <- 
  dca_combined[1,]$dca[[1]]

# Fossil samples
top_500_fossil_samples_pd <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "top_500_yr")
top_500_fossil_turnover <-  
  turnover_combined[2,]$turnover[[1]] 
top_500_fossil_dca <-  
  dca_combined[2,]$dca[[1]] 

# MAP sesMPD
surface_samples_mpd <- 
  base_map +
  geom_point(
    aes(size = ses_mpd,
        colour = climate_zone_revised
        ),
    data = surface_samples_pd) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(-4,5)) + 
  labs(
    size = "sesMPD",
    colour = "Climate zone",
    fill = "Climate zone"
    ) + 
  theme(
    legend.box = "horizontal"
    ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
    )

#ggplot2::ggsave(
#  surface_samples_mpd,
#  filename = paste(
#    "Outputs/Figure/",
#    "Map_surface_samples_MPD_180224.tiff",
#    sep = ""
#    ),
#  height = 20,
#  width = 20,
#  units = "cm",
#  dpi = 400,
#  compression = "lzw"
#  )

top_500_fossil_samples_mpd <- 
  base_map +
  geom_point(
    aes(size = ses_mpd,
        colour = climate_zone_revised
        ),
    data = top_500_fossil_samples_pd) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(-4,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "sesMPD") + 
  theme(
    legend.box = "horizontal"
    ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
  )

final_map_mpd <- 
  ggpubr::ggarrange(
    surface_samples_mpd,
    top_500_fossil_samples_mpd,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    label.x = 0.1
    )
#ggplot2::ggsave(
#  top_500_fossil_samples_mpd,
#  filename = paste(
#    "Outputs/Figure/",
#    "Map_top_500yr_samples_MPD_180224.tiff",
#    sep = ""
#  ),
#  height = 20,
#  width = 20,
#  units = "cm",
#  dpi = 400,
#  compression = "lzw"
#)
ggplot2::ggsave(
  final_map_mpd,
  filename = paste(
    "Outputs/Figure/",
    "Map_sesMPD_210524.tiff",
    sep = ""
  ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )

# Map sesMNTD
surface_samples_mntd <- 
  base_map +
  geom_point(
    aes(size = ses_mntd,
        colour = climate_zone_revised
        ),
    data = surface_samples_pd) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(-4,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "sesMNTD") + 
  theme(
         legend.box = "horizontal"
       ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
    )

#ggplot2::ggsave(
#  surface_samples_mntd,
#  filename = paste(
#    "Outputs/Figure/",
#    "Map_surface_samples_MNTD_180224.tiff",
#    sep = ""
#    ),
#  height = 20,
#  width = 20,
#  units = "cm",
#  dpi = 400,
#  compression = "lzw"
#  )

top_500_fossil_samples_mntd <- 
  base_map +
  geom_point(
    aes(size = ses_mntd,
        colour = climate_zone_revised
        ),
    data = top_500_fossil_samples_pd) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(-4,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "sesMNTD") + 
  theme(
    legend.box = "horizontal"
    ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
    )

final_map_mntd <- 
  ggpubr::ggarrange(
    surface_samples_mntd,
    top_500_fossil_samples_mntd,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    label.x = 0.1
  )

#ggplot2::ggsave(
#  top_500_fossil_samples_mntd,
#  filename = paste(
#    "Outputs/Figure/",
#    "Map_top_500yr_samples_MNTD_180224.tiff",
#    sep = ""
#  ),
#  height = 20,
#  width = 20,
#  units = "cm",
#  dpi = 400,
#  compression = "lzw"
#)
ggplot2::ggsave(
  final_map_mntd,
  filename = paste(
    "Outputs/Figure/",
    "Map_sesMNTD_210524.tiff",
    sep = ""
  ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
)

# Map compositional turnover
surface_sample_turnover <- 
  base_map +
  geom_point(
    aes(size = axis_1,
        colour = climate_zone_revised
        ),
    data = surface_samples_turnover) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(-4,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "DCCA axis 1") +
  theme(
    legend.box = "horizontal"
    ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
    )

#ggplot2::ggsave(
#  surface_samples_turnover,
#  filename = paste(
#    "Outputs/Figure/",
#    "Map_surface_samples_turnover_180224.tiff",
#    sep = ""
#    ),
#  height = 20,
#  width = 20,
#  units = "cm",
#  dpi = 400,
#  compression = "lzw"
#  )


top_500_fossil_sample_turnover <- 
  base_map +
  geom_point(
    aes(size = axis_1,
        colour = climate_zone_revised
    ),
    data = top_500_fossil_turnover) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(0,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "DCCA axis 1") +
  theme(
  legend.box = "horizontal"
  ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
  )

final_map_turnover <- 
  ggpubr::ggarrange(
    surface_sample_turnover,
    top_500_fossil_sample_turnover,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    label.x = 0.1
    )

#ggplot2::ggsave(
#  top_500_fossil_samples_turnover,
#  filename = paste(
#    "Outputs/Figure/",
#    "Map_top_500yr_samples_turnover_180224.tiff",
#    sep = ""
#  ),
#  height = 20,
#  width = 20,
#  units = "cm",
#  dpi = 400,
#  compression = "lzw"
#  )

ggplot2::ggsave(
  final_map_turnover,
  filename = paste(
    "Outputs/Figure/",
    "Map_turnover_210524.tiff",
    sep = ""
  ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
)


# Map compositional variation DCA
surface_sample_dca <- 
  base_map +
  geom_point(
    aes(size = axis_1,
        colour = climate_zone_revised
    ),
    data = surface_samples_dca) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(-4,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "DCA axis 1") +
  theme(
    legend.box = "horizontal"
    ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
    )

top_500_fossil_sample_dca <- 
  base_map +
  geom_point(
    aes(size = axis_1,
        colour = climate_zone_revised
        ),
    data = top_500_fossil_dca) + 
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_size(range = c(0,5)) + 
  labs(colour = "Climate zone",
       fill = "Climate zone",
       size = "DCA axis 1") +
  theme(
    legend.box = "horizontal"
    ) +
  guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
    )

final_map_dca <- 
  ggpubr::ggarrange(
    surface_sample_dca,
    top_500_fossil_sample_dca,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    label.x = 0.1
    )

ggplot2::ggsave(
  final_map_dca,
  filename = paste(
    "Outputs/Figure/",
    "Map_dca_210524.tiff",
    sep = ""
  ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
)






