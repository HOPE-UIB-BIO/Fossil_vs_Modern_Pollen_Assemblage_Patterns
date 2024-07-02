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

# We can not use PD estimates as such for direct correlations because of the 
# difference in the size of data for each data type. Therefore, we fitted GAM model 
# of latitudinal pattern with each data type and predicted the model in a new 
# data frame with similar size for each model

# Predicted GAMs
gam_mod_surface_samples <- 
  read_rds("Outputs/Data/gam_mod_surface_samples_100124.rds") %>% 
  dplyr::select(var, surface_samples = predicted_gam) 

gam_mod_top_500_fossil_samples <- 
  read_rds("Outputs/Data/gam_mod_top_500_fossil_samples_100124.rds") %>% 
  dplyr::select(var, top_500_yr_fossil_samples = predicted_gam)

gam_mod_top_1000_fossil_samples <- 
  read_rds("Outputs/Data/gam_mod_top_1000_fossil_samples_100124.rds") %>% 
  dplyr::select(var, top_1000_yr_fossil_samples = predicted_gam)

cor_data <-
  gam_mod_surface_samples %>%
  dplyr::inner_join(gam_mod_top_500_fossil_samples, by = "var") %>%
  dplyr::inner_join(gam_mod_top_1000_fossil_samples, by = "var") %>%
  dplyr::mutate(
    data_combined = purrr::pmap(
      .l = list(
        surface_samples, # ..1
        top_500_yr_fossil_samples, # ..2
        top_1000_yr_fossil_samples # ..3
      ),
      .f = ~ { 
        data <- 
          ..1 %>%
          dplyr::select(lat,
                        surface_samples = var) %>%
          dplyr::inner_join(..2 %>%
                              dplyr::select(lat,
                                            top_500yr_fossil_samples = var), 
                            by = "lat") %>%
          dplyr::inner_join(..3 %>%
                              dplyr::select(lat, 
                                            top_1000yr_fossil_samples = var), 
                            by = "lat")
      }
    )
  ) %>% 
  dplyr::select(var, data_combined)

mpd <- cor_data[1,]$data_combined[[1]]
mntd <- cor_data[2,]$data_combined[[1]]

mpd_correlation <- 
  ggpubr::ggscatter(
    mpd,
    x = "surface_samples",
    y = "top_500yr_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    point = FALSE,
    xlab = "sesMPD\n(surface pollen data)",
    ylab = "sesMPD\n(fossil pollen data)"
    ) +
  font("xy.text", size = 18, color = color_common) +
  font("xlab", size = 20, color = color_common) +
  font("ylab", size = 20, color = color_common)

mntd_correlation <- 
  ggpubr::ggscatter(
    mntd,
    x = "surface_samples",
    y = "top_500yr_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    point = FALSE,
    xlab = "sesMNTD\n(surface pollen data)",
    ylab = "sesMNTD\n(fossil pollen data)"
    ) +
  font("xy.text", size = 18, color = color_common) +
  font("xlab", size = 20, color = color_common) +
  font("ylab", size = 20, color = color_common)


# Compositional turnover ----
# Load predicted models
gam_mod_turnover <- read_rds("Outputs/Data/gam_mod_turnover_050224.rds")
gam_mod_dca <- read_rds("Outputs/Data/gam_mod_dca_060524.rds")

cor_data <- 
  gam_mod_turnover[1,]$predicted_gam[[1]] %>% 
  dplyr::select(lat, surface_samples = var) %>% 
  dplyr::inner_join(gam_mod_turnover[2,]$predicted_gam[[1]] %>% 
                      dplyr::select(lat, top_500yr_fossil_samples = var),
                    by = "lat") %>% 
  dplyr::inner_join(gam_mod_turnover[3,]$predicted_gam[[1]] %>% 
                      dplyr::select(lat, top_1000yr_fossil_samples = var),
                    by = "lat")

cor_data_dca <- 
  gam_mod_dca[1,]$predicted_gam[[1]] %>% 
  dplyr::select(lat, surface_samples = var) %>% 
  dplyr::inner_join(gam_mod_dca[2,]$predicted_gam[[1]] %>% 
                      dplyr::select(lat, top_500yr_fossil_samples = var),
                    by = "lat")

turnover_correlation <- 
  ggpubr::ggscatter(
    cor_data,
    x = "surface_samples",
    y = "top_500yr_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    point = FALSE,
    xlab = "DCCA axis 1\n(surface pollen data)",
    ylab = "DCCA axis 1\n(fossil pollen data)"
    ) +
  font("xy.text", size = 18, color = color_common) +
  font("xlab", size = 20, color = color_common) +
  font("ylab", size = 20, color = color_common)

dca_correlation <- 
  ggpubr::ggscatter(
    cor_data_dca,
    x = "surface_samples",
    y = "top_500yr_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    point = FALSE,
    xlab = "DCA axis 1\n(surface pollen data)",
    ylab = "DCA axis 1\n(fossil pollen data)"
  ) +
  font("xy.text", size = 18, color = color_common) +
  font("xlab", size = 20, color = color_common) +
  font("ylab", size = 20, color = color_common)


final_plot_correlations <- 
  ggpubr::ggarrange(
    mpd_correlation,
    mntd_correlation,
    turnover_correlation,
    labels = c("(a)", "(b)", "(c)"),
    ncol = 3,
    nrow = 1,
    label.y = 1.01,
    font.label = list(size = 22, color = color_common)
    ) 

final_plot_correlations_dca <- 
  ggpubr::ggarrange(
    mpd_correlation,
    mntd_correlation,
    dca_correlation,
    labels = c("(a)", "(b)", "(c)"),
    ncol = 3,
    nrow = 1,
    label.y = 1.01,
    font.label = list(size = 22, color = color_common)
  ) 

ggsave(final_plot_correlations,
       file = paste(
         "Outputs/Figure/Correlations_",
         "surface_samples_vs_fossil_samples_050224.tiff",
         sep = ""),
       height = 10,
       width = 30,
       units = "cm",
       dpi = 400,
       compression = "lzw"
       )

ggsave(final_plot_correlations_dca,
       file = paste(
         "Outputs/Figure/Correlations_",
         "surface_samples_vs_fossil_samples_060524.tiff",
         sep = ""),
       height = 10,
       width = 30,
       units = "cm",
       dpi = 400,
       compression = "lzw"
       )
