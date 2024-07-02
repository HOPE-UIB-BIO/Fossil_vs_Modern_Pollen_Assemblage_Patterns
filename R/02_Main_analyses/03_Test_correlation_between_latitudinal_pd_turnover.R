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
phylo_div_full <- 
  readr::read_rds("Inputs/Data/phylo_div_full_090123.rds") %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest(.key = "phylodiv") %>% 
  dplyr::ungroup() %>% 
  dplyr::slice(-c(3,4))

turnover_combined <- 
  readr::read_rds("Inputs/Data/turnover_combined_050224.rds") %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest(.key = "turnover") %>% 
  dplyr::ungroup() %>% 
  dplyr::slice(-3)

dca_combined <- 
  read_rds("Inputs/Data/dca_combined_300424.rds") %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest(.key = "dca") %>% 
  ungroup()

#--------------------------------------------------------#
# 3. Randomly sample data rows ----
#--------------------------------------------------------#
# Number of samples in both data types is highly unequal, randomly sampled 10 samples
# for every 0.25o latitude as follows:
# 1. partitioned entire latitudinal gradient into latitudinal bands of 0.25o latitude
# 2. Only selected latitudinal bands in which data are present for each data type.
# 3. Randomly sampled 10 samples from each latitudinal band for both data types using 
# 'dplyr::slice_sample()'. In this function, f > 10 samples are present in any sample, 
#  10 samples are sampled randomly. But if < 10 samples are present in any band, 
#  existing data are just replicated upto 10 samples. This results in equal number of 
#  samples for both data types in each latitudinal band.

# Phylodiversity
phylodiversity_surface_samples <- 
  phylo_div_full[1,]$phylodiv[[1]] %>% 
  dplyr::arrange(lat) %>% 
  dplyr::mutate(lat_grouped = 
                  ceiling(lat/0.25)
                ) %>% 
  dplyr::mutate(lat_grouped = lat_grouped*0.25) %>% 
  dplyr::select(lat_surface_samples = lat, 
                long_surface_samples = long, 
                climate_zone_revised, 
                ses_mpd_surface_samples = ses_mpd, 
                ses_mntd_surface_samples = ses_mntd, 
                lat_grouped)

phylodiversity_fossil_samples_filtered <- 
  phylo_div_full[2,]$phylodiv[[1]] %>% 
  dplyr::arrange(lat) %>% 
  dplyr::mutate(lat_grouped = 
                  ceiling(lat/0.25)
                ) %>% 
  dplyr::mutate(lat_grouped = lat_grouped*0.25) %>%
  dplyr::select(lat_fossil_samples = lat, 
                long_fossil_samples = long, 
                climate_zone_revised, 
                ses_mpd_fossil_samples = ses_mpd, 
                ses_mntd_fossil_samples = ses_mntd, 
                lat_grouped) %>% 
  dplyr::filter(lat_grouped %in% phylodiversity_surface_samples$lat_grouped)


phylodiversity_surface_samples_filtered <- 
  phylodiversity_surface_samples %>% 
  dplyr::filter(lat_grouped %in% phylodiversity_fossil_samples_filtered$lat_grouped)


set.seed(2330)
phylodiversity_surface_samples_sampled <- 
  phylodiversity_surface_samples_filtered %>% 
  dplyr::group_by(lat_grouped) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(lat_grouped)


set.seed(2330)
phylodiversity_fossil_samples_sampled <- 
  phylodiversity_fossil_samples_filtered %>% 
  dplyr::group_by(lat_grouped) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(lat_grouped)

correlation_dat_pd <- bind_cols(
  phylodiversity_surface_samples_sampled, 
  phylodiversity_fossil_samples_sampled
  )

# Turnover
turnover_surface_samples <- 
  turnover_combined[1,]$turnover[[1]] %>% 
  dplyr::arrange(lat) %>% 
  dplyr::mutate(lat_grouped = 
                  ceiling(lat/0.25)
                ) %>% 
  dplyr::mutate(lat_grouped = lat_grouped*0.25) %>% 
  dplyr::select(lat_surface_samples = lat, 
                long_surface_samples = long, 
                climate_zone_revised, 
                turnover_surface_samples = axis_1, 
                lat_grouped)

turnover_fossil_samples_filtered <- 
  turnover_combined[2,]$turnover[[1]] %>% 
  dplyr::arrange(lat) %>% 
  dplyr::mutate(lat_grouped = 
                  ceiling(lat/0.25)
                ) %>% 
  dplyr::mutate(lat_grouped = lat_grouped*0.25) %>%
  dplyr::select(lat_fossil_samples = lat, 
                long_fossil_samples = long, 
                climate_zone_revised, 
                turnover_fossil_samples = axis_1,
                lat_grouped) %>% 
  dplyr::filter(lat_grouped %in% turnover_surface_samples$lat_grouped)


turnover_surface_samples_filtered <- 
  turnover_surface_samples %>% 
  dplyr::filter(lat_grouped %in% turnover_fossil_samples_filtered$lat_grouped)

set.seed(2330)
turnover_surface_samples_sampled <- 
  turnover_surface_samples_filtered %>% 
  dplyr::group_by(lat_grouped) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(lat_grouped)

set.seed(2330)
turnover_fossil_samples_sampled <- 
  turnover_fossil_samples_filtered %>% 
  dplyr::group_by(lat_grouped) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(lat_grouped)

correlation_dat_turnover <- 
  bind_cols(
  turnover_surface_samples_sampled, 
  turnover_fossil_samples_sampled
  )


# Compositional variation DCA
dca_surface_samples <- 
  dca_combined[1,]$dca[[1]] %>% 
  dplyr::arrange(lat) %>% 
  dplyr::mutate(lat_grouped = 
                  ceiling(lat/0.25)
  ) %>% 
  dplyr::mutate(lat_grouped = lat_grouped*0.25) %>% 
  dplyr::select(lat_surface_samples = lat, 
                long_surface_samples = long, 
                climate_zone_revised, 
                dca_surface_samples = axis_1, 
                lat_grouped)

dca_fossil_samples_filtered <- 
  dca_combined[2,]$dca[[1]] %>% 
  dplyr::arrange(lat) %>% 
  dplyr::mutate(lat_grouped = 
                  ceiling(lat/0.25)
  ) %>% 
  dplyr::mutate(lat_grouped = lat_grouped*0.25) %>%
  dplyr::select(lat_fossil_samples = lat, 
                long_fossil_samples = long, 
                climate_zone_revised, 
                dca_fossil_samples = axis_1,
                lat_grouped) %>% 
  dplyr::filter(lat_grouped %in% dca_surface_samples$lat_grouped)


dca_surface_samples_filtered <- 
  dca_surface_samples %>% 
  dplyr::filter(lat_grouped %in% dca_fossil_samples_filtered$lat_grouped)

set.seed(2330)
dca_surface_samples_sampled <- 
  dca_surface_samples_filtered %>% 
  dplyr::group_by(lat_grouped) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(lat_grouped)

set.seed(2330)
dca_fossil_samples_sampled <- 
  dca_fossil_samples_filtered %>% 
  dplyr::group_by(lat_grouped) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(lat_grouped)

correlation_dat_dca <- 
  bind_cols(
    dca_surface_samples_sampled, 
    dca_fossil_samples_sampled
  )


#--------------------------------------------------------#
# 4. Plot the correlation between both data types
#--------------------------------------------------------#
# lmodel2: Model II Regression
# Computes model II simple linear regression using ordinary least squares (OLS), 
# major axis (MA), standard major axis (SMA), and ranged major axis (RMA).

# Citation: Legendre, Pierre (2018). lmodel2: Model II Regression. URL:  https://CRAN.R-project.org/package=lmodel2 

library(lmodel2)

# A. ses_MPD

# Fit major axis model
mod_mpd <- 
  lmodel2(
    ses_mpd_fossil_samples ~ ses_mpd_surface_samples,
    data = correlation_dat_pd,
    range.y = "interval",
    # range.x and range.y = 'relative' not applicable as there are -ve values
    range.x = "interval",
    nperm = 999
  )

# Test different types of major axis models, compare with OLS model
plot_mpd <- 
  plot.lmodel2(
    mod_mpd,
    method = "OLS",
    confidence = TRUE,
    centroid = FALSE,
    col = "black",
    xlab = 'Surface pollen data',
    ylab = 'Fossil pollen data',
    main = "",
    cex.lab = 1.5
  )

lines.lmodel2(mod_mpd,
              method = "MA",
              confidence = TRUE,
              col = "orange",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)
lines.lmodel2(mod_mpd,
              method = "SMA",
              confidence = TRUE,
              col = "red",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)
lines.lmodel2(mod_mpd,
              method = "RMA",
              confidence = TRUE,
              col = "darkgreen",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)

legend("topleft", 
       legend = c("OLS", "MA", "SMA", "RMA"), 
       col = c('black', 'orange','red','darkgreen'), 
       lty = c(1, 1, 1, 1), 
       lwd = 3, 
       bg = "white", 
       bty = "n") # standard major axis (SMA) gives the best fit


# B. ses_MNTD
mod_mntd <- 
  lmodel2(
    ses_mntd_fossil_samples ~ ses_mntd_surface_samples,
    data = correlation_dat_pd,
    range.y = "interval",
    # range.x and range.y = 'relative' not applicable as there are -ve values
    range.x = "interval",
    nperm = 999
  )

# Test different types of major axis models, compare with OLS model
plot_mntd <- 
  plot.lmodel2(
    mod_mntd,
    method = "OLS",
    confidence = TRUE,
    centroid = FALSE,
    col = "black",
    xlab = 'Surface pollen data',
    ylab = 'Fossil pollen data',
    main = "",
    cex.lab = 1.5
  )

lines.lmodel2(mod_mntd,
              method = "MA",
              confidence = TRUE,
              col = "orange",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)
lines.lmodel2(mod_mntd,
              method = "SMA",
              confidence = TRUE,
              col = "red",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)
lines.lmodel2(mod_mntd,
              method = "RMA",
              confidence = TRUE,
              col = "darkgreen",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)

legend("topleft", 
       legend = c("OLS", "MA", "SMA", "RMA"), 
       col = c('black', 'orange','red','darkgreen'), 
       lty = c(1, 1, 1, 1), 
       lwd = 3, 
       bg = "white", 
       bty = "n") # standard major axis (SMA) gives the best fit

# B. Turnover
mod_turnover <- 
  lmodel2(
    turnover_fossil_samples ~ turnover_surface_samples,
    data = correlation_dat_turnover,
    range.y = "relative",
    range.x = "relative",
    nperm = 999
  )

# Test different types of major axis models, compare with OLS model
plot_turnover <- 
  plot.lmodel2(
    mod_turnover,
    method = "OLS",
    confidence = TRUE,
    centroid = FALSE,
    col = "black",
    xlab = 'Surface pollen data',
    ylab = 'Fossil pollen data',
    main = "",
    cex.lab = 1.5
  )

lines.lmodel2(mod_turnover,
              method = "MA",
              confidence = TRUE,
              col = "orange",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)
lines.lmodel2(mod_turnover,
              method = "SMA",
              confidence = TRUE,
              col = "red",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)
lines.lmodel2(mod_turnover,
              method = "RMA",
              confidence = TRUE,
              col = "darkgreen",
              xlab = 'Surface pollen data',
              ylab = 'Fossil pollen data',
              main = "",
              cex.lab = 1.5,
              lwd = 50)

legend("topleft", 
       legend = c("OLS", "MA", "SMA", "RMA"), 
       col = c('black', 'orange','red','darkgreen'), 
       lty = c(1, 1, 1, 1), 
       lwd = 3, 
       bg = "white", 
       bty = "n") # standard major axis (SMA) gives the best fit

# Make plots of the models using ggplot()
library(ggpmisc)

final_plot_mpd <- 
  ggplot(correlation_dat_pd, 
         aes(x = ses_mpd_surface_samples, 
             y = ses_mpd_fossil_samples)
  ) + 
  geom_point() + 
  stat_ma_line(method = "MA",
               range.y = "interval",
               range.x = "interval",
               nperm = 999,
               colour = color_common) +
  stat_ma_eq(use_label(c("R2", "P")),
             method = "MA", # method = "SMA" doesn't give P-value
             range.y = "interval",
             range.x = "interval",
             nperm = 999,
             size = 5) +
  labs(x = "sesMPD\n(surface pollen data)",
       y = "sesMPD\n(fossil pollen data)"
  ) +
  theme_classic() + 
  theme(axis.text = element_text(colour = color_common, size = 16),
        axis.title = element_text(colour = color_common, size = 18)
  )

final_plot_mntd <- 
  ggplot(correlation_dat_pd, 
         aes(x = ses_mntd_surface_samples, 
             y = ses_mntd_fossil_samples)
  ) + 
  geom_point() + 
  stat_ma_line(method = "MA",
               range.y = "interval",
               range.x = "interval",
               nperm = 999,
               colour = color_common) +
  stat_ma_eq(use_label(c("R2", "P")),
             method = "MA", # method = "SMA" doesn't give P-value
             range.y = "interval",
             range.x = "interval",
             nperm = 999,
             size = 5) +
  labs(x = "sesMNTD\n(surface pollen data)",
       y = "sesMNTD\n(fossil pollen data)"
  ) +
  theme_classic() + 
  theme(axis.text = element_text(colour = color_common, size = 16),
        axis.title = element_text(colour = color_common, size = 18)
  )

final_plot_turnover <- 
  ggplot(correlation_dat_turnover, 
         aes(x = turnover_surface_samples, 
             y = turnover_fossil_samples)
  ) + 
  geom_point() + 
  stat_ma_line(method = "MA",
               range.y = "relative",
               range.x = "relative",
               nperm = 999,
               colour = color_common) +
  stat_ma_eq(use_label(c("R2", "P")),
             method = "MA", # method = "SMA" doesn't give P-value
             range.y = "relative",
             range.x = "relative",
             nperm = 999,
             size = 5) +
  labs(x = "DCCA axis 1\n(surface pollen data)",
       y = "DCCA axis 1\n(fossil pollen data)"
  ) +
  theme_classic() + 
  theme(axis.text = element_text(colour = color_common, size = 16),
        axis.title = element_text(colour = color_common, size = 18)
  )

final_composite_plot <- 
  ggpubr::ggarrange(
    final_plot_mpd,
    final_plot_mntd,
    final_plot_turnover,
    labels = c("(a)", "(b)", "(c)"),
    ncol = 3,
    nrow = 1,
    label.y = 1.01,
    font.label = list(size = 18, color = color_common)
  )

ggsave(final_composite_plot,
       file = paste(
         "Outputs/Figure/Correlations_",
         "surface_samples_vs_fossil_samples_050624.tiff",
         sep = ""),
       height = 12,
       width = 35,
       units = "cm",
       dpi = 400,
       compression = "lzw"
)



#-----------------------------------------------------#

# NOT INCLUDED IN THE REVISED VERSION

#-----------------------------------------------------#
mpd_correlation <- 
  ggpubr::ggscatter(
    correlation_dat_pd,
    x = "ses_mpd_surface_samples",
    y = "ses_mpd_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    point = TRUE,
    size = 0.75,
    xlab = "sesMPD\n(surface pollen data)",
    ylab = "sesMPD\n(fossil pollen data)"
  ) +
  font("xy.text", size = 20, color = color_common) +
  font("xlab", size = 22, color = color_common) +
  font("ylab", size = 22, color = color_common)

mntd_correlation <- 
  ggpubr::ggscatter(
    correlation_dat_pd,
    x = "ses_mntd_surface_samples",
    y = "ses_mntd_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    point = TRUE,
    size = 0.75,
    xlab = "sesMNTD\n(surface pollen data)",
    ylab = "sesMNTD\n(fossil pollen data)"
  ) +
  font("xy.text", size = 20, color = color_common) +
  font("xlab", size = 22, color = color_common) +
  font("ylab", size = 22, color = color_common)
  
turnover_correlation <- 
  ggpubr::ggscatter(
    correlation_dat_turnover,
    x = "turnover_surface_samples",
    y = "turnover_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    size = 0.75,
    cor.coeff.args = list(method = "pearson", label.y = 2),
    point = TRUE,
    xlab = "DCCA axis 1\n(surface pollen data)",
    ylab = "DCCA axis 1\n(fossil pollen data)"
  ) +
  font("xy.text", size = 20, color = color_common) +
  font("xlab", size = 22, color = color_common) +
  font("ylab", size = 22, color = color_common)

dca_correlation <- 
  ggpubr::ggscatter(
    correlation_dat_dca,
    x = "dca_surface_samples",
    y = "dca_fossil_samples",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.size = 6,
    size = 0.75,
    cor.coeff.args = list(method = "pearson"),
    point = TRUE,
    xlab = "DCA axis 1\n(surface pollen data)",
    ylab = "DCA axis 1\n(fossil pollen data)"
  ) +
  font("xy.text", size = 20, color = color_common) +
  font("xlab", size = 22, color = color_common) +
  font("ylab", size = 22, color = color_common)

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
         "surface_samples_vs_fossil_samples_210524.tiff",
         sep = ""),
       height = 12,
       width = 35,
       units = "cm",
       dpi = 400,
       compression = "lzw"
       )

ggsave(final_plot_correlations_dca,
       file = paste(
         "Outputs/Figure/Correlations_",
         "surface_samples_vs_fossil_samples_dca_210524.tiff",
         sep = ""),
       height = 12,
       width = 35,
       units = "cm",
       dpi = 400,
       compression = "lzw"
       )
