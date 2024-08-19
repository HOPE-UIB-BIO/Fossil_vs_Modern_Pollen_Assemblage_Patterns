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
# Test correlation of latitudinal patterns in PAPs ----
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
  dplyr::ungroup()

turnover_combined <- 
  readr::read_rds("Inputs/Data/turnover_combined_050224.rds") %>% 
  dplyr::group_by(data_type) %>% 
  tidyr::nest(.key = "turnover") %>% 
  dplyr::ungroup() 

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

# A. Phylodiversity ----
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

# B. Turnover ----
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

#--------------------------------------------------------#
# 4. Test the correlation between both data types ----
#--------------------------------------------------------#
# lmodel2: Model II Regression
# Computes model II simple linear regression using ordinary least squares (OLS), 
# major axis (MA), standard major axis (SMA), and ranged major axis (RMA).

# Citation: Legendre, Pierre (2018). lmodel2: Model II Regression. URL:  https://CRAN.R-project.org/package=lmodel2 

# A. Phylodiversity ----
# ses_MPD ----

# Fit major axis model
mod_mpd <- 
  lmodel2::lmodel2(
    ses_mpd_fossil_samples ~ ses_mpd_surface_samples,
    data = correlation_dat_pd,
    range.y = "interval",
    range.x = "interval",
    nperm = 999
    )

# Test different types of major axis models, compare with OLS model
plot_mpd_ols <- 
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

lmodel2::lines.lmodel2(
  mod_mpd,
  method = "MA",
  confidence = TRUE,
  col = "orange",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

lmodel2::lines.lmodel2(
  mod_mpd,
  method = "SMA",
  confidence = TRUE,
  col = "red",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

lmodel2::lines.lmodel2(
  mod_mpd,
  method = "RMA",
  confidence = TRUE,
  col = "darkgreen",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

legend(
  "topleft",
  legend = c("OLS", "MA", "SMA", "RMA"),
  col = c('black', 'orange', 'red', 'darkgreen'),
  lty = c(1, 1, 1, 1),
  lwd = 3,
  bg = "white",
  bty = "n"
) # standard major axis (SMA) and ranged major axis (RMA) gives the best fit


# ses_MNTD ----
mod_mntd <- 
  lmodel2::lmodel2(
    ses_mntd_fossil_samples ~ ses_mntd_surface_samples,
    data = correlation_dat_pd,
    range.y = "interval",
    range.x = "interval",
    nperm = 999
  )

# Test different types of major axis models, compare with OLS model
plot_mntd_ols <- 
  lmodel2::plot.lmodel2(
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

lmodel2::lines.lmodel2(
  mod_mntd,
  method = "MA",
  confidence = TRUE,
  col = "orange",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

lmodel2::lines.lmodel2(
  mod_mntd,
  method = "SMA",
  confidence = TRUE,
  col = "red",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

lines.lmodel2(
  mod_mntd,
  method = "RMA",
  confidence = TRUE,
  col = "darkgreen",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

legend(
  "topleft",
  legend = c("OLS", "MA", "SMA", "RMA"),
  col = c('black', 'orange', 'red', 'darkgreen'),
  lty = c(1, 1, 1, 1),
  lwd = 3,
  bg = "white",
  bty = "n"
) # standard major axis (SMA) gives the best fit

# B. Turnover ----
mod_turnover <- 
  lmodel2::lmodel2(
    turnover_fossil_samples ~ turnover_surface_samples,
    data = correlation_dat_turnover,
    range.y = "relative",
    range.x = "relative",
    nperm = 999
  )

# Test different types of major axis models, compare with OLS model
plot_turnover_ols <- 
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

lmodel2::lines.lmodel2(
  mod_turnover,
  method = "MA",
  confidence = TRUE,
  col = "orange",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

lmodel2::lines.lmodel2(
  mod_turnover,
  method = "SMA",
  confidence = TRUE,
  col = "red",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

lmodel2::lines.lmodel2(
  mod_turnover,
  method = "RMA",
  confidence = TRUE,
  col = "darkgreen",
  xlab = 'Surface pollen data',
  ylab = 'Fossil pollen data',
  main = "",
  cex.lab = 1.5,
  lwd = 50
)

legend(
  "topleft",
  legend = c("OLS", "MA", "SMA", "RMA"),
  col = c('black', 'orange', 'red', 'darkgreen'),
  lty = c(1, 1, 1, 1),
  lwd = 3,
  bg = "white",
  bty = "n"
) # SMA, RMA and MA gives the best fit


#--------------------------------------------------------#
# 5. Plot the models ----
#--------------------------------------------------------#
final_plot_mpd <- 
  ggplot2::ggplot(correlation_dat_pd,
                  aes(x = ses_mpd_surface_samples,
                      y = ses_mpd_fossil_samples)) +
  ggplot2::geom_point() +
  ggpmisc::stat_ma_line(
    method = "MA",
    range.y = "interval",
    range.x = "interval",
    nperm = 999,
    colour = color_common
  ) +
  ggpmisc::stat_ma_eq(
    use_label(c("R2", "P")),
    method = "MA",
    # method = "SMA" doesn't give P-value
    range.y = "interval",
    range.x = "interval",
    nperm = 999,
    size = 5
  ) +
  
  ggplot2::labs(x = "sesMPD\n(surface pollen data)",
                y = "sesMPD\n(fossil pollen data)") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text = element_text(colour = color_common, size = 16),
    axis.title = element_text(colour = color_common, size = 18)
  )

final_plot_mntd <- 
  ggplot2::ggplot(correlation_dat_pd,
                  aes(x = ses_mntd_surface_samples,
                      y = ses_mntd_fossil_samples)) +
  ggplot2::geom_point() +
  ggpmisc::stat_ma_line(
    method = "MA",
    range.y = "interval",
    range.x = "interval",
    nperm = 999,
    colour = color_common
  ) +
  ggpmisc::stat_ma_eq(
    use_label(c("R2", "P")),
    method = "MA",
    # method = "SMA" doesn't give P-value
    range.y = "interval",
    range.x = "interval",
    nperm = 999,
    size = 5
  ) +
  ggplot2::labs(x = "sesMNTD\n(surface pollen data)",
                y = "sesMNTD\n(fossil pollen data)") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text = element_text(colour = color_common, size = 16),
    axis.title = element_text(colour = color_common, size = 18)
  )

final_plot_turnover <- 
  ggplot2::ggplot(
    correlation_dat_turnover,
    aes(x = turnover_surface_samples,
        y = turnover_fossil_samples)
  ) +
  ggplot2::geom_point() +
  ggpmisc::stat_ma_line(
    method = "MA",
    range.y = "relative",
    range.x = "relative",
    nperm = 999,
    colour = color_common
  ) +
  ggpmisc::stat_ma_eq(
    use_label(c("R2", "P")),
    method = "MA",
    # method = "SMA" doesn't give P-value
    range.y = "relative",
    range.x = "relative",
    nperm = 999,
    size = 5
  ) +
  ggplot2::labs(x = "DCCA axis 1\n(surface pollen data)",
                y = "DCCA axis 1\n(fossil pollen data)") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text = element_text(colour = color_common, size = 16),
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

#--------------------------------------------------------#
# 6. Save final figure ----
#--------------------------------------------------------#
ggplot2::ggsave(final_composite_plot,
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
