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
# Spatial variation in phylogenetic dispersion and compositional turnover ----
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
# 3. Make maps of spatial variation in PAPs ----
#--------------------------------------------------------#

# Base map ----
base_map <-
  phylo_div_full %>%
  dplyr::filter(data_type == "surface_samples") %>%
  ggplot2::ggplot(aes(x = long, y = lat)) +
  ggplot2::coord_fixed(ylim = c(20.00, 80.00),
                       xlim = c(65.00, 135.00)) +
  ggplot2::labs(x = expression(paste('Longitude ', (degree ~ E))),
                y = expression(paste('Latitude ', (degree ~ N))),
                fill = "Climate-zones") +
  ggplot2::theme_classic() +
  ggplot2::borders(colour = color_common,
                   linewidth = 0.2) +
  ggplot2::theme(
    axis.title = element_text(color = color_common, size = 16),
    axis.text = element_text(colour = color_common, size = 12),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) 

# A. Surface samples ----
surface_samples_pd <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "surface_samples")
surface_samples_turnover <- 
  turnover_combined[1,]$turnover[[1]] 

# B. Fossil samples ----
top_500_fossil_samples_pd <- 
  phylo_div_full %>% 
  dplyr::filter(data_type == "top_500_yr")
top_500_fossil_turnover <-  
  turnover_combined[2,]$turnover[[1]] 


# MAP sesMPD ----
surface_samples_mpd <- 
  base_map +
  ggplot2::geom_point(aes(size = ses_mpd,
                          colour = climate_zone_revised),
                      data = surface_samples_pd) +
  ggplot2::scale_colour_manual(values = my_palette) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::scale_size(range = c(-4, 5)) +
  ggplot2::labs(size = "sesMPD",
                colour = "Climate zone",
                fill = "Climate zone") +
  ggplot2::theme(legend.box = "horizontal") +
  ggplot2::guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
  )

top_500_fossil_samples_mpd <- 
  base_map +
  ggplot2::geom_point(aes(size = ses_mpd,
                          colour = climate_zone_revised),
                      data = top_500_fossil_samples_pd) +
  ggplot2::scale_colour_manual(values = my_palette) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::scale_size(range = c(-4, 5)) +
  ggplot2::labs(colour = "Climate zone",
                fill = "Climate zone",
                size = "sesMPD") +
  ggplot2::theme(legend.box = "horizontal") +
  ggplot2::guides(
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

# Map sesMNTD ----
surface_samples_mntd <- 
  base_map +
  ggplot2::geom_point(aes(size = ses_mntd,
                          colour = climate_zone_revised),
                      data = surface_samples_pd) +
  ggplot2::scale_colour_manual(values = my_palette) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::scale_size(range = c(-4, 5)) +
  ggplot2::labs(colour = "Climate zone",
                fill = "Climate zone",
                size = "sesMNTD") +
  ggplot2::theme(legend.box = "horizontal") +
  ggplot2::guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
  )

top_500_fossil_samples_mntd <- 
  base_map +
  ggplot2::geom_point(aes(size = ses_mntd,
                          colour = climate_zone_revised),
                      data = top_500_fossil_samples_pd) +
  ggplot2::scale_colour_manual(values = my_palette) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::scale_size(range = c(-4, 5)) +
  ggplot2::labs(colour = "Climate zone",
                fill = "Climate zone",
                size = "sesMNTD") +
  ggplot2::theme(legend.box = "horizontal") +
  ggplot2::guides(
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

# Map compositional turnover ----
surface_sample_turnover <- 
  base_map +
  ggplot2::geom_point(aes(size = axis_1,
                          colour = climate_zone_revised),
                      data = surface_samples_turnover) +
  ggplot2::scale_colour_manual(values = my_palette) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::scale_size(range = c(-4, 5)) +
  ggplot2::labs(colour = "Climate zone",
                fill = "Climate zone",
                size = "DCCA axis 1") +
  ggplot2::theme(legend.box = "horizontal") +
  ggplot2::guides(
    size =  guide_legend(order = 1),
    fill  = guide_legend(order = 2),
    colour  = guide_legend(order = 2)
  )

top_500_fossil_sample_turnover <- 
  base_map +
  ggplot2::geom_point(aes(size = axis_1,
                          colour = climate_zone_revised),
                      data = top_500_fossil_turnover) +
  ggplot2::scale_colour_manual(values = my_palette) +
  ggplot2::scale_fill_manual(values = my_palette) +
  ggplot2::scale_size(range = c(0, 5)) +
  ggplot2::labs(colour = "Climate zone",
                fill = "Climate zone",
                size = "DCCA axis 1") +
  ggplot2::theme(legend.box = "horizontal") +
  ggplot2::guides(
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
