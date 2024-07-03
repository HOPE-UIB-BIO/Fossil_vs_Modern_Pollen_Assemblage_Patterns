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

gam_mod_surface_samples_pd <- 
  readr::read_rds("Outputs/Data/gam_mod_surface_samples_100124.rds")
gam_mod_top_500_fossil_samples_pd <- 
  readr::read_rds("Outputs/Data/gam_mod_top_500_fossil_samples_100124.rds")
gam_mod_turnover <- 
  readr::read_rds("Outputs/Data/gam_mod_turnover_050224.rds")

#--------------------------------------------------------#
# 3. Plot GAM models ----
#--------------------------------------------------------#
#---------------------------------#
# A. Surface samples ----
#---------------------------------#

# Phylogenetic dispersion ----
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
        
        ggplot2::theme(legend.position = "none",
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


# Compositional turnover ----
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
        
        ggplot2::theme(legend.position = "none",
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

fig_1 <- 
  plot_gam_mod_surface_samples_pd[[1]] + 
  ggplot2::guides(fill = guide_legend(nrow = 3,
                             title.position = "top"),
         colour = guide_legend(nrow = 3,
                               title.position = "top")
         ) +
  ggplot2::theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
        )
fig_2 <- plot_gam_mod_surface_samples_pd[[2]] 
fig_3 <- plot_gam_mod_turnover[[1]] 

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

#---------------------------------#
# B. Top 500-yr fossil pollen samples ----
#---------------------------------#
# Phylogenetic dispersion ----
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
        
        ggplot2::theme(legend.position = "none",
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


# Compositional turnover ----

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
        ggplot2::theme(legend.position = "none",
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

fig_4 <- 
  plot_gam_500_pd[[1]] +
  ggplot2::guides(fill = guide_legend(nrow = 3,
                             title.position = "top"),
         colour = guide_legend(nrow = 3,
                               title.position = "top")
  ) +
  ggplot2::theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
  )
fig_5 <- plot_gam_500_pd[[2]]
fig_6 <- plot_gam_mod_turnover_top_500[[1]] 

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
  ggpubr::annotate_figure(
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

#---------------------------------#
# 4. Save final figure ----
#---------------------------------#
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
