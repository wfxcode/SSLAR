plot_spaniel <- function(data_df, grob, x, y, point_colour, point_size){

  # Inverse Y to flip the coordinates
  data_df$y_inv <- 36 - data_df$y

  data_df[, point_size] <- if_else(data_df[, point_size] == 0, NA_real_, data_df[, point_size])

  tmp_plt <- ggplot2::ggplot(data_df,
                  ggplot2::aes_string(x, "y_inv",
                                      color = point_colour,
                                      size = point_size)) +
    ggplot2::xlim(1, 33) +
    ggplot2::ylim(1, 35) +
    # Layer 1 - Plot image
    ggplot2::annotation_custom(grob,
                               xmin = 1,
                               xmax = 33,
                               ymin = 1,
                               ymax = 35) +
    # Layer 2 - Plot points
    geom_point() +
    # Layer 3 - Join legends all with same name
    labs(color = "Proportion",
         size = "Proportion"
         # alpha = "Proportion"
         ) +
    # Coordinates fixed so x and y keep the proportions
    coord_fixed(1) +
    # Tune colour parameters
    ggplot2::scale_size_continuous(range=c(0, 3), limits = c(0, 1)) +
    ggplot2::scale_color_gradientn(
      colours = heat.colors(10, rev = TRUE),
      limits = c(0, 1)) +
    # Join legends into one
    ggplot2::guides(color = ggplot2::guide_legend(),
                    size = ggplot2::guide_legend()
                    )

  return(tmp_plt)
}

plot_spaniel_B <- function(data_df, grob, x, y, point_colour, point_size){

  # Inverse Y to flip the coordinates
  data_df$y_inv <- 36-data_df$y

  data_df[, point_size] <- if_else(data_df[, point_size] == 0, NA_real_, data_df[, point_size])

  tmp_plt <- ggplot2::ggplot(data_df,
                             ggplot2::aes_string(x, "y_inv",
                                                 color = point_colour,
                                                 size = point_size)) +
    ggplot2::xlim(min(data_df$x)-1, max(data_df$x)+1) +
    ggplot2::ylim(min(data_df$y_inv)-1, max(data_df$y_inv)+1) +
    # Layer 1 - Plot image
    ggplot2::annotation_custom(grob,
                               xmin = min(data_df$x)-1,
                               xmax = max(data_df$x)+1,
                               ymin = min(data_df$y_inv)-1,
                               ymax = max(data_df$y_inv)+1) +
    # Layer 2 - Plot points
    geom_point() +
    # Layer 3 - Join legends all with same name
    labs(color = "Proportion",
         size = "Proportion"
         # alpha = "Proportion"
    ) +
    # Coordinates fixed so x and y keep the proportions
    coord_fixed(1) +
    # Tune colour parameters
    ggplot2::scale_size_continuous(range=c(0, 5), limits = c(0, 1)) +
    ggplot2::scale_color_gradientn(
      colours = heat.colors(10, rev = TRUE),
      limits = c(0, 1)) +
    # Join legends into one
    ggplot2::guides(color = ggplot2::guide_legend(),
                    size = ggplot2::guide_legend()
    )

  return(tmp_plt)
}

