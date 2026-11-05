# ================================================================
# Power curves of Figures 2 and 3 (hyperbolic space, with and without
# outliers), in a common style.
#
# Colours, linetypes and shapes are keyed by the DISPLAY label ("Raw dCov"),
# not by the internal column name (Raw_dCov).
# ================================================================

POWER_LABELS <- c(DistCov = "DistCov", BallCov = "BallCov", HHG = "HHG",
                  Raw_dCov = "Raw dCov", Geo_dCov = "Geo dCov", RankCov = "RankCov")
POWER_PAL <- c(DistCov = "#E41A1C", BallCov = "#377EB8", HHG = "#4DAF4A",
               "Raw dCov" = "#984EA3", "Geo dCov" = "#00BFC4", RankCov = "#FF7F00")
POWER_LTY <- c(DistCov = "solid", BallCov = "dashed", HHG = "dotted",
               "Raw dCov" = "dotdash", "Geo dCov" = "twodash", RankCov = "longdash")
POWER_SHP <- c(DistCov = 16, BallCov = 17, HHG = 15, "Raw dCov" = 3,
               "Geo dCov" = 4, RankCov = 18)

theme_power <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 15),
      legend.key.width = ggplot2::unit(1.0, "cm"),
      legend.key.height = ggplot2::unit(0.7, "cm"),
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(t = 4, b = 0),
      legend.box.margin = ggplot2::margin(t = -2),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray40", size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 8, 2, 6))
}

# Power (%) of the methods `keep` (columns of df) against the column `xvar`
power_plot <- function(df, idvars, xvar, keep, title, subtitle, xlab,
                       show_legend, xbreaks = ggplot2::waiver()) {
  m <- reshape2::melt(df, id.vars = idvars, variable.name = "Method",
                      value.name = "Power")
  m <- m[m$Method %in% keep, ]
  m$Method <- factor(m$Method, levels = keep, labels = POWER_LABELS[keep])
  lab <- POWER_LABELS[keep]
  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data[[xvar]], y = .data$Power,
                                       color = .data$Method,
                                       linetype = .data$Method,
                                       shape = .data$Method)) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_color_manual(values = POWER_PAL[lab]) +
    ggplot2::scale_linetype_manual(values = POWER_LTY[lab]) +
    ggplot2::scale_shape_manual(values = POWER_SHP[lab]) +
    ggplot2::scale_y_continuous(limits = c(0, 104), breaks = seq(0, 100, 20)) +
    ggplot2::scale_x_continuous(breaks = xbreaks) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4,
                        alpha = 0.5) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = "Power (%)",
                  color = NULL, linetype = NULL, shape = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2),
                    linetype = ggplot2::guide_legend(nrow = 2),
                    shape = ggplot2::guide_legend(nrow = 2)) +
    theme_power()
  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none",
                                            plot.margin = ggplot2::margin(6, 8, 2, 6))
  p
}

save_power_plot <- function(p, file, w = 5.6, h = 5.0) {
  ggplot2::ggsave(file, p, width = w, height = h, dpi = 170)
  cat("saved", file, "\n")
}
