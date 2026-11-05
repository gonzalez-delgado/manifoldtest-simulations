# ================================================================
# Power curves of the figures of the paper, in a common style.
#
# Colours, linetypes and shapes are keyed by the DISPLAY label ("Raw dCov"),
# not by the internal column name (Raw_dCov).
#
# The figures are saved close to the size at which they are printed (a panel
# of width 0.48\textwidth is about 2.3 in wide in the article class), so that
# the font sizes below are, up to a factor of about 0.75, those of the page.
# ================================================================

POWER_LABELS <- c(DistCov = "DistCov", BallCov = "BallCov", HHG = "HHG",
                  Raw_dCov = "Raw dCov", Geo_dCov = "Geo dCov", RankCov = "RankCov")
POWER_PAL <- c(DistCov = "#E41A1C", BallCov = "#377EB8", HHG = "#4DAF4A",
               "Raw dCov" = "#984EA3", "Geo dCov" = "#00BFC4", RankCov = "#FF7F00")
POWER_LTY <- c(DistCov = "solid", BallCov = "dashed", HHG = "dotted",
               "Raw dCov" = "dotdash", "Geo dCov" = "twodash", RankCov = "longdash")
POWER_SHP <- c(DistCov = 16, BallCov = 17, HHG = 15, "Raw dCov" = 3,
               "Geo dCov" = 4, RankCov = 18)

theme_power <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = base_size, margin = ggplot2::margin(l = 2, r = 8)),
      legend.key.width = ggplot2::unit(0.55, "cm"),
      legend.key.height = ggplot2::unit(0.4, "cm"),
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(t = 0, b = 0),
      legend.box.margin = ggplot2::margin(t = -4),
      axis.text = ggplot2::element_text(size = base_size - 1),
      axis.title = ggplot2::element_text(size = base_size),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                         size = base_size),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray40",
                                            size = base_size - 2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 10, 6, 4))
}

# Power (%) of the methods `keep` (columns of df) against the column `xvar`.
# legend = "show", "none", or "blank": a blank legend takes the room of a
# visible one, so that a panel placed next to one that carries the legend has
# a plotting area of the same size.
power_plot <- function(df, idvars, xvar, keep, title, subtitle = NULL, xlab,
                       legend = c("show", "none", "blank"),
                       xbreaks = ggplot2::waiver(), legend_rows = 2,
                       base_size = 11) {
  legend <- match.arg(legend)
  m <- reshape2::melt(df, id.vars = idvars, measure.vars = keep,
                      variable.name = "Method", value.name = "Power")
  m$Method <- factor(m$Method, levels = keep, labels = POWER_LABELS[keep])
  lab <- POWER_LABELS[keep]
  # the three legends merge into one, so the blank-legend override is set once
  guide <- ggplot2::guide_legend(nrow = legend_rows)
  guide_col <- ggplot2::guide_legend(
    nrow = legend_rows,
    override.aes = if (legend == "blank") list(alpha = 0) else list())
  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data[[xvar]], y = .data$Power,
                                       color = .data$Method,
                                       linetype = .data$Method,
                                       shape = .data$Method)) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4,
                        alpha = 0.6) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = POWER_PAL[lab]) +
    ggplot2::scale_linetype_manual(values = POWER_LTY[lab]) +
    ggplot2::scale_shape_manual(values = POWER_SHP[lab]) +
    ggplot2::scale_y_continuous(limits = c(0, 104), breaks = seq(0, 100, 20)) +
    ggplot2::scale_x_continuous(breaks = xbreaks) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = "Power (%)",
                  color = NULL, linetype = NULL, shape = NULL) +
    ggplot2::guides(color = guide_col, linetype = guide, shape = guide) +
    theme_power(base_size)
  if (legend == "none") p <- p + ggplot2::theme(legend.position = "none")
  if (legend == "blank")
    p <- p + ggplot2::theme(legend.text = ggplot2::element_text(colour = "transparent"))
  p
}

# Default size: one panel of width 0.48\textwidth
save_power_plot <- function(p, file, w = 3.2, h = 3.0) {
  ggplot2::ggsave(file, p, width = w, height = h, dpi = 300, bg = "white")
  cat("saved", file, "\n")
}

# Several power curves side by side, as one figure: one panel per element of
# `panels`, a shared y axis and a single legend line. Each element is
# list(data, xvar, xlab, breaks), where xlab is a plotmath string shown below
# the panel in place of its x axis title; `title` is common to the panels,
# which are laid out in `nrow` rows.
power_panels <- function(panels, keep, title = NULL, nrow = 1,
                         base_size = 11) {
  labs_x <- vapply(panels, function(pk) pk$xlab, "")
  long <- do.call(rbind, lapply(panels, function(pk) {
    m <- reshape2::melt(pk$data[, c(pk$xvar, keep)], id.vars = pk$xvar,
                        variable.name = "Method", value.name = "Power")
    names(m)[1] <- "x"
    m$panel <- pk$xlab
    m
  }))
  long$panel <- factor(long$panel, levels = labs_x)
  long$Method <- factor(long$Method, levels = keep, labels = POWER_LABELS[keep])
  lab <- POWER_LABELS[keep]
  # the x scales are free: each panel is recognised by the range of its data
  ranges <- lapply(panels, function(pk) range(pk$data[[pk$xvar]]))
  breaks <- function(lim) {
    k <- which.min(vapply(ranges, function(r) sum(abs(r - lim)), 0))
    panels[[k]]$breaks
  }
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$Power,
                                     color = .data$Method,
                                     linetype = .data$Method,
                                     shape = .data$Method)) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4,
                        alpha = 0.6) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = POWER_PAL[lab]) +
    ggplot2::scale_linetype_manual(values = POWER_LTY[lab]) +
    ggplot2::scale_shape_manual(values = POWER_SHP[lab]) +
    ggplot2::scale_y_continuous(limits = c(0, 104), breaks = seq(0, 100, 20)) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::facet_wrap(~ panel, nrow = nrow, scales = "free_x",
                        strip.position = "bottom",
                        labeller = ggplot2::label_parsed) +
    ggplot2::labs(title = title, x = NULL, y = "Power (%)", color = NULL,
                  linetype = NULL, shape = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1),
                    linetype = ggplot2::guide_legend(nrow = 1),
                    shape = ggplot2::guide_legend(nrow = 1)) +
    theme_power(base_size) +
    ggplot2::theme(strip.placement = "outside",
                   strip.text = ggplot2::element_text(size = base_size),
                   panel.spacing = ggplot2::unit(2.5, "lines"))
}
