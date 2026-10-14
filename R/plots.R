# ---- Figures (ggplot2) ----------------------------------------------------------
# wright_map_sqi() and plot_expected_scores() as in the analysis scripts of the
# project, with the axis labels and breaks made arguments.

#' Person-item targeting plot (Wright map) on the SQI scale, in the style of
#' easyRasch2::RMtargeting(): person histogram (top) and one bar per item split
#' into response-category bands (bottom). The x axis is the digital measure
#' (SQI) instead of logits.
#'
#' @param person_sqi  person locations on the digital scale
#' @param thr         data.frame with item, k, sqi_est (threshold on the digital scale)
#' @param item_labels labels for the items, in the order of unique(thr$item)
#' @param note        caption text placed under the plot
#' @param xlim x range (widened to include persons and thresholds)
#' @param binwidth histogram bin width
#' @param category_labels labels of the four response categories
#' @param person_fill histogram colour
#' @param x_label x-axis title
#' @param x_breaks x-axis breaks
#' @return a patchwork object
#' @export
wright_map_sqi <- function(person_sqi, thr, item_labels,
                           xlim = c(15, 90), binwidth = 2,
                           category_labels = c("0 (no difficulty)", "1", "2",
                                               "3 (severe difficulty)"),
                           person_fill = "#0072B2", note = NULL,
                           x_label = "Actigraphy SQI (digital measure; higher = better sleep)",
                           x_breaks = seq(0, 100, by = 10)) {
  xlim  <- range(c(xlim, person_sqi, thr$sqi_est))
  if (!is.null(note))   # wrap each caption line to fit the figure width
    note <- paste(unlist(lapply(strsplit(note, "\n")[[1]], strwrap, width = 120)),
                  collapse = "\n")
  cats  <- factor(0:3, labels = category_labels)
  pal   <- setNames(viridisLite::mako(4, begin = 0.9, end = 0.2), category_labels)
  item_ids <- unique(thr$item)
  n_it  <- length(item_ids)
  shared_x <- scale_x_continuous(limits = xlim, expand = c(0, 0),
                                 breaks = x_breaks,
                                 oob = scales::oob_keep)
  base_theme <- theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, colour = "grey30"),
          plot.title = element_text(size = 11, hjust = 0.5, colour = "grey25"))
  no_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                axis.ticks.x = element_blank())

  # Person histogram with mean and +-1 SD
  m <- mean(person_sqi); s <- sd(person_sqi)
  p1 <- ggplot(data.frame(x = person_sqi), aes(x)) +
    annotate("rect", xmin = m - s, xmax = m + s, ymin = -Inf, ymax = Inf,
             fill = person_fill, alpha = 0.12) +
    geom_histogram(binwidth = binwidth, boundary = 0, fill = person_fill,
                   colour = "white") +
    geom_vline(xintercept = m, linetype = "dashed", linewidth = 0.8,
               colour = "grey20") +
    shared_x +
    labs(y = "Persons", title = sprintf("Mean = %.1f, SD = %.1f", m, s)) +
    base_theme + no_x

  # Category bands. GRM thresholds are ordered and mu() falls with the latent
  # trait, so on the SQI axis the bands run 3 | 2 | 1 | 0 from left to right.
  # The open-ended outer bands fade towards the panel edge.
  bands <- do.call(rbind, lapply(seq_len(n_it), function(j) {
    tj   <- thr[thr$item == item_ids[j], ]
    cuts <- tj$sqi_est[order(tj$k, decreasing = TRUE)]   # k = 3, 2, 1
    y    <- n_it - j + 1
    inner <- data.frame(xmin = cuts[1:2], xmax = cuts[2:3], cat = c(2, 1),
                        alpha = 1)
    fade <- function(from, to, cat) {
      s <- seq(from, to, length.out = 41)
      data.frame(xmin = head(s, -1), xmax = tail(s, -1), cat = cat,
                 alpha = seq(1, 0.1, length.out = 40))
    }
    out <- rbind(fade(cuts[1], xlim[1], 3), inner, fade(cuts[3], xlim[2], 0))
    lo <- pmin(out$xmin, out$xmax); out$xmax <- pmax(out$xmin, out$xmax)
    out$xmin <- lo
    out$y <- y
    out
  }))
  bands$cat <- cats[bands$cat + 1]

  p3 <- ggplot() +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = y - 0.3,
                                ymax = y + 0.3, fill = cat, alpha = alpha)) +
    scale_alpha_identity() +
    shared_x +
    scale_y_continuous(breaks = n_it:1, labels = item_labels,
                       expand = expansion(add = 0.5)) +
    scale_fill_manual(values = pal, limits = category_labels, drop = FALSE,
                      name = "Response category") +
    labs(x = x_label, y = NULL) +
    base_theme +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  (p1 / p3) +
    plot_layout(heights = c(3, 5)) +
    plot_annotation(caption = note,
                    theme = theme(plot.caption = element_text(hjust = 0, size = 9.5)))
}

#' Expected item scores against the digital measure
#'
#' Structural and probabilistic crosswalks, estimated (solid) and true
#' (dashed, columns ending in `_true`), with a rug of the observed measure.
#' @param cw crosswalk table from [crosswalk()] with `_true` columns added
#' @param sqi observed digital measure
#' @param items item names
#' @param x_label x-axis title
#' @export
plot_expected_scores <- function(cw, sqi, items = paste0("C", 1:7),
                                 x_label = "Actigraphy SQI (digital measure)") {
  panels <- c(Estar_ = "Structural: E[C | eta*(SQI)]",
              E_     = "Probabilistic: E[C | SQI], integrating p(eta | SQI)")
  long <- do.call(rbind, lapply(names(panels), function(prefix) {
    do.call(rbind, lapply(items, function(i) rbind(
      data.frame(SQI = cw$SQI, score = cw[[paste0(prefix, i)]], item = i,
                 source = "estimated", panel = panels[[prefix]]),
      data.frame(SQI = cw$SQI, score = cw[[paste0(prefix, i, "_true")]], item = i,
                 source = "true", panel = panels[[prefix]]))))
  }))
  long$panel <- factor(long$panel, levels = panels)
  ggplot(long, aes(SQI, score, colour = item, linetype = source,
                   linewidth = source)) +
    geom_line(na.rm = TRUE) +
    geom_rug(data = data.frame(SQI = sqi), aes(x = SQI), inherit.aes = FALSE,
             sides = "b", alpha = 0.15) +
    facet_wrap(~ panel) +
    scale_linetype_manual(values = c(estimated = "solid", true = "dashed"),
                          name = NULL) +
    scale_linewidth_manual(values = c(estimated = 0.9, true = 0.5), name = NULL) +
    scale_colour_manual(values = setNames(hcl.colors(length(items), "Dark 3"), items),
                        name = "Component") +
    coord_cartesian(xlim = range(cw$SQI), ylim = c(0, 3)) +
    labs(x = x_label, y = "Expected item score (0-3)") +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, colour = "grey30"),
          legend.position = "bottom")
}
