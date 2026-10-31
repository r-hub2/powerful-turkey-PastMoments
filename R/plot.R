## =============================================================================
## Private panel renderers
## =============================================================================

.axis_text_layout <- function(cex.lab, cex.axis, minimum_margin = 4.2) {
    values <- c(cex.lab = cex.lab, cex.axis = cex.axis)
    if (!is.numeric(values) || any(!is.finite(values)) || any(values <= 0)) {
        stop("cex.lab and cex.axis must be strictly positive finite scalars.")
    }

    # `mgp` and `mar` are expressed in text lines, whereas cex.lab and
    # cex.axis enlarge the glyphs without moving their baselines. Increase both
    # the baseline separation and the available margin as the text grows.
    # Keep horizontal tick labels fully below the axis. Values below one text
    # line make enlarged glyphs straddle the axis rather than clear it.
    axis_line <- max(1, 0.7 * unname(values[["cex.axis"]]))
    label_line <- axis_line + max(
        1.45,
        0.85 * unname(values[["cex.axis"]]) +
            0.55 * unname(values[["cex.lab"]])
    )
    margin <- max(
        minimum_margin,
        label_line + 0.85 * unname(values[["cex.lab"]])
    )

    list(mgp = c(label_line, axis_line, 0), margin = margin)
}

.geological_axis_margin <- function(text_layout, cex.axis, periods) {
    n_scales <- if (is.list(periods) && !is.data.frame(periods)) {
        length(periods)
    } else {
        1L
    }
    max(
        7.2,
        text_layout$margin + 2 + 1.5 * cex.axis * max(1L, n_scales)
    )
}

.layout_heights_with_fixed_margins <- function(weights, margin_lines) {
    line_height <- graphics::par("csi") * graphics::par("mex")
    layout_height <- graphics::par("din")[2L] -
        sum(graphics::par("omi")[c(1L, 3L)])
    fixed_height <- margin_lines * line_height
    available_height <- layout_height - sum(fixed_height)
    if (!is.finite(available_height) || available_height <= 0) {
        return(weights)
    }
    fixed_height + available_height * weights / sum(weights)
}

## Plot the reconstructed empirical mean
##
## The solid curve is the conditional expectation. Since the empirical mean is
## conditionally Gaussian, its median is identical to its expectation.
##
## @param x A reconstruction object containing `empirical_mean`.
## @param clades Optional clade names to display.
## @param colors,palette Clade colours or an HCL palette name. When no clade is
##   specified, black is used by default.
## @param interval How to draw an available conditional quantile interval.
## @param interval_probs Two available probabilities delimiting the interval.
## @param interval_alpha Transparency of a filled interval.
## @param halo_style Whether the filled interval is solid or hatched.
## @param hatch_density,hatch_angle,hatch_lwd Hatching controls.
## @param hatched_outline Whether to outline a hatched region.
## @param hatched_outline_lty,hatched_outline_lwd,hatched_outline_alpha Style of
##   the hatched-region outline.
## @param time_start,xlim,ylim Displayed time and value ranges.
## @param xlab,ylab Axis labels.
## @param time_orientation Whether time is measured from the root or before present.
## @param present_time Present time used for the reverse time scale.
## @param draw_x_axis Whether to draw the horizontal axis.
## @param bands,n_bands,band_color Optional alternating time bands.
## @param legend_pos Legend position, or `FALSE`.
## @param expectation_lwd Width of the expectation curve.
## @param cex.lab,cex.axis,cex.legend Text scaling parameters.
## @param ... Additional graphical arguments.
## @return An invisible list containing plotting metadata.
 .plot_empirical_mean_panel <- function(
    x,
    clades = NULL,
    colors = NULL,
    palette = "Dark 3",
    interval = c("ribbon", "lines", "both", "none"),
    interval_probs = NULL,
    interval_alpha = 0.14,
    halo_style = c("solid", "hatched"),
    hatch_density = 18,
    hatch_angle = 45,
    hatch_lwd = 0.7,
    hatched_outline = TRUE,
    hatched_outline_lty = 2,
    hatched_outline_lwd = 0.8,
    hatched_outline_alpha = 0.75,
    time_start = NULL,
    xlim = NULL,
    ylim = NULL,
    xlab = NULL,
    ylab = "Empirical mean",
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    draw_x_axis = TRUE,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    geological_periods = NULL,
    geological_present_time = NULL,
    geological_background_alpha = 0.12,
    legend_pos = "topright",
    show_expectation = TRUE,
    expectation_lwd = 2,
    cex.lab = 1.35,
    cex.axis = 1.0,
    cex.legend = 1.0,
    ...
) {
    time_orientation <- match.arg(time_orientation)
    if (is.null(xlab)) {
        xlab <- if (time_orientation == "from_root") "Time from root" else "Time before present"
    }

    .plot_empirical_quantity(
        x = x,
        quantity = "empirical_mean",
        clades = clades,
        colors = colors,
        palette = palette,
        interval = interval,
        interval_probs = interval_probs,
        interval_alpha = interval_alpha,
        halo_style = halo_style,
        hatch_density = hatch_density,
        hatch_angle = hatch_angle,
        hatch_lwd = hatch_lwd,
        hatched_outline = hatched_outline,
        hatched_outline_lty = hatched_outline_lty,
        hatched_outline_lwd = hatched_outline_lwd,
        hatched_outline_alpha = hatched_outline_alpha,
        show_expectation = show_expectation,
        show_median = FALSE,
        expectation_lwd = expectation_lwd,
        time_start = time_start,
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        time_orientation = time_orientation,
        present_time = present_time,
        draw_x_axis = draw_x_axis,
        bands = bands,
        n_bands = n_bands,
        band_color = band_color,
        geological_periods = geological_periods,
        geological_present_time = geological_present_time,
        geological_background_alpha = geological_background_alpha,
        legend_pos = legend_pos,
        cex.lab = cex.lab,
        cex.axis = cex.axis,
        cex.legend = cex.legend,
        ...
    )
}



## Plot the empirical mean with expected empirical spread
##
## The envelope is the conditional empirical-mean expectation plus or minus a
## selected measure of expected empirical spread among lineages. It is not an
## uncertainty interval for the empirical mean. `spread = "expected_sd"` uses
## the expected empirical standard deviation, `E[sqrt(V) | data]`.
##
## @param x A reconstruction object containing both empirical quantities.
## @param clades Optional clade names to display.
## @param colors,palette Clade colours or an HCL palette name.
## @param envelope How to draw the empirical-spread envelope.
## @param spread Quantity used as the empirical spread around the expected
##   empirical mean. `"expected_sd"` uses the expected empirical standard
##   deviation, and `"none"` suppresses the envelope.
## @param envelope_alpha Transparency of a filled envelope.
## @param halo_style Whether the filled envelope is solid or hatched.
## @param hatch_density,hatch_angle,hatch_lwd Hatching controls.
## @param hatched_outline Whether to outline a hatched region.
## @param hatched_outline_lty,hatched_outline_lwd,hatched_outline_alpha Style of
##   the hatched-region outline.
## @param envelope_lwd,envelope_lty Style of envelope boundary lines.
## @param time_start,xlim,ylim Displayed time and value ranges.
## @param xlab,ylab Axis labels.
## @param time_orientation Whether time is measured from the root or before present.
## @param present_time Present time used for the reverse time scale.
## @param draw_x_axis Whether to draw the horizontal axis.
## @param bands,n_bands,band_color Optional alternating time bands.
## @param legend_pos Legend position, or `FALSE`.
## @param spread_legend_label Optional text identifying the empirical-spread
##   envelope in the clade legend. `NULL` adds no envelope entry.
## @param expectation_lwd Width of the expectation curve.
## @param cex.lab,cex.axis,cex.legend Text scaling parameters.
## @param ylim_expand Relative expansion of an automatically computed y range.
## @param ... Additional graphical arguments.
## @return An invisible list containing plotting metadata and merged data.
.plot_empirical_spread_panel <- function(
    x,
    clades = NULL,
    colors = NULL,
    palette = "Dark 3",
    envelope = c("ribbon", "lines", "both", "none"),
    spread = c("expected_sd", "none"),
    envelope_alpha = 0.14,
    halo_style = c("hatched", "solid"),
    hatch_density = 18,
    hatch_angle = 45,
    hatch_lwd = 0.7,
    hatched_outline = TRUE,
    hatched_outline_lty = 2,
    hatched_outline_lwd = 0.8,
    hatched_outline_alpha = 0.75,
    envelope_lwd = 0.8,
    envelope_lty = 3,
    time_start = NULL,
    xlim = NULL,
    ylim = NULL,
    xlab = NULL,
    ylab = "Empirical mean",
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    draw_x_axis = TRUE,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    geological_periods = NULL,
    geological_present_time = NULL,
    geological_background_alpha = 0.12,
    legend_pos = "topright",
    spread_legend_label = NULL,
    show_expectation = TRUE,
    expectation_lwd = 2,
    cex.lab = 1.35,
    cex.axis = 1.0,
    cex.legend = 1.0,
    ylim_expand = 0.06,
    ...
) {
    .validate_reconstruction_plot_data(x)
    envelope <- match.arg(envelope)
    spread <- match.arg(spread)
    halo_style <- match.arg(halo_style)
    time_orientation <- match.arg(time_orientation)
    if (!is.null(spread_legend_label) &&
        (!is.character(spread_legend_label) ||
         length(spread_legend_label) != 1L ||
         is.na(spread_legend_label))) {
        stop("spread_legend_label must be NULL or a single character string.")
    }
    if (!is.logical(hatched_outline) || length(hatched_outline) != 1L || is.na(hatched_outline)) {
        stop("hatched_outline must be TRUE or FALSE.")
    }
    if (!is.numeric(hatched_outline_lty) || length(hatched_outline_lty) != 1L || !is.finite(hatched_outline_lty)) {
        stop("hatched_outline_lty must be a finite scalar.")
    }
    if (!is.numeric(hatched_outline_lwd) || length(hatched_outline_lwd) != 1L || !is.finite(hatched_outline_lwd) || hatched_outline_lwd <= 0) {
        stop("hatched_outline_lwd must be a strictly positive finite scalar.")
    }
    if (!is.numeric(hatched_outline_alpha) || length(hatched_outline_alpha) != 1L || !is.finite(hatched_outline_alpha) || hatched_outline_alpha < 0 || hatched_outline_alpha > 1) {
        stop("hatched_outline_alpha must be a finite scalar in [0, 1].")
    }

    d_mean <- .subset_reconstruction_plot_data(
        x, "empirical_mean", clades = clades
    )
    clade_names <- unique(d_mean$clade)
    d_var <- .subset_reconstruction_plot_data(
        x, "empirical_variance", clades = clade_names
    )

    mean_part <- d_mean[, c("time", "clade", "n_lineages", "expectation"),
                        drop = FALSE]
    names(mean_part)[names(mean_part) == "n_lineages"] <- "n_lineages_mean"
    names(mean_part)[names(mean_part) == "expectation"] <- "mean_expectation"

    var_columns <- c("time", "clade", "n_lineages", "expectation")
    if ("expected_sd" %in% names(d_var)) {
        var_columns <- c(var_columns, "expected_sd")
    }
    var_part <- d_var[, var_columns, drop = FALSE]
    names(var_part)[names(var_part) == "n_lineages"] <- "n_lineages_variance"
    names(var_part)[names(var_part) == "expectation"] <- "empirical_variance_expectation"

    d <- merge(
        mean_part, var_part,
        by = c("time", "clade"),
        all = FALSE,
        sort = FALSE
    )
    if (nrow(d) == 0L) {
        stop("No matching mean and variance reconstruction rows are available.")
    }

    mismatch <- is.finite(d$n_lineages_mean) &
        is.finite(d$n_lineages_variance) &
        d$n_lineages_mean != d$n_lineages_variance
    if (any(mismatch)) {
        warning(
            "The mean and variance results have different lineage counts for ",
            sum(mismatch), " matched rows."
        )
    }

    if (spread == "expected_sd" && !"expected_sd" %in% names(d)) {
        stop(
            "spread = 'expected_sd' requires a reconstruction containing ",
            "the expected_sd column. Recompute the reconstruction with the ",
            "current version of PastMoments."
        )
    }
    d$envelope <- switch(
        spread,
        expected_sd = d$expected_sd,
        none = 0
    )
    d$lower <- d$mean_expectation - d$envelope
    d$upper <- d$mean_expectation + d$envelope
    d$clade <- factor(d$clade, levels = clade_names)
    d <- d[order(d$clade, d$time), , drop = FALSE]
    d$clade <- as.character(d$clade)
    rownames(d) <- NULL

    cols <- make_clade_colors(clade_names, colors = colors, palette = palette)

    xlim <- .resolve_time_window(d$time, time_start = time_start, xlim = xlim)
    d <- .trim_plot_data_to_window(d, xlim)

    if (is.null(ylim)) {
        values <- d$mean_expectation
        if (envelope != "none" && spread != "none") {
            values <- c(values, d$lower, d$upper)
        }
        values <- values[is.finite(values)]
        if (length(values) == 0L) {
            stop("No finite value is available for plotting.")
        }
        ylim <- range(values)
        if (ylim[1L] == ylim[2L]) {
            delta <- max(1, abs(ylim[1L])) * 0.05
            ylim <- ylim + c(-delta, delta)
        } else {
            ylim <- grDevices::extendrange(ylim, f = ylim_expand)
        }
    }

    if (is.null(xlab)) {
        xlab <- if (time_orientation == "from_root") {
            "Time from root"
        } else {
            "Time before present"
        }
    }

    bands <- .resolve_plot_bands(xlim, bands, n_bands)

    graphics::plot(
        xlim, ylim,
        type = "n",
        axes = FALSE,
        frame.plot = FALSE,
        xlab = if (isTRUE(draw_x_axis)) xlab else "",
        ylab = ylab,
        cex.lab = cex.lab,
        cex.axis = cex.axis,
        xaxs = "i",
        ...
    )

    .draw_geological_background(
        geological_periods,
        present_time = geological_present_time,
        alpha = geological_background_alpha
    )
    .draw_time_bands(bands, band_color)

    for (nm in clade_names) {
        dc <- d[d$clade == nm, , drop = FALSE]
        dc <- dc[order(dc$time), , drop = FALSE]
        col <- unname(cols[[nm]])

        if (spread != "none" && envelope %in% c("ribbon", "both")) {
            .plot_finite_runs(
                dc$time, dc$lower, dc$upper,
                color = grDevices::adjustcolor(col, alpha.f = envelope_alpha),
                halo_style = halo_style,
                hatch_density = hatch_density,
                hatch_angle = hatch_angle,
                hatch_lwd = hatch_lwd
            )
            if (halo_style == "hatched" && isTRUE(hatched_outline) &&
                envelope %in% c("ribbon", "both")) {
                graphics::lines(
                    dc$time, dc$lower,
                    col = grDevices::adjustcolor(col, alpha.f = hatched_outline_alpha),
                    lwd = hatched_outline_lwd, lty = hatched_outline_lty
                )
                graphics::lines(
                    dc$time, dc$upper,
                    col = grDevices::adjustcolor(col, alpha.f = hatched_outline_alpha),
                    lwd = hatched_outline_lwd, lty = hatched_outline_lty
                )
            }
        }
        if (spread != "none" && envelope %in% c("lines", "both")) {
            graphics::lines(
                dc$time, dc$lower,
                col = grDevices::adjustcolor(col, alpha.f = 0.70),
                lwd = envelope_lwd, lty = envelope_lty
            )
            graphics::lines(
                dc$time, dc$upper,
                col = grDevices::adjustcolor(col, alpha.f = 0.70),
                lwd = envelope_lwd, lty = envelope_lty
            )
        }

        if (show_expectation) {
            graphics::lines(
                dc$time, dc$mean_expectation,
                col = col, lwd = expectation_lwd, lty = 1
            )
        }
    }

    if (isTRUE(draw_x_axis)) {
        .draw_time_axis(
            xlim,
            time_orientation = time_orientation,
            present_time = present_time,
            cex.axis = cex.axis
        )
    }
    graphics::axis(2, cex.axis = cex.axis)

    show_spread_legend <- spread != "none" && envelope != "none" &&
        !is.null(spread_legend_label)
    if (!is.null(legend_pos) && !identical(legend_pos, FALSE)) {
        legend_labels <- c(
            clade_names,
            if (show_spread_legend) spread_legend_label
        )
        legend_colors <- c(
            unname(cols[clade_names]),
            if (show_spread_legend) {
                grDevices::adjustcolor("grey40", alpha.f = envelope_alpha)
            }
        )
        graphics::legend(
            legend_pos,
            legend = legend_labels,
            col = legend_colors,
            lty = c(rep(1, length(clade_names)), if (show_spread_legend) NA),
            lwd = c(
                rep(expectation_lwd, length(clade_names)),
                if (show_spread_legend) NA
            ),
            pch = c(rep(NA, length(clade_names)), if (show_spread_legend) 15),
            bty = "n",
            cex = cex.legend
        )
    }

    invisible(list(
        colors = cols,
        xlim = xlim,
        ylim = ylim,
        bands = bands,
        envelope = spread,
        legend_pos = legend_pos,
        show_expectation = show_expectation,
        envelope_label = spread_legend_label,
        data = d
    ))
}

## Plot the reconstructed empirical variance
##
## The expectation is always shown. The median and an interquantile interval are
## shown only when the required variance quantiles are present.
##
## @param x A reconstruction object containing `empirical_variance`.
## @param clades Optional clade names to display.
## @param colors,palette Clade colours or an HCL palette name.
## @param interval How to draw an available quantile interval.
## @param interval_probs Two available probabilities delimiting the interval.
## @param interval_alpha Transparency of a filled interval.
## @param halo_style Whether the filled interval is solid or hatched.
## @param hatch_density,hatch_angle,hatch_lwd Hatching controls.
## @param hatched_outline Whether to outline a hatched region.
## @param hatched_outline_lty,hatched_outline_lwd,hatched_outline_alpha Style of
##   the hatched-region outline.
## @param show_median Whether to draw an available median.
## @param time_start,xlim,ylim Displayed time and value ranges.
## @param xlab,ylab Axis labels.
## @param time_orientation Whether time is measured from the root or before present.
## @param present_time Present time used for the reverse time scale.
## @param draw_x_axis Whether to draw the horizontal axis.
## @param bands,n_bands,band_color Optional alternating time bands.
## @param legend_pos Clade legend position, or `FALSE`. The redundant `all`
##   legend is omitted when no clade is specified.
## @param summary_legend Whether to display a legend for expectation and median.
## @param summary_legend_pos Position of the summary legend.
## @param expectation_lwd Width of the expectation curve.
## @param median_lwd,median_lty Style of the median curve.
## @param cex.lab,cex.axis,cex.legend Text scaling parameters.
## @param ... Additional graphical arguments.
## @return An invisible list containing plotting metadata.
.plot_empirical_variance_panel <- function(
    x,
    clades = NULL,
    colors = NULL,
    palette = "Dark 3",
    interval = c("ribbon", "lines", "both", "none"),
    interval_probs = NULL,
    interval_alpha = 0.14,
    halo_style = c("solid", "hatched"),
    hatch_density = 18,
    hatch_angle = 45,
    hatch_lwd = 0.7,
    hatched_outline = TRUE,
    hatched_outline_lty = 2,
    hatched_outline_lwd = 0.8,
    hatched_outline_alpha = 0.75,
    show_median = TRUE,
    time_start = NULL,
    xlim = NULL,
    ylim = NULL,
    xlab = NULL,
    ylab = "Empirical variance",
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    draw_x_axis = TRUE,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    legend_pos = "topright",
    summary_legend = TRUE,
    summary_legend_pos = "topleft",
    show_expectation = TRUE,
    expectation_lwd = 2,
    median_lwd = 1.3,
    median_lty = 2,
    cex.lab = 1.35,
    cex.axis = 1.0,
    cex.legend = 1.0,
    ...
) {
    time_orientation <- match.arg(time_orientation)
    if (is.null(xlab)) {
        xlab <- if (time_orientation == "from_root") "Time from root" else "Time before present"
    }

    .plot_empirical_quantity(
        x = x,
        quantity = "empirical_variance",
        clades = clades,
        colors = colors,
        palette = palette,
        interval = interval,
        interval_probs = interval_probs,
        interval_alpha = interval_alpha,
        halo_style = halo_style,
        hatch_density = hatch_density,
        hatch_angle = hatch_angle,
        hatch_lwd = hatch_lwd,
        hatched_outline = hatched_outline,
        hatched_outline_lty = hatched_outline_lty,
        hatched_outline_lwd = hatched_outline_lwd,
        hatched_outline_alpha = hatched_outline_alpha,
        show_expectation = show_expectation,
        show_median = show_median,
        expectation_lwd = expectation_lwd,
        median_lwd = median_lwd,
        median_lty = median_lty,
        time_start = time_start,
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        time_orientation = time_orientation,
        present_time = present_time,
        draw_x_axis = draw_x_axis,
        bands = bands,
        n_bands = n_bands,
        band_color = band_color,
        legend_pos = legend_pos,
        summary_legend = summary_legend,
        summary_legend_pos = summary_legend_pos,
        cex.lab = cex.lab,
        cex.axis = cex.axis,
        cex.legend = cex.legend,
        ...
    )
}


.plot_disparity_level_panel <- function(
    x,
    clades = NULL,
    colors = NULL,
    palette = "Dark 3",
    time_start = NULL,
    xlim = NULL,
    xlab = NULL,
    ylab = "Disparity level",
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    draw_x_axis = TRUE,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    geological_periods = NULL,
    geological_present_time = NULL,
    geological_background_alpha = 0.12,
    legend_pos = "topright",
    cex.lab = 1.35,
    cex.axis = 1,
    cex.legend = 1,
    curve_lwd = 2,
    reference_lwd = 1,
    reference_lty = 2,
    ...
) {
    time_orientation <- match.arg(time_orientation)
    d <- .subset_reconstruction_plot_data(
        x, "empirical_variance", clades = clades
    )
    if (!"disparity_level" %in% names(d)) {
        stop(
            "The reconstruction has no Brownian-reference probabilities. ",
            "Recompute it with quantities including 'disparity_level'.",
            call. = FALSE
        )
    }
    clade_names <- unique(as.character(d$clade))
    cols <- make_clade_colors(clade_names, colors = colors, palette = palette)
    xlim <- .resolve_time_window(d$time, time_start = time_start, xlim = xlim)
    d <- .trim_plot_data_to_window(d, xlim)
    bands <- .resolve_plot_bands(xlim, bands, n_bands)
    if (is.null(xlab)) {
        xlab <- if (time_orientation == "from_root") {
            "Time from root"
        } else {
            "Time before present"
        }
    }

    disparity_ylim <- c(-0.05, 1.05)
    graphics::plot(
        xlim, disparity_ylim, type = "n", axes = FALSE, frame.plot = FALSE,
        xlab = if (isTRUE(draw_x_axis)) xlab else "",
        ylab = ylab,
        cex.lab = cex.lab, xaxs = "i", yaxs = "i", ...
    )
    .draw_geological_background(
        geological_periods,
        present_time = geological_present_time,
        alpha = geological_background_alpha
    )
    .draw_time_bands(bands, band_color)
    graphics::abline(
        h = 0.5, col = "grey45", lty = reference_lty,
        lwd = reference_lwd
    )
    for (nm in clade_names) {
        dc <- d[d$clade == nm, , drop = FALSE]
        dc <- dc[order(dc$time), , drop = FALSE]
        graphics::lines(
            dc$time, dc$disparity_level,
            col = unname(cols[[nm]]), lwd = curve_lwd
        )
    }
    probability_ticks <- c(0, 0.5, 1)
    graphics::axis(
        2,
        at = probability_ticks,
        labels = c("0", "0.5", "1"),
        cex.axis = cex.axis
    )
    if (isTRUE(draw_x_axis)) {
        .draw_time_axis(
            xlim,
            time_orientation = time_orientation,
            present_time = present_time,
            cex.axis = cex.axis
        )
    }
    if (!is.null(legend_pos) && !identical(legend_pos, FALSE)) {
        graphics::legend(
            legend_pos,
            legend = clade_names,
            col = unname(cols[clade_names]),
            lty = 1,
            lwd = curve_lwd,
            bty = "n",
            cex = cex.legend
        )
    }
    invisible(list(
        colors = cols, xlim = xlim, ylim = disparity_ylim, bands = bands,
        legend_pos = legend_pos, data = d
    ))
}


.plot_variance_reference_panels <- function(
    x,
    variance_args,
    reference_args,
    cex.lab,
    cex.axis,
    geological_periods = NULL,
    geological_present_time = NULL,
    geological_background_alpha = 0.12
) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({
        graphics::layout(matrix(1L))
        graphics::par(oldpar)
    }, add = TRUE)
    text_layout <- .axis_text_layout(cex.lab, cex.axis)
    geological_scale <- !is.null(geological_periods)
    geological_margin <- .geological_axis_margin(
        text_layout, cex.axis, geological_periods
    )
    if (geological_scale) {
        outer_margin <- oldpar$oma
        outer_margin[1L] <- max(outer_margin[1L], geological_margin)
        graphics::par(oma = outer_margin)
        layout_heights <- .layout_heights_with_fixed_margins(c(1, 1), c(0.8, 0))
    } else {
		layout_heights <- c(1, 1)
    }
    graphics::layout(
        matrix(c(1, 2), ncol = 1L), heights = layout_heights
    )
    for (args_name in c("variance_args", "reference_args")) {
        args <- get(args_name)
        args$geological_periods <- geological_periods
        args$geological_present_time <- geological_present_time
        args$geological_background_alpha <- geological_background_alpha
        assign(args_name, args)
    }

    graphics::par(
        mar = c(0, text_layout$margin, 0.8, 0.8),
        mgp = text_layout$mgp
    )
    variance_args$draw_x_axis <- FALSE
    variance <- do.call(.plot_empirical_variance_panel, variance_args)

    graphics::par(
        mar = c(
            if (geological_scale) 0 else text_layout$margin,
            text_layout$margin, 0, 0.8
        ),
        mgp = text_layout$mgp
    )
    reference_args$draw_x_axis <- !geological_scale
    reference <- do.call(.plot_disparity_level_panel, reference_args)
    if (geological_scale) {
        .draw_geological_time_axis(
            reference$xlim, geological_present_time, geological_periods,
            cex.axis = cex.axis, cex.lab = cex.lab
        )
    }
    invisible(list(
        variance = variance,
        disparity_level = reference,
        xlim = variance$xlim,
        geological_scale = geological_scale,
        geological_background_alpha = geological_background_alpha
    ))
}

.plot_mean_disparity_panels <- function(
    mean_fun,
    mean_args,
    disparity_args,
    cex.lab,
    cex.axis,
    geological_periods = NULL,
    geological_present_time = NULL,
    geological_background_alpha = 0.12
) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({
        graphics::layout(matrix(1L))
        graphics::par(oldpar)
    }, add = TRUE)
    text_layout <- .axis_text_layout(cex.lab, cex.axis)
    geological_scale <- !is.null(geological_periods)
    geological_margin <- .geological_axis_margin(
        text_layout, cex.axis, geological_periods
    )
    if (geological_scale) {
        outer_margin <- oldpar$oma
        outer_margin[1L] <- max(outer_margin[1L], geological_margin)
        graphics::par(oma = outer_margin)
    }
    layout_heights <- if (geological_scale) {
        .layout_heights_with_fixed_margins(c(1, 1), c(0.8, 0))
    } else {
        c(1, 1)
    }
    graphics::layout(
        matrix(c(1, 2), ncol = 1L), heights = layout_heights
    )
    for (args_name in c("mean_args", "disparity_args")) {
        args <- get(args_name)
        args$geological_periods <- geological_periods
        args$geological_present_time <- geological_present_time
        args$geological_background_alpha <- geological_background_alpha
        assign(args_name, args)
    }

    graphics::par(
        mar = c(0, text_layout$margin, 0.8, 0.8),
        mgp = text_layout$mgp
    )
    mean_args$draw_x_axis <- FALSE
    mean <- do.call(mean_fun, mean_args)

    graphics::par(
        mar = c(
            if (geological_scale) 0 else text_layout$margin,
            text_layout$margin, 0, 0.8
        ),
        mgp = text_layout$mgp
    )
    disparity_args$draw_x_axis <- !geological_scale
    disparity <- do.call(
        .plot_disparity_level_panel, disparity_args
    )
    if (geological_scale) {
        .draw_geological_time_axis(
            disparity$xlim, geological_present_time, geological_periods,
            cex.axis = cex.axis, cex.lab = cex.lab
        )
    }
    invisible(list(
        mean = mean,
        disparity_level = disparity,
        xlim = mean$xlim,
        geological_scale = geological_scale,
        geological_background_alpha = geological_background_alpha
    ))
}


## Plot empirical mean and variance in aligned panels
##
## @param x A reconstruction object containing both empirical quantities.
## @param clades Optional clade names to display.
## @param colors,palette Clade colours or an HCL palette name.
## @param mean_envelope_style,variance_envelope_style Interval display modes.
## @param mean_envelope_probs,variance_envelope_probs Requested interval
##   probabilities for the two quantities.
## @param halo_style Common solid or hatched style for both intervals.
## @param mean_halo_style,variance_halo_style Panel-specific styles. An
##   explicitly supplied value takes priority over `halo_style`.
## @param halo_alpha Common transparency for both intervals.
## @param mean_halo_alpha,variance_halo_alpha Panel-specific transparencies. An
##   explicitly supplied value takes priority over `halo_alpha`.
## @param hatch_density,hatch_angle,hatch_lwd Hatching controls.
## @param hatched_outline Whether to outline hatched regions.
## @param hatched_outline_lty,hatched_outline_lwd,hatched_outline_alpha Style of
##   hatched-region outlines.
## @param variance_display Central variance curves to draw.
## @param time_start,xlim,mean_ylim,variance_ylim Displayed ranges.
## @param time_orientation Whether time is measured from the root or before present.
## @param present_time Present time used for the reverse time scale.
## @param bands,n_bands,band_color Optional alternating time bands.
## @param legend_pos Clade legend position, or `FALSE`.
## @param cex.lab,cex.axis,cex.legend Text scaling parameters.
## @param panel_heights Relative heights of the mean and variance panels.
## @param ... Additional graphical arguments.
## @return An invisible list containing metadata for both panels.
.plot_empirical_moment_panels <- function(
    x,
    clades = NULL,
    colors = NULL,
    palette = "Dark 3",
    mean_display = c("expectation", "none"),
    mean_envelope = c("conditional_quantiles", "expected_sd", "none"),
    mean_envelope_legend_label = NULL,
    mean_envelope_style = c("ribbon", "lines", "both", "none"),
    variance_envelope_style = c("ribbon", "lines", "both", "none"),
    mean_envelope_probs = NULL,
    variance_envelope_probs = NULL,
    mean_halo_style = c("solid", "hatched"),
    variance_halo_style = c("solid", "hatched"),
    halo_style = NULL,
    mean_halo_alpha = 0.14,
    variance_halo_alpha = 0.14,
    halo_alpha = NULL,
    hatch_density = 18,
    hatch_angle = 45,
    hatch_lwd = 0.7,
    hatched_outline = TRUE,
    hatched_outline_lty = 2,
    hatched_outline_lwd = 0.8,
    hatched_outline_alpha = 0.75,
    variance_display = c("both", "expectation", "median", "none"),
    mean_expectation_lwd = 2,
    mean_envelope_lwd = 0.8,
    mean_envelope_lty = 3,
    variance_expectation_lwd = 2,
    variance_median_lwd = 1.3,
    variance_median_lty = 2,
    time_start = NULL,
    xlim = NULL,
    mean_ylim = NULL,
    variance_ylim = NULL,
    mean_ylab = "Empirical mean",
    variance_ylab = "Empirical variance",
    disparity_level_ylab = "Disparity level",
    disparity_level_lwd = 2,
    disparity_level_reference_lwd = 1,
    disparity_level_reference_lty = 2,
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    legend_pos = "topright",
    cex.lab = 1.35,
    cex.axis = 1.0,
    cex.legend = 1.0,
    panel_heights = c(1, 1),
    show_disparity = FALSE,
    geological_periods = NULL,
    geological_present_time = NULL,
    geological_background_alpha = 0.12,
    ...
) {
    mean_halo_style_missing <- missing(mean_halo_style)
    variance_halo_style_missing <- missing(variance_halo_style)
    mean_halo_alpha_missing <- missing(mean_halo_alpha)
    variance_halo_alpha_missing <- missing(variance_halo_alpha)

    .validate_reconstruction_plot_data(x)
    time_orientation <- match.arg(time_orientation)
    mean_display <- match.arg(mean_display)
    mean_envelope <- match.arg(mean_envelope)
    variance_display <- match.arg(variance_display)
    mean_envelope_style <- match.arg(mean_envelope_style)
    variance_envelope_style <- match.arg(variance_envelope_style)
    if (!is.null(halo_style)) {
        if (mean_halo_style_missing) {
            mean_halo_style <- halo_style
        }
        if (variance_halo_style_missing) {
            variance_halo_style <- halo_style
        }
    }
    if (!is.null(halo_alpha)) {
        if (mean_halo_alpha_missing) {
            mean_halo_alpha <- halo_alpha
        }
        if (variance_halo_alpha_missing) {
            variance_halo_alpha <- halo_alpha
        }
    }
    mean_halo_style <- match.arg(mean_halo_style)
    variance_halo_style <- match.arg(variance_halo_style)
    if (!is.numeric(mean_halo_alpha) || length(mean_halo_alpha) != 1L ||
        !is.finite(mean_halo_alpha) || mean_halo_alpha < 0 || mean_halo_alpha > 1) {
        stop("mean_halo_alpha must be a finite scalar in [0, 1].")
    }
    if (!is.numeric(variance_halo_alpha) || length(variance_halo_alpha) != 1L ||
        !is.finite(variance_halo_alpha) || variance_halo_alpha < 0 ||
        variance_halo_alpha > 1) {
        stop("variance_halo_alpha must be a finite scalar in [0, 1].")
    }

    quantities <- unique(as.character(x$quantity))
    has_mean <- "empirical_mean" %in% quantities
    has_variance <- "empirical_variance" %in% quantities

    if (!has_mean || !has_variance) {
        stop("x must contain both empirical_mean and empirical_variance.")
    }

    if (is.null(clades)) {
        clades <- unique(as.character(x$clade))
    }
    cols <- make_clade_colors(clades, colors = colors, palette = palette)

    xlim <- .resolve_time_window(x$time, time_start = time_start, xlim = xlim)
    bands <- .resolve_plot_bands(xlim, bands, n_bands)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({
        graphics::layout(matrix(1L))
        graphics::par(oldpar)
    }, add = TRUE)

    if (!is.logical(show_disparity) ||
        length(show_disparity) != 1L ||
        is.na(show_disparity)) {
        stop("show_disparity must be TRUE or FALSE.")
    }
    if (show_disparity &&
        !"disparity_level" %in% names(x)) {
        stop(
            "show_disparity requires a reconstruction computed with ",
            "quantities including 'disparity_level'."
        )
    }
    layout_heights <- if (show_disparity) {
        c(panel_heights, 1)
    } else {
        panel_heights
    }
    text_layout <- .axis_text_layout(cex.lab, cex.axis)
    geological_scale <- !is.null(geological_periods)
    geological_margin <- .geological_axis_margin(
        text_layout, cex.axis, geological_periods
    )
    if (geological_scale) {
        outer_margin <- oldpar$oma
        outer_margin[1L] <- max(outer_margin[1L], geological_margin)
        graphics::par(oma = outer_margin)
    }
    graphics::layout(
        matrix(seq_along(layout_heights), ncol = 1L),
        heights = layout_heights
    )

    graphics::par(
        mar = c(0, text_layout$margin, 1, 0.8),
        mgp = text_layout$mgp
    )
    mean_expectation <- mean_display == "expectation"
    variance_expectation <- variance_display %in% c("both", "expectation")
    variance_median <- variance_display %in% c("both", "median")

    if (mean_envelope == "expected_sd") {
        mean_info <- .plot_empirical_spread_panel(
            x,
            clades = clades, colors = cols,
            envelope = mean_envelope_style, spread = mean_envelope,
            spread_legend_label = mean_envelope_legend_label,
            envelope_alpha = mean_halo_alpha,
            envelope_lwd = mean_envelope_lwd,
            envelope_lty = mean_envelope_lty,
            show_expectation = mean_expectation,
            expectation_lwd = mean_expectation_lwd,
            halo_style = mean_halo_style,
            hatch_density = hatch_density, hatch_angle = hatch_angle,
            hatch_lwd = hatch_lwd, hatched_outline = hatched_outline,
            hatched_outline_lty = hatched_outline_lty,
            hatched_outline_lwd = hatched_outline_lwd,
            hatched_outline_alpha = hatched_outline_alpha,
            time_start = time_start, xlim = xlim, ylim = mean_ylim,
            ylab = mean_ylab,
            time_orientation = time_orientation, present_time = present_time,
            draw_x_axis = FALSE, bands = bands, band_color = band_color,
            geological_periods = geological_periods,
            geological_present_time = geological_present_time,
            geological_background_alpha = geological_background_alpha,
            legend_pos = legend_pos, cex.lab = cex.lab,
            cex.axis = cex.axis, cex.legend = cex.legend, ...
        )
    } else {
        mean_info <- .plot_empirical_mean_panel(
            x,
            clades = clades,
            colors = cols,
            interval = if (mean_envelope == "conditional_quantiles") {
                mean_envelope_style
            } else {
                "none"
            },
            interval_probs = mean_envelope_probs,
            interval_alpha = mean_halo_alpha,
            show_expectation = mean_expectation,
            expectation_lwd = mean_expectation_lwd,
            halo_style = mean_halo_style,
            hatch_density = hatch_density,
            hatch_angle = hatch_angle,
            hatch_lwd = hatch_lwd,
            hatched_outline = hatched_outline,
            hatched_outline_lty = hatched_outline_lty,
            hatched_outline_lwd = hatched_outline_lwd,
            hatched_outline_alpha = hatched_outline_alpha,
            time_start = time_start,
            xlim = xlim,
            ylim = mean_ylim,
            ylab = mean_ylab,
            time_orientation = time_orientation,
            present_time = present_time,
            draw_x_axis = FALSE,
            bands = bands,
            band_color = band_color,
            geological_periods = geological_periods,
            geological_present_time = geological_present_time,
            geological_background_alpha = geological_background_alpha,
            legend_pos = legend_pos,
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            cex.legend = cex.legend,
            ...
        )
    }

    graphics::par(
        mar = c(
            if (show_disparity) {
                0
            } else if (geological_scale) {
                0
            } else {
                text_layout$margin
            },
            text_layout$margin, 0, 0.8
        ),
        mgp = text_layout$mgp
    )
    variance_info <- .plot_empirical_variance_panel(
        x,
        clades = clades,
        colors = cols,
        interval = variance_envelope_style,
        interval_probs = variance_envelope_probs,
        interval_alpha = variance_halo_alpha,
        halo_style = variance_halo_style,
        hatch_density = hatch_density,
        hatch_angle = hatch_angle,
        hatch_lwd = hatch_lwd,
        hatched_outline = hatched_outline,
        hatched_outline_lty = hatched_outline_lty,
        hatched_outline_lwd = hatched_outline_lwd,
        hatched_outline_alpha = hatched_outline_alpha,
        show_expectation = variance_expectation,
        show_median = variance_median,
        expectation_lwd = variance_expectation_lwd,
        median_lwd = variance_median_lwd,
        median_lty = variance_median_lty,
        time_start = time_start,
        xlim = xlim,
        ylim = variance_ylim,
        ylab = variance_ylab,
        time_orientation = time_orientation,
        present_time = present_time,
        draw_x_axis = !show_disparity && !geological_scale,
        bands = bands,
        band_color = band_color,
        geological_periods = geological_periods,
        geological_present_time = geological_present_time,
        geological_background_alpha = geological_background_alpha,
        legend_pos = FALSE,
        summary_legend = variance_display != "none",
        cex.lab = cex.lab,
        cex.axis = cex.axis,
        cex.legend = cex.legend,
        ...
    )
    if (!show_disparity && geological_scale) {
        .draw_geological_time_axis(
            variance_info$xlim, geological_present_time, geological_periods,
            cex.axis = cex.axis, cex.lab = cex.lab
        )
    }
    reference_info <- NULL
    if (show_disparity) {
        graphics::par(
            mar = c(
                if (geological_scale) 0 else text_layout$margin,
                text_layout$margin, 0, 0.8
            ),
            mgp = text_layout$mgp
        )
        reference_info <- .plot_disparity_level_panel(
            x, clades = clades, colors = cols,
            ylab = disparity_level_ylab,
            curve_lwd = disparity_level_lwd,
            reference_lwd = disparity_level_reference_lwd,
            reference_lty = disparity_level_reference_lty,
            time_start = time_start, xlim = xlim,
            time_orientation = time_orientation, present_time = present_time,
            draw_x_axis = !geological_scale,
            bands = bands, band_color = band_color,
            geological_periods = geological_periods,
            geological_present_time = geological_present_time,
            geological_background_alpha = geological_background_alpha,
            legend_pos = FALSE,
            cex.lab = cex.lab, cex.axis = cex.axis,
            cex.legend = cex.legend, ...
        )
        if (geological_scale) {
            .draw_geological_time_axis(
                reference_info$xlim, geological_present_time,
                geological_periods,
                cex.axis = cex.axis, cex.lab = cex.lab
            )
        }
    }

    invisible(list(
        colors = cols,
        xlim = xlim,
        bands = bands,
        halo_style = c(mean = mean_halo_style, variance = variance_halo_style),
        halo_alpha = c(mean = mean_halo_alpha, variance = variance_halo_alpha),
        geological_scale = geological_scale,
        geological_background_alpha = geological_background_alpha,
        mean = mean_info,
        variance = variance_info,
        disparity_level = reference_info
    ))
}


## =============================================================================
## Tree-colouring helpers
## =============================================================================

.compute_descendant_tips_for_plot <- function(tree) {
    n_tip <- length(tree$tip.label)
    n_nodes <- max(tree$edge)
    children <- vector("list", n_nodes)

    for (e in seq_len(nrow(tree$edge))) {
        parent <- tree$edge[e, 1]
        child <- tree$edge[e, 2]
        children[[parent]] <- c(children[[parent]], child)
    }

    descendants <- vector("list", n_nodes)

    visit <- function(node) {
        if (!is.null(descendants[[node]])) {
            return(descendants[[node]])
        }

        if (node <= n_tip) {
            descendants[[node]] <<- node
            return(node)
        }

        out <- integer(0)
        for (child in children[[node]]) {
            out <- c(out, visit(child))
        }
        descendants[[node]] <<- unique(out)
        descendants[[node]]
    }

    root <- setdiff(tree$edge[, 1], tree$edge[, 2])
    if (length(root) != 1L) {
        stop("Could not identify a unique root.")
    }
    visit(root[[1L]])
    descendants
}


.get_clade_mrca_for_plot <- function(tree, tip_indices) {
    if (length(tip_indices) == 0L) {
        return(NA_integer_)
    }
    if (length(tip_indices) == 1L) {
        return(tip_indices[[1L]])
    }

    get_mrca <- .get_mrca_function(tree)
    mrca <- tip_indices[[1L]]
    for (tip in tip_indices[-1L]) {
        mrca <- get_mrca(mrca, tip)
    }
    mrca
}


.edge_clade_membership_for_plot <- function(tree, clades,
                                            include_stem = FALSE) {
    clades <- .normalize_clades(tree, clades)
    membership <- matrix(
        FALSE,
        nrow = nrow(tree$edge),
        ncol = length(clades),
        dimnames = list(NULL, names(clades))
    )

    descendants <- .compute_descendant_tips_for_plot(tree)

    for (j in seq_along(clades)) {
        idx <- match(clades[[j]], tree$tip.label)
        missing <- clades[[j]][is.na(idx)]
        if (length(missing) > 0L) {
            warning(
                "Tips absent from tree for clade ", names(clades)[j], ": ",
                paste(missing, collapse = ", ")
            )
        }
        idx <- idx[!is.na(idx)]
        if (length(idx) == 0L) {
            next
        }

        mrca <- .get_clade_mrca_for_plot(tree, idx)

        for (e in seq_len(nrow(tree$edge))) {
            child <- tree$edge[e, 2]
            desc <- descendants[[child]]
            inside <- length(desc) > 0L && all(desc %in% idx)

            if (!inside) {
                next
            }
            if (!include_stem && length(idx) > 1L && child == mrca) {
                next
            }

            membership[e, j] <- TRUE
        }
    }

    membership
}


.resolve_membership_colors_for_plot <- function(membership, clade_colors,
                                                 other_color = "grey80",
                                                 overlap_color = "black",
                                                 overlap_policy = c("mixed", "first", "last")) {
    overlap_policy <- match.arg(overlap_policy)
    if (is.null(dim(membership))) {
        membership <- matrix(membership, nrow = 1L)
    }

    out <- rep(other_color, nrow(membership))

    for (i in seq_len(nrow(membership))) {
        idx <- which(membership[i, ])
        if (length(idx) == 0L) {
            next
        }
        if (length(idx) == 1L) {
            out[i] <- clade_colors[idx]
        } else if (overlap_policy == "first") {
            out[i] <- clade_colors[idx[1L]]
        } else if (overlap_policy == "last") {
            out[i] <- clade_colors[idx[length(idx)]]
        } else {
            out[i] <- overlap_color
        }
    }

    out
}


.assign_tree_colors <- function(tree, clades, clade_colors,
                                include_stem = FALSE,
                                other_edge_color = "grey80",
                                other_tip_color = "grey50",
                                overlap_color = "black",
                                overlap_policy = c("mixed", "first", "last")) {
    overlap_policy <- match.arg(overlap_policy)
    clades <- .normalize_clades(tree, clades)

    edge_membership <- .edge_clade_membership_for_plot(
        tree, clades, include_stem = include_stem
    )
    edge_colors <- .resolve_membership_colors_for_plot(
        edge_membership,
        clade_colors,
        other_color = other_edge_color,
        overlap_color = overlap_color,
        overlap_policy = overlap_policy
    )

    tip_membership <- sapply(clades, function(z) tree$tip.label %in% z)
    if (is.null(dim(tip_membership))) {
        tip_membership <- matrix(tip_membership, ncol = 1L)
    }
    tip_colors <- .resolve_membership_colors_for_plot(
        tip_membership,
        clade_colors,
        other_color = other_tip_color,
        overlap_color = overlap_color,
        overlap_policy = overlap_policy
    )
    names(tip_colors) <- tree$tip.label

    list(edge = edge_colors, tip = tip_colors)
}


.get_last_phylo_plot <- function() {
    env <- get(".PlotPhyloEnv", envir = asNamespace("ape"))
    get("last_plot.phylo", envir = env)
}


.redraw_phylogram_edges <- function(tree, edge_colors, edge_lwd = 1) {
    last <- .get_last_phylo_plot()
    xx <- last$xx
    yy <- last$yy

    if (length(edge_colors) == 1L) {
        edge_colors <- rep(edge_colors, nrow(tree$edge))
    }
    if (length(edge_lwd) == 1L) {
        edge_lwd <- rep(edge_lwd, nrow(tree$edge))
    }

    for (e in seq_len(nrow(tree$edge))) {
        parent <- tree$edge[e, 1]
        child <- tree$edge[e, 2]

        graphics::segments(
            xx[parent], yy[child], xx[child], yy[child],
            col = edge_colors[e], lwd = edge_lwd[e], lend = 1
        )
        graphics::segments(
            xx[parent], yy[parent], xx[parent], yy[child],
            col = edge_colors[e], lwd = edge_lwd[e], lend = 1
        )
    }

    invisible(last)
}


.geological_ranks <- c(
    "stage", "sub-epoch", "epoch", "sub-system", "period", "era", "eon",
    "super-eonothem"
)


.validate_geological_options <- function(rank, axis_args) {
    if (!is.character(rank) || length(rank) == 0L || anyNA(rank) ||
        any(!rank %in% .geological_ranks) || anyDuplicated(rank)) {
        stop(
            "geological_rank must be one or more unique values chosen ",
            "from: ", paste(.geological_ranks, collapse = ", "), "."
        )
    }
    if (!is.list(axis_args) ||
        length(axis_args) > 0L &&
            (is.null(names(axis_args)) || anyNA(names(axis_args)) ||
             any(names(axis_args) == "") || anyDuplicated(names(axis_args)))) {
        stop("geological_axis_args must be a named list with unique names.")
    }
    protected <- intersect(
        names(axis_args),
        c(
            "side", "intervals", "height", "tick_at", "tick_labels", "title",
            "phylo", "cex.axis"
        )
    )
    if (length(protected) > 0L) {
        stop(
            "geological_axis_args cannot override layout-controlled argument",
            if (length(protected) > 1L) "s: " else ": ",
            paste(protected, collapse = ", "), "."
        )
    }
    invisible(TRUE)
}


.get_geological_periods <- function(xlim, present_time, rank = "period",
                                    axis_args = list()) {
    .validate_geological_options(rank, axis_args)
    displayed_range <- sort(present_time - xlim)
    scale_data <- new.env(parent = emptyenv())
    utils::data(
        "GTS2020", package = "palaeoverse", envir = scale_data
    )
    gts <- scale_data$GTS2020
    endpoint_tol <- sqrt(.Machine$double.eps) *
        max(1, diff(displayed_range))
    intervals <- lapply(rank, function(current_rank) {
        selected <- gts[
            gts$rank == current_rank &
                gts$max_ma > displayed_range[1L] &
                gts$min_ma < displayed_range[2L],
            ,
            drop = FALSE
        ]
        if (nrow(selected) == 0L) {
            warning(
                "No geological ", current_rank,
                " interval overlaps the displayed time range; that rank has ",
                "been omitted.",
                call. = FALSE
            )
            return(NULL)
        }

        # axis_geo() 1.5.0 filters intervals with strict inequalities. Clip
        # the visible ends infinitesimally inwards so a view contained within
        # one interval is still rendered.
        selected$min_ma <- pmax(
            selected$min_ma, displayed_range[1L] + endpoint_tol
        )
        selected$max_ma <- pmin(
            selected$max_ma, displayed_range[2L] - endpoint_tol
        )
        selected
    })
    keep <- !vapply(intervals, is.null, logical(1L))
    intervals <- intervals[keep]
    rank <- rank[keep]
    if (length(intervals) == 0L) {
        warning(
            "No requested geological rank overlaps the displayed time range; ",
            "the geological scale has been omitted.",
            call. = FALSE
        )
        return(NULL)
    }
    periods <- if (length(intervals) == 1L) intervals[[1L]] else intervals
    attr(periods, "geological_rank") <- rank
    attr(periods, "axis_args") <- axis_args
    periods
}


.draw_geological_background <- function(periods, present_time, alpha = 0.12) {
    if (is.null(periods) || alpha <= 0) {
        return(invisible(NULL))
    }
    # With several scale rows, use the first requested rank for
    # panel backgrounds rather than overlaying several incompatible colours.
    background_periods <- if (is.list(periods) &&
        !is.data.frame(periods)) periods[[1L]] else periods
    usr <- graphics::par("usr")
    graphics::rect(
        xleft = pmax(usr[1L], present_time - background_periods$max_ma),
        xright = pmin(usr[2L], present_time - background_periods$min_ma),
        ybottom = usr[3L],
        ytop = usr[4L],
        col = grDevices::adjustcolor(
            background_periods$colour, alpha.f = alpha
        ),
        border = NA
    )
    invisible(NULL)
}

.geological_strip_height <- function(periods, cex.axis, cex.lab) {
    axis_args <- attr(periods, "axis_args", exact = TRUE)
    label_sizes <- unlist(
        axis_args$lab_size, recursive = TRUE, use.names = FALSE
    )
    label_sizes <- label_sizes[
        is.numeric(label_sizes) & is.finite(label_sizes) & label_sizes > 0
    ]
    largest_cex <- max(c(cex.axis, cex.lab, label_sizes), na.rm = TRUE)

    # axis_geo() expects a fraction of the current plot height. Convert a
    # font-based physical height to that fraction so the strip remains the
    # same size when the device or panel layout changes.
    strip_inches <- 2 * graphics::strheight(
        "Mg", units = "inches", cex = largest_cex
    )
    strip_inches / graphics::par("pin")[2L]
}


.draw_geological_time_axis <- function(xlim, present_time, periods,
                                       cex.axis = 1, cex.lab = 1.35) {
    if (is.null(periods)) {
        return(invisible(FALSE))
    }
    displayed_range <- sort(present_time - xlim)
    tick_at <- pretty(displayed_range, n = 6L)
    tick_at <- tick_at[
        tick_at >= displayed_range[1L] & tick_at <= displayed_range[2L]
    ]
    tick_labels <- format(round(tick_at, 4L), trim = TRUE)
    text_layout <- .axis_text_layout(cex.lab, cex.axis)
    usr <- graphics::par("usr")
    n_scales <- if (is.list(periods) && !is.data.frame(periods)) {
        length(periods)
    } else {
        1L
    }
    strip_height <- .geological_strip_height(periods, cex.axis, cex.lab)
    strip_bottom <- usr[3L] - n_scales * strip_height * diff(usr[3:4])
    plot_periods <- if (n_scales == 1L) list(periods) else periods
    plot_periods <- lapply(plot_periods, function(intervals) {
        period_max_ma <- intervals$max_ma
        intervals$max_ma <- present_time - intervals$min_ma
        intervals$min_ma <- present_time - period_max_ma
        intervals
    })
    if (n_scales == 1L) {
        plot_periods <- plot_periods[[1L]]
    }

    axis_args <- utils::modifyList(
        list(
        side = 1,
        intervals = plot_periods,
        height = rep(list(strip_height), n_scales),
        lab = TRUE,
        lab_size = cex.axis,
        abbr = FALSE,
        center_end_labels = TRUE,
        autofit = TRUE,
        # axis_geo() places enlarged tick labels too close to the strip. Let it
        # draw the geological intervals only; ticks and title are added below
        # with the same typography rules as the other panels.
        tick_at = numeric(),
        tick_labels = FALSE,
        title = NULL,
        phylo = FALSE,
        cex.axis = cex.axis
        ),
        attr(periods, "axis_args", exact = TRUE),
        keep.null = TRUE
    )
    do.call(palaeoverse::axis_geo, axis_args)
    graphics::axis(
        side = 1,
        pos = strip_bottom,
        at = present_time - tick_at,
        labels = tick_labels,
        cex.axis = cex.axis,
        mgp = text_layout$mgp,
        xpd = NA
    )
    strip_lines <- abs(
        graphics::grconvertY(usr[3L], from = "user", to = "lines") -
            graphics::grconvertY(strip_bottom, from = "user", to = "lines")
    )
    graphics::title(
        xlab = "Time before present (Ma)",
        line = strip_lines + text_layout$mgp[1L],
        cex.lab = cex.lab,
        xpd = NA
    )
    invisible(TRUE)
}


## =============================================================================
## Tree + trait + reconstruction figure
## =============================================================================

## Plot the tree, tip traits, and reconstructed empirical moments
##
## @param x A reconstruction object.
## @param tree The rooted phylogenetic tree used for reconstruction.
## @param trait Observed tip values.
## @param clades Clades to colour; by default, those stored in `x`.
## @param colors,palette Clade colours or an HCL palette name.
## @param include_stem Whether to colour the stem ancestral to a focal crown.
## @param other_edge_color,other_tip_color Colours outside focal clades.
## @param overlap_color,overlap_policy Display of overlapping memberships.
## @param mean_display Whether to draw the expected empirical mean.
## @param mean_envelope Envelope around the mean: conditional quantiles,
##   expected empirical standard deviation, or none.
## @param mean_envelope_legend_label Optional text identifying the empirical
##   spread in the clade legend.
## @param show_variance_panel Whether to include the variance panel.
## @param mean_envelope_style,variance_envelope_style Interval display modes.
## @param mean_envelope_probs,variance_envelope_probs Requested interval
##   probabilities.
## @param halo_style Common solid or hatched style for both moment panels.
## @param mean_halo_style,variance_halo_style Panel-specific styles. An
##   explicitly supplied value takes priority over `halo_style`.
## @param halo_alpha Common transparency for both moment panels.
## @param mean_halo_alpha,variance_halo_alpha Panel-specific transparencies. An
##   explicitly supplied value takes priority over `halo_alpha`.
## @param hatch_density,hatch_angle,hatch_lwd Hatching controls.
## @param hatched_outline Whether to outline hatched regions.
## @param hatched_outline_lty,hatched_outline_lwd,hatched_outline_alpha Style of
##   hatched-region outlines.
## @param variance_display Central variance curves to draw.
## @param time_start,xlim,mean_ylim,variance_ylim Displayed ranges.
## @param time_orientation Whether time is measured from the root or before present.
## @param present_time Present time used for the reverse time scale.
## @param bands,n_bands,band_color Optional alternating time bands.
## @param tree_lwd,trait_lwd Widths of tree edges and terminal trait bars.
## @param trait_axis Whether to draw the trait-value axis.
## @param trait_axis_label Label for the trait-value axis, or `NULL` for none.
## @param geological_scale Whether to add geological periods below the last
##   moment panel.
## @param legend_pos Clade legend position, or `FALSE`.
## @param cex.lab,cex.axis,cex.legend Text scaling parameters.
## @param tree_cex.lab,tree_cex.axis Text scaling parameters for the terminal
##   trait axis.
## @param tree_fraction Optional fraction of figure height assigned to the tree.
## @param tree_height,moment_height Relative heights when `tree_fraction` is
##   `NULL`.
## @param left_margin,right_margin Panel margins.
## @param ... Additional arguments forwarded to the moment-panel functions.
## @return An invisible list containing plotting and layout metadata.
.plot_reconstruction_with_tree <- function(
    x,
    tree,
    trait,
    panels = c("tree", "mean", "variance"),
    clades = NULL,
    colors = NULL,
    palette = "Dark 3",
    include_stem = FALSE,
    other_edge_color = "grey80",
    other_tip_color = "grey50",
    overlap_color = "black",
    overlap_policy = c("mixed", "first", "last"),
    mean_display = c("expectation", "none"),
    mean_envelope = c("conditional_quantiles", "expected_sd", "none"),
    mean_envelope_legend_label = NULL,
    mean_envelope_style = c("ribbon", "lines", "both", "none"),
    variance_envelope_style = c("ribbon", "lines", "both", "none"),
    mean_envelope_probs = NULL,
    variance_envelope_probs = NULL,
    mean_halo_style = NULL,
    variance_halo_style = c("solid", "hatched"),
    halo_style = NULL,
    mean_halo_alpha = 0.14,
    variance_halo_alpha = 0.14,
    halo_alpha = NULL,
    hatch_density = 18,
    hatch_angle = 45,
    hatch_lwd = 0.7,
    hatched_outline = TRUE,
    hatched_outline_lty = 2,
    hatched_outline_lwd = 0.8,
    hatched_outline_alpha = 0.75,
    variance_display = c("both", "expectation", "median", "none"),
    mean_expectation_lwd = 2,
    mean_envelope_lwd = 0.8,
    mean_envelope_lty = 3,
    variance_expectation_lwd = 2,
    variance_median_lwd = 1.3,
    variance_median_lty = 2,
    time_start = NULL,
    xlim = NULL,
    mean_ylim = NULL,
    variance_ylim = NULL,
    mean_ylab = "Empirical mean",
    variance_ylab = "Empirical variance",
    disparity_level_ylab = "Disparity level",
    disparity_level_lwd = 2,
    disparity_level_reference_lwd = 1,
    disparity_level_reference_lty = 2,
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    tree_lwd = 1,
    trait_lwd = 1,
    trait_axis = TRUE,
    trait_axis_label = "Trait value",
    geological_scale = FALSE,
    geological_rank = "period",
    geological_axis_args = list(),
    geological_background_alpha = 0.12,
    legend_pos = "topright",
    cex.lab = 1.35,
    cex.axis = 1.0,
    cex.legend = 1.0,
    tree_cex.lab = cex.lab,
    tree_cex.axis = cex.axis,
    show_disparity = "disparity_level" %in% panels,
    tree_fraction = NULL,
    tree_height = 1.2,
    moment_height = 1.1,
    left_margin = 4.2,
    right_margin = 0.6,
    ...
) {
    mean_halo_style_missing <- missing(mean_halo_style)
    variance_halo_style_missing <- missing(variance_halo_style)
    mean_halo_alpha_missing <- missing(mean_halo_alpha)
    variance_halo_alpha_missing <- missing(variance_halo_alpha)

    ## Additional graphical arguments in `...` are forwarded only to the
    ## moment-panel plotting functions. Function-specific arguments above are
    ## intercepted here and are therefore never passed accidentally to plot().
    .validate_reconstruction_plot_data(x)
    .validate_tree(tree)
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required to plot the phylogenetic tree.")
    }
    trait <- .align_trait_to_tree(tree, trait, "trait")

    panels <- unique(as.character(panels))
    if (length(panels) == 0L || anyNA(panels) ||
        any(!panels %in% c("tree", "mean", "variance", "disparity_level")) ||
        !"tree" %in% panels) {
        stop(
            "Tree layouts require panels containing 'tree' and optionally ",
            "'mean', 'variance', or 'disparity_level'."
        )
    }

    overlap_policy <- match.arg(overlap_policy)
    mean_display <- match.arg(mean_display)
    mean_envelope <- match.arg(mean_envelope)
    variance_display <- match.arg(variance_display)
    time_orientation <- match.arg(time_orientation)
    mean_envelope_style <- match.arg(mean_envelope_style)
    variance_envelope_style <- match.arg(variance_envelope_style)
    if (!is.logical(geological_scale) || length(geological_scale) != 1L ||
        is.na(geological_scale)) {
        stop("geological_scale must be TRUE or FALSE.")
    }
    .validate_geological_options(geological_rank, geological_axis_args)
    if (!is.logical(show_disparity) ||
        length(show_disparity) != 1L ||
        is.na(show_disparity)) {
        stop("show_disparity must be TRUE or FALSE.")
    }
    if (geological_scale && time_orientation != "before_present") {
        stop("geological_scale requires time_orientation = 'before_present'.")
    }
    if (geological_scale &&
        !requireNamespace("palaeoverse", quietly = TRUE)) {
        warning(
            "geological_scale was requested, but package 'palaeoverse' is ",
            "not installed; the geological scale has been omitted. Install ",
            "it with install.packages('palaeoverse').",
            call. = FALSE
        )
        geological_scale <- FALSE
    }
    if (!is.numeric(geological_background_alpha) ||
        length(geological_background_alpha) != 1L ||
        !is.finite(geological_background_alpha) ||
        geological_background_alpha < 0 ||
        geological_background_alpha > 1) {
        stop("geological_background_alpha must be a finite scalar in [0, 1].")
    }
    if (!is.null(halo_style)) {
        if (mean_halo_style_missing) {
            mean_halo_style <- halo_style
        }
        if (variance_halo_style_missing) {
            variance_halo_style <- halo_style
        }
    }
    if (is.null(mean_halo_style)) {
        mean_halo_style <- if (mean_envelope == "expected_sd") {
            "hatched"
        } else {
            "solid"
        }
    } else {
        mean_halo_style <- match.arg(mean_halo_style, c("solid", "hatched"))
    }
    variance_halo_style <- match.arg(variance_halo_style)
    if (!is.null(halo_alpha)) {
        if (mean_halo_alpha_missing) {
            mean_halo_alpha <- halo_alpha
        }
        if (variance_halo_alpha_missing) {
            variance_halo_alpha <- halo_alpha
        }
    }
    if (!is.numeric(mean_halo_alpha) || length(mean_halo_alpha) != 1L ||
        !is.finite(mean_halo_alpha) || mean_halo_alpha < 0 || mean_halo_alpha > 1) {
        stop("mean_halo_alpha must be a finite scalar in [0, 1].")
    }
    if (!is.numeric(variance_halo_alpha) || length(variance_halo_alpha) != 1L ||
        !is.finite(variance_halo_alpha) || variance_halo_alpha < 0 ||
        variance_halo_alpha > 1) {
        stop("variance_halo_alpha must be a finite scalar in [0, 1].")
    }
    if (!is.null(trait_axis_label) &&
        (!is.character(trait_axis_label) || length(trait_axis_label) != 1L ||
         is.na(trait_axis_label))) {
        stop("trait_axis_label must be NULL or a single non-missing character string.")
    }

    quantities <- unique(as.character(x$quantity))
    has_mean <- "mean" %in% panels && "empirical_mean" %in% quantities
    has_variance <- "variance" %in% panels &&
        "empirical_variance" %in% quantities

    if (length(quantities) == 0L) {
        stop("x contains no plottable empirical quantity.")
    }
    if (has_mean && mean_envelope == "expected_sd" &&
        !all(c("empirical_mean", "empirical_variance") %in% quantities)) {
        stop(
            "The selected mean_envelope requires both empirical_mean and ",
            "empirical_variance in x$data."
        )
    }
    plot_variance <- has_variance
    has_reference <- "disparity_level" %in% panels
    if (has_reference && !"disparity_level" %in% names(x)) {
        stop(
            "The requested disparity_level panel requires a reconstruction ",
            "computed with quantities including 'disparity_level'."
        )
    }
    if ("mean" %in% panels && !has_mean) {
        stop("The requested mean panel is absent from x.")
    }
    if ("variance" %in% panels && !has_variance) {
        stop("The requested variance panel is absent from x.")
    }

    if (is.null(clades)) {
        stop("clades must be supplied for a tree plot.")
    }
    clades <- .normalize_clades(tree, clades)

    plotted_clades <- unique(as.character(x$clade))
    missing_results <- setdiff(names(clades), plotted_clades)
    if (length(missing_results) > 0L) {
        warning(
            "No reconstruction rows are available for clades: ",
            paste(missing_results, collapse = ", "), "."
        )
        clades <- clades[setdiff(names(clades), missing_results)]
    }
    if (length(clades) == 0L) {
        stop("No clade remains after matching clades to reconstruction results.")
    }

    clade_names <- names(clades)
    implicit_all <- identical(clade_names, "all") &&
        setequal(clades[[1L]], tree$tip.label)
    if (implicit_all && is.null(colors)) {
        cols <- stats::setNames("black", "all")
    } else {
        cols <- make_clade_colors(clade_names, colors = colors, palette = palette)
    }
    clade_legend_pos <- if (implicit_all &&
        is.null(mean_envelope_legend_label)) FALSE else legend_pos

    xlim <- .resolve_time_window(x$time, time_start = time_start, xlim = xlim)
    if (is.null(present_time)) {
        present_time <- max(.compute_node_depths(tree))
    }
    geological_periods <- if (geological_scale) {
        .get_geological_periods(
            xlim, present_time, geological_rank, geological_axis_args
        )
    } else {
        NULL
    }
    if (geological_scale && is.null(geological_periods)) {
        geological_scale <- FALSE
    }
    bands <- .resolve_plot_bands(xlim, bands, n_bands)

    cw <- stats::reorder(tree, order = "cladewise")
    tree_colors <- .assign_tree_colors(
        cw,
        clades,
        clade_colors = cols,
        include_stem = include_stem,
        other_edge_color = other_edge_color,
        other_tip_color = other_tip_color,
        overlap_color = overlap_color,
        overlap_policy = overlap_policy
    )

    n_tip <- length(cw$tip.label)
    tree_ylim <- c(0.5, n_tip + 0.5)
    n_moment_panels <- as.integer(has_mean) + as.integer(plot_variance) +
        as.integer(has_reference)
    n_rows <- 1L + n_moment_panels

    if (n_moment_panels == 0L) {
        layout_heights <- tree_height
    } else if (is.null(tree_fraction)) {
        layout_heights <- c(tree_height, rep(moment_height, n_moment_panels))
    } else {
        if (!is.numeric(tree_fraction) || length(tree_fraction) != 1L ||
            !is.finite(tree_fraction) || tree_fraction <= 0 || tree_fraction >= 1) {
            stop("tree_fraction must be a finite scalar strictly between 0 and 1.")
        }
        layout_heights <- c(
            tree_fraction,
            rep((1 - tree_fraction) / n_moment_panels, n_moment_panels)
        )
    }

    layout_matrix <- matrix(0L, nrow = n_rows, ncol = 2L)
    layout_matrix[1L, ] <- c(1L, 2L)
    if (n_moment_panels > 0L) {
        panel_ids <- seq.int(3L, length.out = n_moment_panels)
        layout_matrix[seq.int(2L, n_rows), 1L] <- panel_ids
    }

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({
        graphics::layout(matrix(1L))
        graphics::par(oldpar)
    }, add = TRUE)

    moment_text_layout <- .axis_text_layout(cex.lab, cex.axis)
    trait_text_layout <- .axis_text_layout(tree_cex.lab, tree_cex.axis)
    trait_label_margin <- if (!isTRUE(trait_axis)) {
        0
    } else if (is.null(trait_axis_label)) {
        trait_text_layout$mgp[2L] + 0.85 * tree_cex.axis
    } else {
        trait_text_layout$margin
    }
    tree_bottom_margin <- if (n_moment_panels > 0L) 0 else trait_label_margin
    if (geological_scale && n_moment_panels > 0L) {
        geological_margin <- .geological_axis_margin(
            moment_text_layout, cex.axis, geological_periods
        )
        outer_margin <- oldpar$oma
        outer_margin[1L] <- max(outer_margin[1L], geological_margin)
        graphics::par(oma = outer_margin)
    }
    graphics::layout(
        layout_matrix,
        widths = c(5, 1),
        heights = layout_heights
    )

    ## --- Tree panel -------------------------------------------------------
    graphics::par(
        # Both panels in this layout row must use the same bottom margin to
        # preserve the vertical alignment of tips and terminal trait bars.
        mar = c(
            tree_bottom_margin, max(left_margin, moment_text_layout$margin),
            0.8, right_margin
        ),
        mgp = moment_text_layout$mgp, xaxs = "i", yaxs = "i"
    )

    ape::plot.phylo(
        cw,
        direction = "rightwards",
        show.tip.label = FALSE,
        edge.color = tree_colors$edge,
        edge.width = tree_lwd,
        x.lim = xlim,
        y.lim = tree_ylim,
        no.margin = FALSE
    )

    last <- .get_last_phylo_plot()
    tip_y <- last$yy[seq_along(cw$tip.label)]
    names(tip_y) <- cw$tip.label

    .draw_geological_background(
        geological_periods,
        present_time = present_time,
        alpha = geological_background_alpha
    )
    if (!is.null(bands)) {
        .draw_time_bands(bands, band_color)
    }
    if (!is.null(geological_periods) || !is.null(bands)) {
        .redraw_phylogram_edges(cw, tree_colors$edge, tree_lwd)
    }

    ## --- Terminal trait bars ---------------------------------------------
    graphics::par(
        mar = c(tree_bottom_margin, 0.2, 0.8, 1.0),
        mgp = trait_text_layout$mgp, xaxs = "i", yaxs = "i"
    )

    trait_ordered <- trait[cw$tip.label]
    xr <- range(trait_ordered, na.rm = TRUE)
    if (!all(is.finite(xr))) {
        stop("trait contains no finite value.")
    }
    if (xr[1L] == xr[2L]) {
        delta <- max(1, abs(xr[1L])) * 0.05
        xr <- xr + c(-delta, delta)
    }
    baseline <- if (xr[1L] < 0 && xr[2L] > 0) 0 else xr[1L]
    trait_xlim <- grDevices::extendrange(xr, f = 0.04)

    graphics::plot(
        trait_xlim, tree_ylim,
        type = "n", axes = FALSE,
        xlab = "", ylab = "",
        xaxs = "i", yaxs = "i"
    )
    graphics::segments(
        baseline,
        tip_y[cw$tip.label],
        trait_ordered,
        tip_y[cw$tip.label],
        col = unname(tree_colors$tip[cw$tip.label]),
        lwd = trait_lwd
    )
    graphics::abline(v = baseline, col = "grey60", lwd = 0.6)
    if (isTRUE(trait_axis)) {
        # With moment panels, the cells below the trait panel are deliberately
        # unused. Draw the annotations into that free area so the two top-row
        # panels can keep identical zero bottom margins without moving the
        # trait axis inside the plot.
        annotation_xpd <- if (n_moment_panels > 0L) NA else FALSE
        graphics::axis(
            1, cex.axis = tree_cex.axis, mgp = trait_text_layout$mgp,
            xpd = annotation_xpd
        )
        if (!is.null(trait_axis_label)) {
            # Use title(), exactly as plot.default() does for the mean and
            # variance panels. This avoids a device-dependent size difference
            # between mtext(cex = ...) and title(cex.lab = ...).
            graphics::title(
                xlab = trait_axis_label, line = trait_text_layout$mgp[1L],
                cex.lab = tree_cex.lab, adj = 1, xpd = annotation_xpd
            )
        }
    }

    ## --- Moment panels ----------------------------------------------------
    panel_index <- 0L
    mean_info <- variance_info <- reference_info <- NULL
    mean_expectation <- mean_display == "expectation"
    variance_expectation <- variance_display %in% c("both", "expectation")
    variance_median <- variance_display %in% c("both", "median")

    if (has_mean) {
        panel_index <- panel_index + 1L
        is_last <- panel_index == n_moment_panels
        graphics::par(
            mar = c(
                if (is_last && geological_scale) {
                    0
                } else if (is_last) {
                    moment_text_layout$margin
                } else 0,
                max(left_margin, moment_text_layout$margin), 0, right_margin
            ),
            mgp = moment_text_layout$mgp
        )

        if (mean_envelope == "expected_sd") {
            mean_info <- .plot_empirical_spread_panel(
                x,
                clades = clade_names,
                colors = cols,
                envelope = mean_envelope_style,
                spread = mean_envelope,
                spread_legend_label = mean_envelope_legend_label,
                envelope_alpha = mean_halo_alpha,
                envelope_lwd = mean_envelope_lwd,
                envelope_lty = mean_envelope_lty,
                show_expectation = mean_expectation,
                expectation_lwd = mean_expectation_lwd,
                halo_style = mean_halo_style,
                hatch_density = hatch_density,
                hatch_angle = hatch_angle,
                hatch_lwd = hatch_lwd,
                hatched_outline = hatched_outline,
                hatched_outline_lty = hatched_outline_lty,
                hatched_outline_lwd = hatched_outline_lwd,
                hatched_outline_alpha = hatched_outline_alpha,
                time_start = time_start,
                xlim = xlim,
                ylim = mean_ylim,
                ylab = mean_ylab,
                time_orientation = time_orientation,
                present_time = present_time,
                draw_x_axis = is_last && !geological_scale,
                bands = bands,
                band_color = band_color,
                geological_periods = geological_periods,
                geological_present_time = present_time,
                geological_background_alpha = geological_background_alpha,
                legend_pos = clade_legend_pos,
                cex.lab = cex.lab,
                cex.axis = cex.axis,
                cex.legend = cex.legend,
                ...
            )
        } else {
            mean_info <- .plot_empirical_mean_panel(
                x,
                clades = clade_names,
                colors = cols,
                interval = if (mean_envelope == "conditional_quantiles") {
                    mean_envelope_style
                } else {
                    "none"
                },
                interval_probs = mean_envelope_probs,
                interval_alpha = mean_halo_alpha,
                show_expectation = mean_expectation,
                expectation_lwd = mean_expectation_lwd,
                halo_style = mean_halo_style,
                hatch_density = hatch_density,
                hatch_angle = hatch_angle,
                hatch_lwd = hatch_lwd,
                hatched_outline = hatched_outline,
                hatched_outline_lty = hatched_outline_lty,
                hatched_outline_lwd = hatched_outline_lwd,
                hatched_outline_alpha = hatched_outline_alpha,
                time_start = time_start,
                xlim = xlim,
                ylim = mean_ylim,
                ylab = mean_ylab,
                time_orientation = time_orientation,
                present_time = present_time,
                draw_x_axis = is_last && !geological_scale,
                bands = bands,
                band_color = band_color,
                geological_periods = geological_periods,
                geological_present_time = present_time,
                geological_background_alpha = geological_background_alpha,
                legend_pos = clade_legend_pos,
                cex.lab = cex.lab,
                cex.axis = cex.axis,
                cex.legend = cex.legend,
                ...
            )
        }
        if (is_last && geological_scale) {
            .draw_geological_time_axis(
                xlim, present_time, geological_periods,
                cex.axis = cex.axis, cex.lab = cex.lab
            )
        }
    }

    if (plot_variance) {
        panel_index <- panel_index + 1L
        is_last <- panel_index == n_moment_panels
        graphics::par(
            mar = c(
                if (is_last && geological_scale) {
                    0
                } else if (is_last) {
                    moment_text_layout$margin
                } else 0,
                max(left_margin, moment_text_layout$margin), 0, right_margin
            ),
            mgp = moment_text_layout$mgp
        )

        variance_info <- .plot_empirical_variance_panel(
            x,
            clades = clade_names,
            colors = cols,
            interval = variance_envelope_style,
            interval_probs = variance_envelope_probs,
            interval_alpha = variance_halo_alpha,
            halo_style = variance_halo_style,
            hatch_density = hatch_density,
            hatch_angle = hatch_angle,
            hatch_lwd = hatch_lwd,
            hatched_outline = hatched_outline,
            hatched_outline_lty = hatched_outline_lty,
            hatched_outline_lwd = hatched_outline_lwd,
            hatched_outline_alpha = hatched_outline_alpha,
            show_expectation = variance_expectation,
            show_median = variance_median,
            expectation_lwd = variance_expectation_lwd,
            median_lwd = variance_median_lwd,
            median_lty = variance_median_lty,
            time_start = time_start,
            xlim = xlim,
            ylim = variance_ylim,
            ylab = variance_ylab,
            time_orientation = time_orientation,
            present_time = present_time,
            draw_x_axis = is_last && !geological_scale,
            bands = bands,
            band_color = band_color,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha,
            legend_pos = if (panel_index == 1L) clade_legend_pos else FALSE,
            summary_legend = variance_display != "none",
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            cex.legend = cex.legend,
            ...
        )
        if (is_last && geological_scale) {
            .draw_geological_time_axis(
                xlim, present_time, geological_periods,
                cex.axis = cex.axis, cex.lab = cex.lab
            )
        }
    }

    if (has_reference) {
        panel_index <- panel_index + 1L
        is_last <- panel_index == n_moment_panels
        graphics::par(
            mar = c(
                if (is_last && geological_scale) {
                    0
                } else if (is_last) {
                    moment_text_layout$margin
                } else 0,
                max(left_margin, moment_text_layout$margin), 0, right_margin
            ),
            mgp = moment_text_layout$mgp
        )
        reference_info <- .plot_disparity_level_panel(
            x,
            ylab = disparity_level_ylab,
            curve_lwd = disparity_level_lwd,
            reference_lwd = disparity_level_reference_lwd,
            reference_lty = disparity_level_reference_lty,
            clades = clade_names,
            colors = cols,
            time_start = time_start,
            xlim = xlim,
            time_orientation = time_orientation,
            present_time = present_time,
            draw_x_axis = is_last && !geological_scale,
            bands = bands,
            band_color = band_color,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha,
            legend_pos = if (panel_index == 1L) clade_legend_pos else FALSE,
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            cex.legend = cex.legend,
            ...
        )
        if (is_last && geological_scale) {
            .draw_geological_time_axis(
                xlim, present_time, geological_periods,
                cex.axis = cex.axis, cex.lab = cex.lab
            )
        }
    }

    invisible(list(
        colors = cols,
        edge_colors = tree_colors$edge,
        tip_colors = tree_colors$tip,
        tip_y = tip_y,
        xlim = xlim,
        tree_ylim = tree_ylim,
        trait_axis_label = trait_axis_label,
        moment_text_layout = moment_text_layout,
        trait_text_layout = trait_text_layout,
        clade_legend_pos = clade_legend_pos,
        layout_heights = layout_heights,
        tree_fraction = if (n_moment_panels == 0L) {
            1
        } else if (is.null(tree_fraction)) {
            layout_heights[1L] / sum(layout_heights)
        } else {
            tree_fraction
        },
        geological_scale = geological_scale,
        geological_background_alpha = geological_background_alpha,
        bands = bands,
        panels = panels,
        halo_style = c(mean = mean_halo_style, variance = variance_halo_style),
        halo_alpha = c(mean = mean_halo_alpha, variance = variance_halo_alpha),
        mean = mean_info,
        variance = variance_info,
        disparity_level = reference_info
    ))
}


## =============================================================================
## S3 plot method
## =============================================================================

#' Plot an empirical-moment reconstruction
#'
#' This is the single public plotting interface for reconstructions. It can
#' display the phylogenetic tree with terminal trait bars, the reconstructed
#' empirical mean, the reconstructed empirical variance, or any compatible
#' combination. The tree, aligned tip trait and clades are read from `x`.
#' Shared graphical settings are inherited by the moment panels; panel-specific
#' settings take precedence.
#'
#' @param x An object of class `empirical_moment_reconstruction`.
#' @param panels Character vector selecting `"tree"`, `"mean"`, `"variance"`,
#'   and/or `"disparity_level"`. Defaults to the available tree and moment panels;
#'   the optional disparity-level panel must be requested explicitly. Panels are
#'   drawn in canonical order.
#' @param colors,palette Clade colours or an HCL palette name. `colors` may be
#'   a vector named by clade; unspecified colours are taken from `palette`,
#'   whose default is `"Dark 3"`.
#' @param mean_display Whether to draw the expected empirical mean:
#'   `"expectation"` (default) or `"none"`.
#' @param mean_envelope Envelope around the mean: conditional quantiles,
#'   expected empirical standard deviation, or none. Accepted values are
#'   `"conditional_quantiles"` (default), `"expected_sd"`, and `"none"`.
#' @param mean_envelope_legend_label Optional text identifying the
#'   lineage-dispersion envelope in the clade legend. `NULL` (default) adds no
#'   envelope entry. Its position is controlled by `legend_pos`.
#' @param mean_envelope_style,variance_envelope_style Envelope display modes:
#'   `"ribbon"` (default), `"lines"`, `"both"`, or `"none"`.
#' @param mean_envelope_probs,variance_envelope_probs Requested interval
#'   probabilities as a numeric vector of length two, for example
#'   `c(0.025, 0.975)`. `NULL` selects the outermost available interval.
#' @param variance_display Central variance curves to draw: expectation,
#'   median, both (default), or neither.
#' @param common_style Named list of graphical settings inherited by both
#'   moment panels. Every supported name and its default are listed in Details.
#' @param mean_style,variance_style,disparity_level_style,tree_style Named lists of panel-specific
#'   graphical settings. Every supported name and its default are listed in
#'   Details. Unknown names cause an error.
#' @param time_start,xlim Displayed time range. `time_start` is a nonnegative
#'   lower truncation (default `0`); `xlim` is an increasing numeric vector of
#'   length two and defaults to the available time range.
#' @param time_orientation Whether time is measured `"from_root"` (default) or
#'   `"before_present"`.
#' @param present_time Present time used for the reverse time scale. `NULL`
#'   uses the greatest node depth in the tree.
#' @param bands,n_bands,band_color Optional alternating time bands. `bands` is
#'   an explicit set of boundaries; otherwise `n_bands` equally spaced bands
#'   are generated (`0` by default). `band_color` defaults to transparent grey.
#' @param legend_pos Clade legend position accepted by [graphics::legend()], or
#'   `FALSE`/`NULL` to suppress it. Default `"topright"`.
#' @param geological_scale Logical. If `TRUE`, add geological periods below the
#'   last time-based panel using [palaeoverse::axis_geo()]. This requires
#'   `time_orientation = "before_present"` and the suggested package
#'   \pkg{palaeoverse}; the tree panel need not be displayed. Default `FALSE`.
#' @param geological_rank One or more geological ranks displayed as stacked
#'   scale rows, in the supplied order. Accepted values are `"stage"`,
#'   `"sub-epoch"`, `"epoch"`, `"sub-system"`, `"period"` (default), `"era"`,
#'   `"eon"`, and `"super-eonothem"`.
#' @param geological_axis_args Named list of additional options passed to
#'   [palaeoverse::axis_geo()], for example `list(abbr = TRUE, rot = 45)`.
#'   Arguments controlling integration with the surrounding plot (`side`,
#'   `intervals`, `height`, `tick_at`, `tick_labels`, `title`, `phylo`, and
#'   `cex.axis`) cannot be overridden.
#' @param geological_background_alpha Opacity of the geological-period colours
#'   drawn behind every time-based panel when `geological_scale = TRUE`.
#'   Must lie in `[0, 1]`; `0` suppresses the coloured backgrounds. Default
#'   `0.12`.
#' @param ... Additional low-level graphical arguments forwarded to moment
#'   panels and ultimately to [graphics::plot.default()], such as `main`,
#'   `sub`, `las`, or `bty`. They do not style the tree panel.
#' @return Plotting and layout metadata, invisibly.
#' @details
#' The following sections form an exhaustive reference to the options accepted
#' by this method. Options passed in a `*_style` list are deliberately separated
#' from top-level options; unknown names in a style list cause an error.
#'
#' ## General options
#'
#' - `panels = NULL`: panels to display, chosen from `"tree"`, `"mean"`,
#'   `"variance"`, and `"disparity_level"`. `NULL` selects every available tree and
#'   moment panel, but not the optional disparity-level panel. Irrespective of the
#'   supplied order, panels are drawn in canonical
#'   tree--mean--variance--disparity_level order.
#' - `colors = NULL`: clade colours. A named vector matches colours to clade
#'   names; missing or unnamed colours are completed from `palette`.
#' - `palette = "Dark 3"`: HCL palette used to generate missing clade colours.
#' - `legend_pos = "topright"`: position of the clade legend, using any value
#'   accepted by [graphics::legend()]; `FALSE` or `NULL` suppresses it.
#' - `common_style = list()`: styles shared by the mean and variance panels.
#'   Its complete set of names is listed under `common_style` below.
#' - `mean_style = list()`, `variance_style = list()`,
#'   `disparity_level_style = list()`, `tree_style = list()`:
#'   styles local to one panel. Their complete sets of names are listed below.
#' - `...`: additional arguments for [graphics::plot.default()] in the moment
#'   panels, for example `main`, `sub`, `las`, or `bty`. These arguments do not
#'   style the tree panel.
#'
#' ## Time axis, bands, and geological scale
#'
#' - `time_start = NULL`: nonnegative lower time cutoff. `NULL` is interpreted
#'   as `0`.
#' - `xlim = NULL`: increasing two-element displayed time range. `NULL` uses
#'   the available range after applying `time_start`.
#' - `time_orientation = "from_root"`: axis direction; use `"before_present"`
#'   for ages measured backward from the present.
#' - `present_time = NULL`: present used by `"before_present"`. `NULL` uses the
#'   greatest node depth in the stored tree.
#' - `bands = NULL`: explicit numeric boundaries for alternating vertical time
#'   bands.
#' - `n_bands = 0`: number of equally spaced alternating bands generated when
#'   `bands` is `NULL`; zero draws no bands.
#' - `band_color = adjustcolor("grey", alpha.f = 0.10)`: band colour.
#' - `geological_scale = FALSE`: add a geological-period axis below the last
#'   time-based panel. This requires `time_orientation = "before_present"` and
#'   the suggested package \pkg{palaeoverse}, but does not require the tree
#'   panel.
#' - `geological_rank = "period"`: one or more stratigraphic ranks displayed
#'   as stacked rows in the supplied order; for example
#'   `c("stage", "period", "era")`. The first rank supplies the transparent
#'   geological backgrounds drawn behind the time-based panels.
#' - `geological_axis_args = list()`: additional named options for
#'   [palaeoverse::axis_geo()], such as `list(abbr = TRUE, rot = 45)`.
#' - `geological_background_alpha = 0.12`: opacity in `[0, 1]` of geological
#'   period colours behind all time-based panels; zero suppresses them.
#'
#' With `geological_scale = TRUE`, numerical ticks and the title
#' `"Time before present (Ma)"` are placed outside the geological-period row.
#' The range is derived from `xlim` and `present_time`; no geological age limit
#' is hard-coded. The terminal-trait inset is not shaded because its horizontal
#' axis represents the trait rather than time.
#'
#' ## Mean-panel options
#'
#' - `mean_display = "expectation"`: draw the expected empirical mean;
#'   `"none"` suppresses the central curve.
#' - `mean_envelope = "conditional_quantiles"`: envelope quantity. Other
#'   choices are `"expected_sd"` and `"none"`.
#' - `mean_envelope_legend_label = NULL`: add an `"expected_sd"` entry to the
#'   clade legend when a character string is supplied; `NULL` adds no entry.
#' - `mean_envelope_style = "ribbon"`: envelope rendering; `"lines"` draws
#'   boundaries only, `"both"` draws ribbon and boundaries, and `"none"`
#'   suppresses both.
#' - `mean_envelope_probs = NULL`: two interval probabilities, such as
#'   `c(0.025, 0.975)`; `NULL` selects the outermost available interval.
#' - `mean_style = list()`: panel styles listed under `mean_style` below.
#'
#' ## Variance-panel options
#'
#' - `variance_display = "both"`: central curves to draw. Other choices are
#'   `"expectation"`, `"median"`, and `"none"`.
#' - `variance_envelope_style = "ribbon"`: envelope rendering; `"lines"` draws
#'   boundaries only, `"both"` draws ribbon and boundaries, and `"none"`
#'   suppresses both.
#' - `variance_envelope_probs = NULL`: two interval probabilities, such as
#'   `c(0.025, 0.975)`; `NULL` selects the outermost available interval.
#' - `variance_style = list()`: panel styles listed under `variance_style`
#'   below.
#' - `panels = "disparity_level"`: display
#'   `P(reconstructed variance > Brownian variance)` on `[0, 1]`, with a
#'   horizontal reference at `0.5`. This is a descriptive probability, not a
#'   p-value, and requires a reconstruction computed with
#'   `quantities` including `"disparity_level"`.
#'
#' A requested median or quantile envelope is silently omitted when the
#' reconstruction does not contain the required quantiles. Halo options affect
#' ribbons only and therefore have no visible effect when the corresponding
#' envelope style is `"lines"` or `"none"`.
#'
#' ## `common_style`
#'
#' These settings are inherited by both moment panels. A value in `mean_style`
#' or `variance_style` takes precedence where the same name is supported.
#'
#' - `halo_style = "solid"`: filled-envelope style, `"solid"` or `"hatched"`.
#' - `halo_alpha = 0.14`: filled-envelope opacity in `[0, 1]`.
#' - `hatch_density = 18`, `hatch_angle = 45`, `hatch_lwd = 0.7`: density,
#'   angle in degrees, and line width of hatching.
#' - `hatched_outline = TRUE`: draw borders around hatched regions.
#' - `hatched_outline_lty = 2`, `hatched_outline_lwd = 0.8`,
#'   `hatched_outline_alpha = 0.75`: line type, width, and opacity of those
#'   borders.
#' - `expectation_lwd = 2`: expectation-curve width.
#' - `cex.lab = 1.35`, `cex.axis = 1`, `cex.legend = 1`: label, axis-text,
#'   and legend-text scaling. These settings also control the terminal-trait
#'   axis (`cex.lab` for its title and `cex.axis` for its tick labels). Axis
#'   spacing and plot margins expand automatically when either value grows.
#'
#' ## `mean_style`
#'
#' - `ylab = "Empirical mean"`: vertical-axis label.
#' - `ylim = NULL`: fixed two-element y range; `NULL` computes it from data.
#' - `expectation_lwd = 2`: expected-mean curve width.
#' - `envelope_lwd = 0.8`, `envelope_lty = 3`: boundary style for an
#'   `"expected_sd"` envelope.
#' - `ylim_expand = 0.06`: relative padding added to an automatically computed
#'   y range for an empirical-spread envelope.
#'
#' ## `variance_style`
#'
#' - `ylab = "Empirical variance"`: vertical-axis label.
#' - `ylim = NULL`: fixed two-element y range; `NULL` computes it from data.
#' - `expectation_lwd = 2`: expected-variance curve width.
#' - `median_lwd = 1.3`, `median_lty = 2`: median-curve width and line type.
#'
#' ## `disparity_level_style`
#'
#' - `ylab = "Disparity level"`: vertical-axis label.
#' - `curve_lwd = 2`: width of the probability curves.
#' - `reference_lwd = 1`, `reference_lty = 2`: width and line type of the
#'   horizontal reference at `0.5`.
#'
#' ## `tree_style`
#'
#' These options affect the tree panel and its terminal-trait inset only.
#'
#' - `tree_lwd = 1`, `trait_lwd = 1`: widths of tree edges and terminal trait
#'   bars.
#' - `trait_axis = TRUE`: show the terminal-trait axis.
#' - `trait_axis_label = "Trait value"`: its label; `NULL` suppresses the label.
#' - `cex.lab`, `cex.axis`: text scaling for the terminal-trait axis title and
#'   tick labels, respectively. By default, each inherits the corresponding
#'   value from `common_style`.
#' - `tree_fraction = NULL`: fraction of total figure height assigned to the
#'   tree. When specified it must lie strictly between 0 and 1; the moment
#'   panels divide the remainder equally.
#' - `tree_height = 1.2`, `moment_height = 1.1`: relative tree and per-moment
#'   panel heights used only when `tree_fraction` is `NULL`.
#' - `left_margin = 4.2`, `right_margin = 0.6`: tree-panel margins in text
#'   lines.
#' - `include_stem = FALSE`: whether the edge entering a selected clade's MRCA
#'   receives the clade colour.
#' - `other_edge_color = "grey80"`, `other_tip_color = "grey50"`: colours for
#'   edges and tips outside selected clades.
#' - `overlap_color = "black"`: colour used for membership in multiple clades
#'   when `overlap_policy = "mixed"`.
#' - `overlap_policy = "mixed"`: treatment of overlapping clades; `"mixed"`
#'   uses `overlap_color`, while `"first"` or `"last"` uses the corresponding
#'   clade colour.
#'
#' @examples
#' \dontrun{
#' # Give the tree exactly 40% of the total figure height.
#' plot(reconstruction, tree_style = list(tree_fraction = 0.4))
#'
#' # Relative sizing is available when no exact fraction is requested.
#' plot(
#'   reconstruction,
#'   tree_style = list(tree_height = 2, moment_height = 1)
#' )
#'
#' # Common settings are inherited and can be overridden per moment panel.
#' plot(
#'   reconstruction,
#'   common_style = list(halo_style = "hatched", halo_alpha = 0.2),
#'   mean_style = list(halo_alpha = 0.35),
#'   variance_style = list(median_lty = 3)
#' )
#'
#' # Stack stages, periods, and eras, and customise their labels. The first
#' # rank supplies the transparent geological backgrounds.
#' plot(
#'   reconstruction,
#'   time_orientation = "before_present",
#'   geological_scale = TRUE,
#'   geological_rank = c("stage", "period", "era"),
#'   geological_axis_args = list(abbr = TRUE, rot = 45)
#' )
#'
#' # Display the optional Brownian probability-of-superiority curve.
#' plot(reconstruction_with_reference, panels = c("variance", "disparity_level"))
#' }
#' @rdname plot.empirical_moment_reconstruction
#' @export
plot.empirical_moment_reconstruction <- function(
    x,
    panels = NULL,
    colors = NULL,
    palette = "Dark 3",
    mean_display = c("expectation", "none"),
    mean_envelope = c("conditional_quantiles", "expected_sd", "none"),
    mean_envelope_legend_label = NULL,
    mean_envelope_style = c("ribbon", "lines", "both", "none"),
    variance_envelope_style = c("ribbon", "lines", "both", "none"),
    mean_envelope_probs = NULL,
    variance_envelope_probs = NULL,
    variance_display = c("both", "expectation", "median", "none"),
    common_style = list(),
    mean_style = list(),
    variance_style = list(),
    disparity_level_style = list(),
    tree_style = list(),
    time_start = NULL,
    xlim = NULL,
    time_orientation = c("from_root", "before_present"),
    present_time = NULL,
    bands = NULL,
    n_bands = 0L,
    band_color = grDevices::adjustcolor("grey", alpha.f = 0.10),
    legend_pos = "topright",
    geological_scale = FALSE,
    geological_rank = "period",
    geological_axis_args = list(),
    geological_background_alpha = 0.12,
    ...
) {
    dots <- list(...)
    call_names <- names(as.list(sys.call())[-1L])
    biological_inputs <- intersect(
        c(names(dots), call_names),
        c("tree", "trait", "clades")
    )
    if (length(biological_inputs) > 0L) {
        stop(
            "tree, trait, and clades are stored in x and cannot be supplied ",
            "to plot(). Modify the corresponding components of x instead."
        )
    }

    .validate_reconstruction_object(x)
    tree <- x$tree
    trait <- x$trait
    clade_definitions <- x$clades
    clades <- names(clade_definitions)
    implicit_all <- identical(clades, "all") &&
        setequal(clade_definitions[[1L]], tree$tip.label)
    if (implicit_all && is.null(colors)) {
        colors <- stats::setNames("black", "all")
    }
    x <- x$data
    quantities <- unique(as.character(x$quantity))

    mean_display <- match.arg(mean_display)
    mean_envelope <- match.arg(mean_envelope)
    variance_display <- match.arg(variance_display)
    mean_envelope_style <- match.arg(mean_envelope_style)
    variance_envelope_style <- match.arg(variance_envelope_style)
    time_orientation <- match.arg(time_orientation)
    if (!is.null(mean_envelope_legend_label) &&
        (!is.character(mean_envelope_legend_label) ||
         length(mean_envelope_legend_label) != 1L ||
         is.na(mean_envelope_legend_label))) {
        stop(
            "mean_envelope_legend_label must be NULL or a single ",
            "character string."
        )
    }
    if (!is.logical(geological_scale) || length(geological_scale) != 1L ||
        is.na(geological_scale)) {
        stop("geological_scale must be TRUE or FALSE.")
    }
    .validate_geological_options(geological_rank, geological_axis_args)
    available <- c(
        "tree",
        if ("empirical_mean" %in% quantities) "mean",
        if ("empirical_variance" %in% quantities) "variance"
    )
    if (is.null(panels)) {
        panels <- available
    }
    panels <- unique(as.character(panels))
    if (length(panels) == 0L || anyNA(panels) ||
        any(!panels %in% c("tree", "mean", "variance", "disparity_level"))) {
        stop(
            "panels must select at least one of 'tree', 'mean', 'variance', ",
            "or 'disparity_level'."
        )
    }
    unavailable <- setdiff(panels, c(
        available,
        if ("disparity_level" %in% names(x)) "disparity_level"
    ))
    if (length(unavailable) > 0L) {
        if ("disparity_level" %in% unavailable) {
            stop(
                "The requested disparity_level panel requires a reconstruction ",
                "computed with quantities including 'disparity_level'."
            )
        }
        stop("Requested panel is unavailable: ", paste(unavailable, collapse = ", "))
    }
    panels <- intersect(c("tree", "mean", "variance", "disparity_level"), panels)
    show_disparity <- "disparity_level" %in% panels
    if (geological_scale && time_orientation != "before_present") {
        stop("geological_scale requires time_orientation = 'before_present'.")
    }
    if (!is.numeric(geological_background_alpha) ||
        length(geological_background_alpha) != 1L ||
        !is.finite(geological_background_alpha) ||
        geological_background_alpha < 0 ||
        geological_background_alpha > 1) {
        stop("geological_background_alpha must be a finite scalar in [0, 1].")
    }
    if (geological_scale &&
        !requireNamespace("palaeoverse", quietly = TRUE)) {
        warning(
            "geological_scale was requested, but package 'palaeoverse' is ",
            "not installed; the geological scale has been omitted. Install ",
            "it with install.packages('palaeoverse').",
            call. = FALSE
        )
        geological_scale <- FALSE
    }
    xlim <- .resolve_time_window(x$time, time_start = time_start, xlim = xlim)
    if (is.null(present_time)) {
        present_time <- max(.compute_node_depths(tree))
    }
    geological_periods <- if (geological_scale) {
        .get_geological_periods(
            xlim, present_time, geological_rank, geological_axis_args
        )
    } else {
        NULL
    }
    if (geological_scale && is.null(geological_periods)) {
        geological_scale <- FALSE
    }

    style_lists <- list(
        common_style = common_style, mean_style = mean_style,
        variance_style = variance_style,
        disparity_level_style = disparity_level_style,
        tree_style = tree_style
    )
    for (nm in names(style_lists)) {
        value <- style_lists[[nm]]
        if (!is.list(value) || is.null(names(value)) && length(value) > 0L ||
            any(names(value) == "")) {
            stop(nm, " must be a named list.")
        }
    }
    common_allowed <- c(
        "halo_style", "halo_alpha", "hatch_density", "hatch_angle",
        "hatch_lwd", "hatched_outline", "hatched_outline_lty",
        "hatched_outline_lwd", "hatched_outline_alpha", "expectation_lwd",
        "cex.lab", "cex.axis", "cex.legend"
    )
    mean_allowed <- c("halo_style", "halo_alpha", "ylim", "ylab", "expectation_lwd",
                      "envelope_lwd", "envelope_lty", "ylim_expand")
    variance_allowed <- c("halo_style", "halo_alpha", "ylim", "ylab",
                          "expectation_lwd", "median_lwd", "median_lty")
    disparity_level_allowed <- c(
        "ylab", "curve_lwd", "reference_lwd", "reference_lty"
    )
    tree_allowed <- c(
        "tree_lwd", "trait_lwd", "trait_axis", "trait_axis_label",
        "cex.lab", "cex.axis",
        "tree_fraction", "tree_height", "moment_height", "left_margin",
        "right_margin", "include_stem", "other_edge_color", "other_tip_color",
        "overlap_color", "overlap_policy"
    )
    allowed <- list(
        common_style = common_allowed, mean_style = mean_allowed,
        variance_style = variance_allowed,
        disparity_level_style = disparity_level_allowed,
        tree_style = tree_allowed
    )
    for (nm in names(style_lists)) {
        unknown <- setdiff(names(style_lists[[nm]]), allowed[[nm]])
        if (length(unknown) > 0L) {
            stop("Unknown ", nm, " setting: ", paste(unknown, collapse = ", "))
        }
    }

    common_defaults <- list(
        halo_style = "solid", halo_alpha = 0.14, hatch_density = 18,
        hatch_angle = 45, hatch_lwd = 0.7, hatched_outline = TRUE,
        hatched_outline_lty = 2, hatched_outline_lwd = 0.8,
        hatched_outline_alpha = 0.75, expectation_lwd = 2,
        cex.lab = 1.35, cex.axis = 1, cex.legend = 1
    )
    common <- utils::modifyList(common_defaults, common_style)
    mean <- utils::modifyList(
        c(common, list(ylim = NULL, ylab = "Empirical mean", envelope_lwd = 0.8,
                       envelope_lty = 3, ylim_expand = 0.06)),
        mean_style
    )
    variance <- utils::modifyList(
        c(common, list(ylim = NULL, ylab = "Empirical variance",
                       median_lwd = 1.3, median_lty = 2)),
        variance_style
    )
    disparity_level <- utils::modifyList(
        list(
            ylab = "Disparity level", curve_lwd = 2,
            reference_lwd = 1, reference_lty = 2
        ),
        disparity_level_style
    )
    tree_defaults <- list(
        tree_lwd = 1, trait_lwd = 1, trait_axis = TRUE,
        trait_axis_label = "Trait value", cex.lab = common$cex.lab,
        cex.axis = common$cex.axis, tree_fraction = NULL,
        tree_height = 1.2, moment_height = 1.1, left_margin = 4.2,
        right_margin = 0.6, include_stem = FALSE,
        other_edge_color = "grey80", other_tip_color = "grey50",
        overlap_color = "black", overlap_policy = "mixed"
    )
    tree_options <- utils::modifyList(tree_defaults, tree_style)

    base_args <- list(
        x = x, clades = clades, colors = colors, palette = palette,
        time_start = time_start, xlim = xlim, time_orientation = time_orientation,
        present_time = present_time, bands = bands, n_bands = n_bands,
        band_color = band_color, legend_pos = legend_pos,
        hatch_density = common$hatch_density, hatch_angle = common$hatch_angle,
        hatch_lwd = common$hatch_lwd,
        hatched_outline = common$hatched_outline,
        hatched_outline_lty = common$hatched_outline_lty,
        hatched_outline_lwd = common$hatched_outline_lwd,
        hatched_outline_alpha = common$hatched_outline_alpha,
        cex.lab = common$cex.lab, cex.axis = common$cex.axis,
        cex.legend = common$cex.legend
    )

    if ("tree" %in% panels) {
        args <- c(base_args, list(
            tree = tree, trait = trait, panels = panels,
            show_disparity = show_disparity,
            mean_display = mean_display, mean_envelope = mean_envelope,
            mean_envelope_legend_label = mean_envelope_legend_label,
            mean_envelope_style = mean_envelope_style, variance_envelope_style = variance_envelope_style,
            mean_envelope_probs = mean_envelope_probs,
            variance_envelope_probs = variance_envelope_probs,
            mean_halo_style = mean$halo_style,
            variance_halo_style = variance$halo_style,
            mean_halo_alpha = mean$halo_alpha,
            variance_halo_alpha = variance$halo_alpha,
            variance_display = variance_display,
            mean_expectation_lwd = mean$expectation_lwd,
            mean_envelope_lwd = mean$envelope_lwd,
            mean_envelope_lty = mean$envelope_lty,
            variance_expectation_lwd = variance$expectation_lwd,
            variance_median_lwd = variance$median_lwd,
            variance_median_lty = variance$median_lty,
            mean_ylim = mean$ylim, variance_ylim = variance$ylim,
            mean_ylab = mean$ylab, variance_ylab = variance$ylab,
            disparity_level_ylab = disparity_level$ylab,
            disparity_level_lwd = disparity_level$curve_lwd,
            disparity_level_reference_lwd = disparity_level$reference_lwd,
            disparity_level_reference_lty = disparity_level$reference_lty,
            tree_lwd = tree_options$tree_lwd,
            trait_lwd = tree_options$trait_lwd,
            trait_axis = tree_options$trait_axis,
            trait_axis_label = tree_options$trait_axis_label,
            geological_scale = geological_scale,
            geological_rank = geological_rank,
            geological_axis_args = geological_axis_args,
            geological_background_alpha = geological_background_alpha,
            tree_fraction = tree_options$tree_fraction,
            tree_height = tree_options$tree_height,
            moment_height = tree_options$moment_height,
            left_margin = tree_options$left_margin,
            right_margin = tree_options$right_margin,
            include_stem = tree_options$include_stem,
            other_edge_color = tree_options$other_edge_color,
            other_tip_color = tree_options$other_tip_color,
            overlap_color = tree_options$overlap_color,
            overlap_policy = tree_options$overlap_policy
        ), dots)
        args$clades <- clade_definitions
        args$tree_cex.lab <- tree_options$cex.lab
        args$tree_cex.axis <- tree_options$cex.axis
        return(do.call(.plot_reconstruction_with_tree, args))
    }

    if (all(c("mean", "variance") %in% panels)) {
        args <- c(base_args, list(
            mean_display = mean_display, mean_envelope = mean_envelope,
            mean_envelope_legend_label = mean_envelope_legend_label,
            mean_envelope_style = mean_envelope_style, variance_envelope_style = variance_envelope_style,
            mean_envelope_probs = mean_envelope_probs,
            variance_envelope_probs = variance_envelope_probs,
            mean_halo_style = mean$halo_style,
            variance_halo_style = variance$halo_style,
            mean_halo_alpha = mean$halo_alpha,
            variance_halo_alpha = variance$halo_alpha,
            variance_display = variance_display,
            mean_expectation_lwd = mean$expectation_lwd,
            mean_envelope_lwd = mean$envelope_lwd,
            mean_envelope_lty = mean$envelope_lty,
            variance_expectation_lwd = variance$expectation_lwd,
            variance_median_lwd = variance$median_lwd,
            variance_median_lty = variance$median_lty,
            mean_ylim = mean$ylim, variance_ylim = variance$ylim,
            mean_ylab = mean$ylab, variance_ylab = variance$ylab,
            disparity_level_ylab = disparity_level$ylab,
            disparity_level_lwd = disparity_level$curve_lwd,
            disparity_level_reference_lwd = disparity_level$reference_lwd,
            disparity_level_reference_lty = disparity_level$reference_lty,
            show_disparity = show_disparity,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha
        ), dots)
        return(do.call(.plot_empirical_moment_panels, args))
    }

    if (setequal(panels, c("variance", "disparity_level"))) {
        variance_args <- c(base_args, list(
            interval = variance_envelope_style,
            interval_probs = variance_envelope_probs,
            interval_alpha = variance$halo_alpha,
            halo_style = variance$halo_style,
            show_expectation =
                variance_display %in% c("both", "expectation"),
            show_median = variance_display %in% c("both", "median"),
            summary_legend = variance_display != "none",
            ylim = variance$ylim,
            ylab = variance$ylab,
            expectation_lwd = variance$expectation_lwd,
            median_lwd = variance$median_lwd,
            median_lty = variance$median_lty
        ), dots)
        reference_args <- list(
            x = x, clades = clades, colors = colors, palette = palette,
            ylab = disparity_level$ylab,
            curve_lwd = disparity_level$curve_lwd,
            reference_lwd = disparity_level$reference_lwd,
            reference_lty = disparity_level$reference_lty,
            time_start = time_start, xlim = xlim,
            time_orientation = time_orientation, present_time = present_time,
            bands = bands, n_bands = n_bands, band_color = band_color,
            legend_pos = FALSE,
            cex.lab = common$cex.lab, cex.axis = common$cex.axis,
            cex.legend = common$cex.legend
        )
        return(.plot_variance_reference_panels(
            x = x,
            variance_args = variance_args,
            reference_args = reference_args,
            cex.lab = common$cex.lab,
            cex.axis = common$cex.axis,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha
        ))
    }

    if (setequal(panels, c("mean", "disparity_level"))) {
        mean_args <- c(base_args, list(
            interval = if (mean_envelope == "conditional_quantiles") {
                mean_envelope_style
            } else {
                "none"
            },
            interval_probs = mean_envelope_probs,
            interval_alpha = mean$halo_alpha,
            halo_style = mean$halo_style,
            ylim = mean$ylim,
            ylab = mean$ylab,
            show_expectation = mean_display == "expectation",
            expectation_lwd = mean$expectation_lwd
        ), dots)
        mean_fun <- .plot_empirical_mean_panel
        if (mean_envelope == "expected_sd") {
            mean_fun <- .plot_empirical_spread_panel
            mean_args$interval <- NULL
            mean_args$interval_probs <- NULL
            mean_args$interval_alpha <- NULL
            mean_args$envelope <- mean_envelope_style
            mean_args$spread <- mean_envelope
            mean_args$spread_legend_label <- mean_envelope_legend_label
            mean_args$envelope_alpha <- mean$halo_alpha
            mean_args$envelope_lwd <- mean$envelope_lwd
            mean_args$envelope_lty <- mean$envelope_lty
            mean_args$ylim_expand <- mean$ylim_expand
        }
        disparity_args <- list(
            x = x, clades = clades, colors = colors, palette = palette,
            ylab = disparity_level$ylab,
            curve_lwd = disparity_level$curve_lwd,
            reference_lwd = disparity_level$reference_lwd,
            reference_lty = disparity_level$reference_lty,
            time_start = time_start, xlim = xlim,
            time_orientation = time_orientation, present_time = present_time,
            bands = bands, n_bands = n_bands, band_color = band_color,
            legend_pos = FALSE,
            cex.lab = common$cex.lab, cex.axis = common$cex.axis,
            cex.legend = common$cex.legend
        )
        return(.plot_mean_disparity_panels(
            mean_fun = mean_fun,
            mean_args = mean_args,
            disparity_args = disparity_args,
            cex.lab = common$cex.lab,
            cex.axis = common$cex.axis,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha
        ))
    }

    # A single moment panel otherwise inherits the device's fixed margins.
    # Enlarge them here so increased axis and label sizes remain separated and
    # are not clipped.
    text_layout <- .axis_text_layout(common$cex.lab, common$cex.axis)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    current_mar <- graphics::par("mar")
    geological_margin <- .geological_axis_margin(
        text_layout, common$cex.axis, geological_periods
    )
    graphics::par(
        mar = c(
            if (geological_scale) {
                geological_margin
            } else {
                max(current_mar[1L], text_layout$margin)
            },
            max(current_mar[2L], text_layout$margin),
            current_mar[3L], current_mar[4L]
        ),
        mgp = text_layout$mgp
    )
    finish_single_panel <- function(info) {
        if (geological_scale) {
            .draw_geological_time_axis(
                info$xlim, present_time, geological_periods,
                cex.axis = common$cex.axis, cex.lab = common$cex.lab
            )
        }
        info$geological_scale <- geological_scale
        info$geological_background_alpha <- geological_background_alpha
        invisible(info)
    }

    if (identical(panels, "disparity_level")) {
        args <- c(list(
            x = x, clades = clades, colors = colors, palette = palette,
            ylab = disparity_level$ylab,
            curve_lwd = disparity_level$curve_lwd,
            reference_lwd = disparity_level$reference_lwd,
            reference_lty = disparity_level$reference_lty,
            time_start = time_start, xlim = xlim,
            time_orientation = time_orientation, present_time = present_time,
            bands = bands, n_bands = n_bands, band_color = band_color,
            legend_pos = legend_pos,
            draw_x_axis = !geological_scale,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha,
            cex.lab = common$cex.lab, cex.axis = common$cex.axis,
            cex.legend = common$cex.legend
        ), dots)
        return(finish_single_panel(
            do.call(.plot_disparity_level_panel, args)
        ))
    }

    if (identical(panels, "mean")) {
        args <- c(base_args, list(
            interval = if (mean_envelope == "conditional_quantiles") {
                mean_envelope_style
            } else {
                "none"
            },
            interval_probs = mean_envelope_probs, interval_alpha = mean$halo_alpha,
            halo_style = mean$halo_style, ylim = mean$ylim,
            ylab = mean$ylab,
            show_expectation = mean_display == "expectation",
            expectation_lwd = mean$expectation_lwd,
            draw_x_axis = !geological_scale,
            geological_periods = geological_periods,
            geological_present_time = present_time,
            geological_background_alpha = geological_background_alpha
        ), dots)
        if (mean_envelope == "expected_sd") {
            args$interval <- NULL
            args$interval_probs <- NULL
            args$interval_alpha <- NULL
            args$envelope <- mean_envelope_style
            args$spread <- mean_envelope
            args$spread_legend_label <- mean_envelope_legend_label
            args$envelope_alpha <- mean$halo_alpha
            args$envelope_lwd <- mean$envelope_lwd
            args$envelope_lty <- mean$envelope_lty
            args$ylim_expand <- mean$ylim_expand
            return(finish_single_panel(
                do.call(.plot_empirical_spread_panel, args)
            ))
        }
        return(finish_single_panel(do.call(.plot_empirical_mean_panel, args)))
    }

    args <- c(base_args, list(
        interval = variance_envelope_style, interval_probs = variance_envelope_probs,
        interval_alpha = variance$halo_alpha, halo_style = variance$halo_style,
        show_expectation = variance_display %in% c("both", "expectation"),
        show_median = variance_display %in% c("both", "median"),
        summary_legend = variance_display != "none", ylim = variance$ylim,
        ylab = variance$ylab,
        expectation_lwd = variance$expectation_lwd,
        median_lwd = variance$median_lwd, median_lty = variance$median_lty,
        draw_x_axis = !geological_scale,
        geological_periods = geological_periods,
        geological_present_time = present_time,
        geological_background_alpha = geological_background_alpha
    ), dots)
    finish_single_panel(do.call(.plot_empirical_variance_panel, args))
}
