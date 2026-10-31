## =============================================================================
## Validation and data extraction
## =============================================================================

.validate_reconstruction_plot_data <- function(x) {
    required <- c(
        "time", "clade", "n_lineages", "quantity",
        "expectation", "variance"
    )

    if (!is.data.frame(x)) {
        stop("x must be a data.frame returned by a reconstruction function.")
    }

    missing <- setdiff(required, names(x))
    if (length(missing) > 0L) {
        stop(
            "x is missing required columns: ",
            paste(missing, collapse = ", "), "."
        )
    }

    allowed <- c("empirical_mean", "empirical_variance")
    bad <- setdiff(unique(as.character(x$quantity)), allowed)
    if (length(bad) > 0L) {
        stop("Unknown quantity values: ", paste(bad, collapse = ", "), ".")
    }

    invisible(TRUE)
}


.subset_reconstruction_plot_data <- function(x, quantity, clades = NULL) {
    .validate_reconstruction_plot_data(x)

    d <- x[x$quantity == quantity, , drop = FALSE]
    if (nrow(d) == 0L) {
        stop("x does not contain quantity '", quantity, "'.")
    }

    available <- unique(as.character(d$clade))

    if (is.null(clades)) {
        clades <- available
    } else {
        clades <- as.character(clades)
        missing <- setdiff(clades, available)
        if (length(missing) > 0L) {
            stop(
                "Requested clades absent from x: ",
                paste(missing, collapse = ", "), "."
            )
        }
    }

    d <- d[d$clade %in% clades, , drop = FALSE]
    d$clade <- factor(d$clade, levels = clades)
    d <- d[order(d$clade, d$time), , drop = FALSE]
    d$clade <- as.character(d$clade)
    rownames(d) <- NULL
    d
}


.quantile_column_probabilities <- function(x) {
    qcols <- grep("^q[0-9]+(_[0-9]+)?$", names(x), value = TRUE)
    if (length(qcols) == 0L) {
        return(setNames(numeric(0), character(0)))
    }

    txt <- sub("^q", "", qcols)
    txt <- sub("_", ".", txt, fixed = TRUE)
    probs <- suppressWarnings(as.numeric(txt) / 100)

    has_values <- vapply(qcols, function(nm) any(is.finite(x[[nm]])), logical(1))
    keep <- is.finite(probs) & probs > 0 & probs < 1 & has_values
    probs <- probs[keep]
    qcols <- qcols[keep]

    ord <- order(probs)
    setNames(probs[ord], qcols[ord])
}


.match_quantile_column <- function(qprobs, probability, tol = 1e-10) {
    if (length(qprobs) == 0L) {
        return(NA_character_)
    }

    idx <- which(abs(unname(qprobs) - probability) <= tol)
    if (length(idx) == 0L) {
        return(NA_character_)
    }

    names(qprobs)[idx[[1L]]]
}


.select_interval_columns <- function(x, interval_probs = NULL) {
    qprobs <- .quantile_column_probabilities(x)

    if (length(qprobs) < 2L) {
        return(list(lower = NA_character_, upper = NA_character_, probs = numeric(0)))
    }

    if (is.null(interval_probs)) {
        selected <- range(unname(qprobs))
    } else {
        interval_probs <- sort(as.numeric(interval_probs))
        if (length(interval_probs) != 2L ||
            any(!is.finite(interval_probs)) ||
            any(interval_probs <= 0) || any(interval_probs >= 1) ||
            interval_probs[1L] >= interval_probs[2L]) {
            stop("interval_probs must contain two increasing probabilities in (0, 1).")
        }
        selected <- interval_probs
    }

    lower <- .match_quantile_column(qprobs, selected[1L])
    upper <- .match_quantile_column(qprobs, selected[2L])

    if (is.na(lower) || is.na(upper)) {
        stop(
            "Requested interval probabilities are unavailable. Available probabilities: ",
            paste(format(unname(qprobs), trim = TRUE), collapse = ", "), "."
        )
    }

    list(lower = lower, upper = upper, probs = selected)
}


.select_median_column <- function(x) {
    .match_quantile_column(.quantile_column_probabilities(x), 0.5)
}


## =============================================================================
## Colours, ranges and drawing helpers
## =============================================================================

#' Assign colours to reconstructed clades
#'
#' @param clades Character vector of unique clade names.
#' @param colors Optional vector of colours, either unnamed or named by clade.
#' @param palette Name of an HCL palette used when \code{colors} is \code{NULL}.
#' @param alpha Alpha transparency passed to \code{grDevices::hcl.colors}.
#' @return A named character vector of colours.
#' @export
make_clade_colors <- function(clades, colors = NULL,
                              palette = "Dark 3", alpha = 1) {
    clades <- as.character(clades)
    if (length(clades) == 0L || anyNA(clades) || any(!nzchar(clades))) {
        stop("clades must contain non-empty names.")
    }
    if (anyDuplicated(clades)) {
        stop("clade names must be unique.")
    }

    if (is.null(colors)) {
        out <- grDevices::hcl.colors(
            length(clades), palette = palette, alpha = alpha
        )
        names(out) <- clades
        return(out)
    }

    if (is.null(names(colors))) {
        if (length(colors) < length(clades)) {
            stop("Not enough colors for the requested clades.")
        }
        out <- colors[seq_along(clades)]
        names(out) <- clades
        return(out)
    }

    missing <- setdiff(clades, names(colors))
    if (length(missing) > 0L) {
        stop("Missing colors for clades: ", paste(missing, collapse = ", "), ".")
    }

    colors[clades]
}


.plot_finite_runs <- function(
    time,
    lower,
    upper,
    color,
    border = NA,
    halo_style = c("solid", "hatched"),
    hatch_density = 18,
    hatch_angle = 45,
    hatch_lwd = 0.7
) {
    halo_style <- match.arg(halo_style)

    if (!is.numeric(hatch_density) || length(hatch_density) != 1L ||
        !is.finite(hatch_density) || hatch_density <= 0) {
        stop("hatch_density must be a strictly positive finite scalar.")
    }
    if (!is.numeric(hatch_angle) || length(hatch_angle) != 1L ||
        !is.finite(hatch_angle)) {
        stop("hatch_angle must be a finite scalar.")
    }
    if (!is.numeric(hatch_lwd) || length(hatch_lwd) != 1L ||
        !is.finite(hatch_lwd) || hatch_lwd <= 0) {
        stop("hatch_lwd must be a strictly positive finite scalar.")
    }

    ok <- is.finite(time) & is.finite(lower) & is.finite(upper)
    if (!any(ok)) {
        return(invisible(NULL))
    }

    idx <- which(ok)
    groups <- cumsum(c(TRUE, diff(idx) != 1L))

    for (g in unique(groups)) {
        ii <- idx[groups == g]
        if (length(ii) < 2L) {
            next
        }

        if (halo_style == "solid") {
            graphics::polygon(
                c(time[ii], rev(time[ii])),
                c(lower[ii], rev(upper[ii])),
                col = color,
                border = border
            )
        } else {
            ## With a non-NULL density, polygon() draws hatch lines rather than
            ## a solid fill. The alpha channel of `color` is retained on devices
            ## supporting transparency.
            graphics::polygon(
                c(time[ii], rev(time[ii])),
                c(lower[ii], rev(upper[ii])),
                col = color,
                border = border,
                density = hatch_density,
                angle = hatch_angle,
                lwd = hatch_lwd
            )
        }
    }

    invisible(NULL)
}


.compute_plot_ylim <- function(d, interval, interval_columns,
                               show_median = FALSE, median_column = NA_character_,
                               ylim_expand = 0.06) {
    values <- d$expectation

    if (interval != "none" && !is.na(interval_columns$lower)) {
        values <- c(
            values,
            d[[interval_columns$lower]],
            d[[interval_columns$upper]]
        )
    }

    if (show_median && !is.na(median_column)) {
        values <- c(values, d[[median_column]])
    }

    values <- values[is.finite(values)]
    if (length(values) == 0L) {
        stop("No finite values are available for plotting.")
    }

    r <- range(values)
    if (r[1L] == r[2L]) {
        delta <- max(1, abs(r[1L])) * 0.05
        return(r + c(-delta, delta))
    }

    grDevices::extendrange(r, f = ylim_expand)
}


.make_time_bands <- function(xlim, n_bands = 0L) {
    if (is.null(n_bands) || identical(n_bands, 0L) || identical(n_bands, 0)) {
        return(NULL)
    }

    if (!is.numeric(n_bands) || length(n_bands) != 1L ||
        !is.finite(n_bands) || n_bands < 1 || n_bands != as.integer(n_bands)) {
        stop("n_bands must be a nonnegative integer.")
    }

    n_bands <- as.integer(n_bands)
    breaks <- seq(xlim[1L], xlim[2L], length.out = n_bands + 1L)
    idx <- seq_len(n_bands)
    idx <- idx[idx %% 2L == 1L]

    data.frame(
        start = breaks[idx],
        end = breaks[idx + 1L],
        stringsAsFactors = FALSE
    )
}


.validate_bands <- function(bands) {
    if (is.null(bands)) {
        return(NULL)
    }
    if (!is.data.frame(bands) || !all(c("start", "end") %in% names(bands))) {
        stop("bands must be NULL or a data.frame with columns 'start' and 'end'.")
    }
    if (any(!is.finite(bands$start)) || any(!is.finite(bands$end)) ||
        any(bands$end < bands$start)) {
        stop("bands must contain finite intervals with end >= start.")
    }
    bands
}


.draw_time_bands <- function(bands, color) {
    if (is.null(bands) || nrow(bands) == 0L) {
        return(invisible(NULL))
    }

    usr <- graphics::par("usr")
    for (i in seq_len(nrow(bands))) {
        graphics::rect(
            bands$start[i], usr[3L], bands$end[i], usr[4L],
            col = color, border = NA
        )
    }

    invisible(NULL)
}


.draw_time_axis <- function(xlim, time_orientation = c("from_root", "before_present"),
                            present_time = NULL, side = 1L, n = 6L,
                            digits = 4L, ...) {
    time_orientation <- match.arg(time_orientation)

    if (time_orientation == "from_root") {
        at <- pretty(xlim, n = n)
        at <- at[at >= xlim[1L] & at <= xlim[2L]]
        labels <- format(round(at, digits), trim = TRUE)
    } else {
        if (is.null(present_time)) {
            present_time <- xlim[2L]
        }
        if (!is.numeric(present_time) || length(present_time) != 1L ||
            !is.finite(present_time)) {
            stop("present_time must be a finite scalar.")
        }

        displayed_range <- sort(present_time - xlim)
        displayed_ticks <- pretty(displayed_range, n = n)
        displayed_ticks <- displayed_ticks[
            displayed_ticks >= displayed_range[1L] &
                displayed_ticks <= displayed_range[2L]
        ]
        at <- present_time - displayed_ticks
        tick_order <- order(at)
        at <- at[tick_order]
        labels <- format(round(displayed_ticks[tick_order], digits), trim = TRUE)
    }

    graphics::axis(side = side, at = at, labels = labels, ...)
    invisible(list(at = at, labels = labels))
}


.resolve_plot_bands <- function(xlim, bands, n_bands) {
    bands <- .validate_bands(bands)
    if (is.null(bands)) {
        bands <- .make_time_bands(xlim, n_bands = n_bands)
    }
    bands
}


## =============================================================================
## Shared curve plotter
## =============================================================================

.resolve_time_window <- function(time, time_start = NULL, xlim = NULL) {
    finite_time <- time[is.finite(time)]
    if (length(finite_time) == 0L) {
        stop("No finite time value is available for plotting.")
    }

    if (is.null(time_start)) {
        time_start <- 0
    } else if (!is.numeric(time_start) || length(time_start) != 1L ||
               !is.finite(time_start) || time_start < 0) {
        stop("time_start must be a finite nonnegative scalar or NULL.")
    }

    if (is.null(xlim)) {
        xlim <- c(time_start, max(finite_time))
    } else {
        if (!is.numeric(xlim) || length(xlim) != 2L ||
            any(!is.finite(xlim)) || xlim[1L] >= xlim[2L]) {
            stop("xlim must contain two increasing finite values.")
        }
        ## time_start is a semantic lower truncation and therefore cannot be
        ## undone by supplying a smaller lower x limit.
        xlim[1L] <- max(xlim[1L], time_start)
    }

    if (xlim[1L] >= xlim[2L]) {
        stop("time_start must be smaller than the upper plotting limit.")
    }
    if (!any(finite_time >= xlim[1L] & finite_time <= xlim[2L])) {
        stop("No reconstruction time falls inside the requested plotting window.")
    }

    xlim
}

.trim_plot_data_to_window <- function(d, xlim) {
    keep <- is.finite(d$time) & d$time >= xlim[1L] & d$time <= xlim[2L]
    d <- d[keep, , drop = FALSE]
    if (nrow(d) == 0L) {
        stop("No reconstruction row falls inside the requested plotting window.")
    }
    rownames(d) <- NULL
    d
}


.plot_empirical_quantity <- function(
    x,
    quantity,
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
    show_expectation = TRUE,
    show_median = FALSE,
    expectation_lwd = 2,
    median_lwd = 1.3,
    median_lty = 2,
    interval_lwd = 0.8,
    interval_lty = 3,
    time_start = NULL,
    xlim = NULL,
    ylim = NULL,
    xlab = NULL,
    ylab = NULL,
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
    summary_legend = TRUE,
    summary_legend_pos = "topleft",
    cex.lab = 1.35,
    cex.axis = 1.0,
    cex.legend = 1.0,
    ylim_expand = 0.06,
    ...
) {
    interval <- match.arg(interval)
    halo_style <- match.arg(halo_style)
    time_orientation <- match.arg(time_orientation)
    if (!is.logical(hatched_outline) || length(hatched_outline) != 1L ||
        is.na(hatched_outline)) {
        stop("hatched_outline must be TRUE or FALSE.")
    }
    if (!is.numeric(hatched_outline_lty) || length(hatched_outline_lty) != 1L ||
        !is.finite(hatched_outline_lty)) {
        stop("hatched_outline_lty must be a finite scalar.")
    }
    if (!is.numeric(hatched_outline_lwd) || length(hatched_outline_lwd) != 1L ||
        !is.finite(hatched_outline_lwd) || hatched_outline_lwd <= 0) {
        stop("hatched_outline_lwd must be a strictly positive finite scalar.")
    }
    if (!is.numeric(hatched_outline_alpha) || length(hatched_outline_alpha) != 1L ||
        !is.finite(hatched_outline_alpha) || hatched_outline_alpha < 0 ||
        hatched_outline_alpha > 1) {
        stop("hatched_outline_alpha must be a finite scalar in [0, 1].")
    }

    d <- .subset_reconstruction_plot_data(x, quantity, clades = clades)
    clade_names <- unique(d$clade)
    cols <- make_clade_colors(clade_names, colors = colors, palette = palette)

    interval_columns <- .select_interval_columns(d, interval_probs)
    if (interval != "none" && is.na(interval_columns$lower)) {
        interval <- "none"
    }

    median_column <- .select_median_column(d)
    if (show_median && is.na(median_column)) {
        show_median <- FALSE
    }

    xlim <- .resolve_time_window(d$time, time_start = time_start, xlim = xlim)
    d <- .trim_plot_data_to_window(d, xlim)

    if (is.null(ylim)) {
        ylim <- .compute_plot_ylim(
            d,
            interval = interval,
            interval_columns = interval_columns,
            show_median = show_median,
            median_column = median_column,
            ylim_expand = ylim_expand
        )
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

        if (interval %in% c("ribbon", "both")) {
            .plot_finite_runs(
                dc$time,
                dc[[interval_columns$lower]],
                dc[[interval_columns$upper]],
                color = grDevices::adjustcolor(col, alpha.f = interval_alpha),
                halo_style = halo_style,
                hatch_density = hatch_density,
                hatch_angle = hatch_angle,
                hatch_lwd = hatch_lwd
            )
            if (halo_style == "hatched" && isTRUE(hatched_outline) &&
                interval %in% c("ribbon", "both")) {
                graphics::lines(
                    dc$time, dc[[interval_columns$lower]],
                    col = grDevices::adjustcolor(col, alpha.f = hatched_outline_alpha),
                    lwd = hatched_outline_lwd, lty = hatched_outline_lty
                )
                graphics::lines(
                    dc$time, dc[[interval_columns$upper]],
                    col = grDevices::adjustcolor(col, alpha.f = hatched_outline_alpha),
                    lwd = hatched_outline_lwd, lty = hatched_outline_lty
                )
            }
        }

        if (interval %in% c("lines", "both")) {
            graphics::lines(
                dc$time, dc[[interval_columns$lower]],
                col = grDevices::adjustcolor(col, alpha.f = 0.65),
                lwd = interval_lwd, lty = interval_lty
            )
            graphics::lines(
                dc$time, dc[[interval_columns$upper]],
                col = grDevices::adjustcolor(col, alpha.f = 0.65),
                lwd = interval_lwd, lty = interval_lty
            )
        }

        if (show_expectation) {
            graphics::lines(
                dc$time, dc$expectation,
                col = col, lwd = expectation_lwd, lty = 1
            )
        }

        if (show_median) {
            graphics::lines(
                dc$time, dc[[median_column]],
                col = col, lwd = median_lwd, lty = median_lty
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

    if (!is.null(legend_pos) && !identical(legend_pos, FALSE)) {
        graphics::legend(
            legend_pos,
            legend = clade_names,
            col = unname(cols[clade_names]),
            lty = 1,
            lwd = expectation_lwd,
            bty = "n",
            cex = cex.legend
        )
    }

    if ((show_expectation || show_median) && isTRUE(summary_legend) &&
        !is.null(summary_legend_pos) && !identical(summary_legend_pos, FALSE)) {
        summary_labels <- c(
            if (show_expectation) "Expectation",
            if (show_median) "Median"
        )
        graphics::legend(
            summary_legend_pos,
            legend = summary_labels,
            col = "black",
            lty = c(
                if (show_expectation) 1,
                if (show_median) median_lty
            ),
            lwd = c(
                if (show_expectation) expectation_lwd,
                if (show_median) median_lwd
            ),
            bty = "n",
            cex = cex.legend
        )
    }

    invisible(list(
        colors = cols,
        xlim = xlim,
        ylim = ylim,
        interval = interval,
        interval_columns = interval_columns,
        median_column = median_column,
        show_expectation = show_expectation,
        show_median = show_median,
        bands = bands,
        legend_pos = legend_pos
    ))
}
