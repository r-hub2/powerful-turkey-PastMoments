test_that("basic plotting methods run", {
    tree <- example_tree()
    trait <- example_trait()
    result <- reconstruct_empirical_moments(
        tree,
        trait,
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        mean_probs = c(0.25, 0.5, 0.75),
        variance_probs = NULL
    )

    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    expect_invisible(plot(result, panels = "mean"))
    expect_invisible(plot(
        result, panels = "variance", variance_envelope_style = "none"
    ))
    expect_invisible(plot(
        result, panels = c("mean", "variance"),
        variance_envelope_style = "none"
    ))
})


test_that("clade colours are named consistently", {
    cols <- make_clade_colors(c("A", "B"))
    expect_named(cols, c("A", "B"))
    expect_length(cols, 2L)
})


test_that("time-before-present axis uses round displayed values", {
    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.13, 9.47), ylim = c(0, 1))

    present_time <- 9.73
    axis_info <- PastMoments:::.draw_time_axis(
        c(0.13, 9.47),
        time_orientation = "before_present",
        present_time = present_time
    )
    displayed_ticks <- as.numeric(axis_info$labels)

    expect_equal(axis_info$at + displayed_ticks, rep(present_time, length(displayed_ticks)))
    expect_true(all(displayed_ticks == round(displayed_ticks)))
})

test_that("geological margin accounts for every stacked scale", {
    text_layout <- PastMoments:::.axis_text_layout(3, 2.5)
    one_scale <- data.frame(min_ma = 0, max_ma = 1)
    two_scales <- list(one_scale, one_scale)

    one_margin <- PastMoments:::.geological_axis_margin(
        text_layout, 2.5, one_scale
    )
    two_margin <- PastMoments:::.geological_axis_margin(
        text_layout, 2.5, two_scales
    )

    expect_equal(two_margin - one_margin, 1.5 * 2.5)
    expect_gt(two_margin, text_layout$margin)
})

test_that("geological strip height follows typography, not device height", {
    periods <- data.frame(min_ma = 0, max_ma = 1)
    attr(periods, "axis_args") <- list(lab_size = 1.6)
    physical_heights <- numeric(2L)

    for (i in seq_along(physical_heights)) {
        filename <- tempfile(fileext = ".pdf")
        grDevices::pdf(filename, width = 7, height = c(6, 11)[i])
        graphics::plot.new()
        graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
        fraction <- PastMoments:::.geological_strip_height(
            periods, cex.axis = 1, cex.lab = 1.35
        )
        physical_heights[i] <- fraction * graphics::par("pin")[2L]
        grDevices::dev.off()
        unlink(filename)
    }

    expect_equal(physical_heights[1L], physical_heights[2L])
    expect_gt(physical_heights[1L], 0)
})

test_that("common halo options apply to both moment panels", {
    tree <- example_tree()
    result <- reconstruct_empirical_moments(
        tree,
        example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        variance_probs = NULL
    )

    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    info <- plot(
        result,
        variance_envelope_style = "none",
        common_style = list(halo_style = "hatched", halo_alpha = 0.25)
    )
    expect_equal(info$halo_style, c(mean = "hatched", variance = "hatched"))
    expect_equal(info$halo_alpha, c(mean = 0.25, variance = 0.25))

    overridden <- plot(
        result,
        variance_envelope_style = "none",
        common_style = list(halo_style = "hatched", halo_alpha = 0.25),
        mean_style = list(halo_alpha = 0.4),
        variance_style = list(halo_style = "solid")
    )
    expect_equal(overridden$halo_style,
                 c(mean = "hatched", variance = "solid"))
    expect_equal(overridden$halo_alpha, c(mean = 0.4, variance = 0.25))

    with_tree <- plot(
        result,
        variance_envelope_style = "none",
        common_style = list(halo_style = "hatched", halo_alpha = 0.25)
    )
    expect_equal(with_tree$halo_style,
                 c(mean = "hatched", variance = "hatched"))
    expect_equal(with_tree$halo_alpha, c(mean = 0.25, variance = 0.25))
    expect_equal(with_tree$layout_heights, c(1.2, 1.1, 1.1))
    expect_equal(with_tree$trait_axis_label, "Trait value")
    expect_equal(with_tree$colors, c(all = "black"))
    expect_false(with_tree$clade_legend_pos)

    without_tree <- plot(result, panels = "mean")
    expect_equal(without_tree$colors, c(all = "black"))
})


test_that("trait axis label is validated", {
    tree <- example_tree()
    result <- reconstruct_empirical_moments(
        tree,
        example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        variance_probs = NULL
    )

    expect_error(
        plot(
            result,
            tree_style = list(trait_axis_label = c("a", "b"))
        ),
        "trait_axis_label"
    )
})

test_that("tree axis typography can override common typography", {
    result <- reconstruct_empirical_moments(
        example_tree(),
        example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        variance_probs = NULL
    )
    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename, width = 9, height = 9)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    info <- plot(
        result,
        variance_envelope_style = "none",
        common_style = list(cex.lab = 1.2, cex.axis = 0.9),
        tree_style = list(cex.lab = 2.1, cex.axis = 1.7)
    )

    expect_equal(
        info$moment_text_layout,
        PastMoments:::.axis_text_layout(1.2, 0.9)
    )
    expect_equal(
        info$trait_text_layout,
        PastMoments:::.axis_text_layout(2.1, 1.7)
    )
})

test_that("axis typography enlarges spacing and works in every layout", {
    normal <- PastMoments:::.axis_text_layout(1.35, 1)
    large <- PastMoments:::.axis_text_layout(2.5, 2)
    expect_gte(normal$mgp[2L], 1)
    expect_gt(large$mgp[1L], normal$mgp[1L])
    expect_gt(large$mgp[2L], normal$mgp[2L])
    expect_gt(large$margin, normal$margin)

    result <- reconstruct_empirical_moments(
        example_tree(),
        example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        variance_probs = NULL
    )
    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename, width = 9, height = 9)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    large_text <- list(cex.lab = 2.5, cex.axis = 2, cex.legend = 1.8)
    expect_invisible(plot(
        result, panels = "mean", common_style = large_text
    ))
    expect_invisible(plot(
        result, panels = c("mean", "variance"),
        variance_envelope_style = "none", common_style = large_text
    ))
    expect_invisible(plot(
        result, variance_envelope_style = "none",
        common_style = large_text,
        tree_style = list(trait_axis_label = "Measured trait")
    ))
})

test_that("disparity is an explicit optional panel", {
    result <- reconstruct_empirical_moments(
        example_tree(), example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1, variance_probs = NULL,
        quantities = c("mean", "variance", "disparity_level")
    )
    without_reference <- reconstruct_empirical_moments(
        example_tree(), example_trait(), times = c(0.5, 1.5),
        sigma2 = 1, quantities = "variance"
    )
    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename, width = 9, height = 11)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    default <- plot(result, variance_envelope_style = "none")
    expect_null(default$disparity_level)
    expect_equal(default$layout_heights, c(1.2, 1.1, 1.1))

    tree_panels <- plot(
        result, panels = c("tree", "mean", "variance", "disparity_level"),
        variance_envelope_style = "none"
    )
    expect_equal(tree_panels$layout_heights, c(1.2, 1.1, 1.1, 1.1))
    expect_equal(tree_panels$disparity_level$ylim, c(-0.05, 1.05))

    variance_panels <- plot(
        result, panels = c("variance", "disparity_level"),
        variance_envelope_style = "none"
    )
    expect_equal(variance_panels$disparity_level$ylim, c(-0.05, 1.05))

    moment_panels <- plot(
        result, panels = c("mean", "variance", "disparity_level"),
        variance_envelope_style = "none"
    )
    expect_equal(moment_panels$disparity_level$ylim, c(-0.05, 1.05))

    expect_error(
        plot(without_reference, panels = "disparity_level"),
        "quantities including 'disparity_level'"
    )
    disparity <- plot(
        result,
        panels = "disparity_level",
        xlab = "Time before present (Ma)"
    )
    expect_equal(disparity$ylim, c(-0.05, 1.05))
    expect_equal(disparity$colors, c(all = "black"))

    mean_disparity <- plot(result, panels = c("mean", "disparity_level"))
    expect_equal(mean_disparity$disparity_level$ylim, c(-0.05, 1.05))
})

test_that("the first non-tree panel carries the clade legend", {
    result <- reconstruct_empirical_moments(
        example_tree(), example_trait(),
        times = seq(1, 2, length.out = 8),
        clades = list(left = c("a", "b"), right = c("c", "d")),
        sigma2 = 1,
        variance_probs = NULL,
        quantities = c("mean", "variance", "disparity_level")
    )
    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename, width = 9, height = 11)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    variance_first <- plot(
        result,
        panels = c("tree", "variance"),
        variance_envelope_style = "none"
    )
    expect_identical(variance_first$variance$legend_pos, "topright")

    disparity_first <- plot(
        result,
        panels = c("tree", "disparity_level")
    )
    expect_identical(disparity_first$disparity_level$legend_pos, "topright")

    mean_first <- plot(
        result,
        panels = c("tree", "mean", "variance", "disparity_level"),
        variance_envelope_style = "none"
    )
    expect_identical(mean_first$mean$legend_pos, "topright")
    expect_false(mean_first$variance$legend_pos)
    expect_false(mean_first$disparity_level$legend_pos)
})


test_that("empirical-spread plotting choices are explicit and compatible", {
    tree <- example_tree()
    result <- reconstruct_empirical_moments(
        tree,
        example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        variance_probs = NULL
    )

    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    expect_error(
        plot(
            result, panels = "mean",
            mean_envelope = "sqrt_expected_variance",
            legend_pos = FALSE
        ),
        "'arg' should be one of"
    )

    expected <- plot(
        result,
        panels = "mean",
        mean_envelope = "expected_sd",
        mean_envelope_legend_label = "Expected empirical SD",
        mean_envelope_style = "none",
        legend_pos = FALSE
    )
    expect_identical(expected$envelope_label, "Expected empirical SD")
    expect_equal(expected$data$envelope, expected$data$expected_sd)
    expect_true(all(expected$data$envelope >= 0, na.rm = TRUE))
    labelled <- plot(
        result,
        panels = "mean",
        mean_envelope = "expected_sd",
        mean_envelope_legend_label = "Expected empirical SD"
    )
    expect_identical(labelled$envelope_label, "Expected empirical SD")
    expect_error(
        plot(
            result,
            panels = "mean",
            mean_envelope = "expected_sd",
            mean_envelope_legend_label = TRUE
        ),
        "must be NULL or a single character string"
    )

    suppressed <- plot(
        result,
        panels = "mean",
        mean_envelope = "none",
        mean_envelope_style = "ribbon",
        legend_pos = FALSE
    )
    expect_identical(suppressed$interval, "none")
    expect_equal(suppressed$data$lower, suppressed$data$mean_expectation)
    expect_equal(suppressed$data$upper, suppressed$data$mean_expectation)

    with_tree <- plot(
        result,
        mean_envelope = "expected_sd",
        mean_envelope_style = "none",
        variance_envelope_style = "none"
    )
    expect_identical(with_tree$mean$envelope, "expected_sd")
})


test_that("central curves and envelopes are selected independently", {
    result <- reconstruct_empirical_moments(
        example_tree(), example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1,
        variance_probs = c(0.25, 0.5, 0.75)
    )

    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    mean_envelope_only <- plot(
        result, panels = "mean", mean_display = "none",
        mean_envelope = "conditional_quantiles"
    )
    expect_false(mean_envelope_only$show_expectation)
    expect_identical(mean_envelope_only$interval, "ribbon")

    variance_median_only <- plot(
        result, panels = "variance", variance_display = "median",
        variance_envelope_style = "none"
    )
    expect_false(variance_median_only$show_expectation)
    expect_true(variance_median_only$show_median)

    variance_envelope_only <- plot(
        result, panels = "variance", variance_display = "none",
        variance_envelope_style = "ribbon"
    )
    expect_false(variance_envelope_only$show_expectation)
    expect_false(variance_envelope_only$show_median)
    expect_identical(variance_envelope_only$interval, "ribbon")
})


test_that("panel selection and the reduced public API are enforced", {
    tree <- example_tree()
    result <- reconstruct_empirical_moments(
        tree, example_trait(), times = seq(0.25, 2, length.out = 8),
        sigma2 = 1, variance_probs = NULL
    )
    exports <- getNamespaceExports("PastMoments")
    expect_false(any(c(
        "plot_empirical_mean", "plot_empirical_variance",
        "plot_empirical_moments", "plot_empirical_moments_with_tree",
        "plot_empirical_mean_with_empirical_sd"
    ) %in% exports))

    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    tree_only <- plot(result, panels = "tree")
    expect_identical(tree_only$panels, "tree")
    expect_identical(tree_only$tree_fraction, 1)

    selected <- plot(
        result,
        panels = c("tree", "variance"), variance_envelope_style = "none"
    )
    expect_identical(selected$panels, c("tree", "variance"))
    expect_null(selected$mean)
    expect_true(is.list(selected$variance))

    modified <- result
    modified$trait <- modified$trait[-1L]
    expect_error(plot(modified, panels = "tree"), "missing values")
    expect_error(
        plot(result, tree = tree),
        "stored in x"
    )
    expect_error(
        plot(result, mean_style = list(unknown = 1)),
        "Unknown mean_style"
    )
    expect_error(
        plot(result, disparity_level_style = list(unknown = 1)),
        "Unknown disparity_level_style"
    )

    expect_silent(plot(
        result,
        panels = c("tree", "mean", "variance", "disparity_level"),
        mean_style = list(ylab = "Custom mean"),
        variance_style = list(ylab = "Custom variance"),
        disparity_level_style = list(
            ylab = "Custom disparity", curve_lwd = 3,
            reference_lwd = 2, reference_lty = 3
        )
    ))
})


test_that("geological scale validates its plotting context", {
    result <- reconstruct_empirical_moments(
        example_tree(), example_trait(),
        times = seq(0.25, 2, length.out = 8),
        sigma2 = 1, variance_probs = NULL
    )

    filename <- tempfile(fileext = ".pdf")
    grDevices::pdf(filename)
    on.exit({
        grDevices::dev.off()
        unlink(filename)
    }, add = TRUE)

    expect_error(
        plot(result, geological_scale = NA),
        "must be TRUE or FALSE"
    )
    expect_error(
        plot(result, geological_background_alpha = 2),
        "must be a finite scalar in \\[0, 1\\]"
    )
    expect_error(
        plot(result, geological_rank = "age"),
        "geological_rank must be one or more"
    )
    expect_error(
        plot(result, geological_rank = c("period", "period")),
        "unique values"
    )
    expect_error(
        plot(result, geological_axis_args = list(side = 2)),
        "cannot override layout-controlled argument"
    )
    expect_error(
        plot(result, geological_axis_args = list(1)),
        "must be a named list"
    )
    expect_error(
        plot(result, geological_scale = TRUE),
        "requires time_orientation"
    )

    if (!requireNamespace("palaeoverse", quietly = TRUE)) {
        expect_warning(
            plot(
                result, time_orientation = "before_present",
                geological_scale = TRUE
            ),
            "package 'palaeoverse' is not installed"
        )
    } else {
        info <- plot(
            result, time_orientation = "before_present",
            geological_scale = TRUE
        )
        expect_true(info$geological_scale)
        expect_equal(info$geological_background_alpha, 0.12)

        stage_info <- plot(
            result, time_orientation = "before_present",
            geological_scale = TRUE,
            geological_rank = c("era", "period", "stage"),
            geological_axis_args = list(abbr = TRUE, rot = 45, lab_size = 0.8)
        )
        expect_true(stage_info$geological_scale)

        moment_info <- plot(
            result, panels = c("mean", "variance"),
            time_orientation = "before_present", geological_scale = TRUE
        )
        expect_true(moment_info$geological_scale)

        disparity_result <- reconstruct_empirical_moments(
            example_tree(), example_trait(),
            times = seq(0.25, 2, length.out = 8),
            sigma2 = 1, variance_probs = NULL,
            quantities = c("mean", "variance", "disparity_level")
        )
        disparity_info <- plot(
            disparity_result, panels = "disparity_level",
            time_orientation = "before_present", geological_scale = TRUE
        )
        expect_true(disparity_info$geological_scale)

        variance_disparity_info <- plot(
            disparity_result, panels = c("variance", "disparity_level"),
            time_orientation = "before_present", geological_scale = TRUE
        )
        expect_true(variance_disparity_info$geological_scale)
        expect_true(is.list(variance_disparity_info$variance))
        expect_true(is.list(variance_disparity_info$disparity_level))
    }
})
