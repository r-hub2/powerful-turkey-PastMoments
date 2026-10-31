test_that("identical continuous spectral distributions give one half", {
    answer <- .variance_brownian_superiority(
        conditional_weights = c(0.3, 1.7),
        conditional_noncentrality = c(0, 0),
        unconditional_weights = c(0.3, 1.7),
        n = 3
    )
    expect_equal(answer$probability, 0.5, tolerance = 2e-5)
    expect_true(answer$convergence)
    expect_identical(answer$method, "davies")
})


test_that("the probability uses the upper tail and deterministic shifts", {
    larger <- .variance_brownian_superiority(
        conditional_weights = 0.01,
        unconditional_weights = 1,
        conditional_shift = 100,
        n = 2
    )
    smaller <- .variance_brownian_superiority(
        conditional_weights = 1,
        unconditional_weights = 0.01,
        unconditional_shift = 100,
        n = 2
    )
    expect_gt(larger$probability, 0.999)
    expect_lt(smaller$probability, 0.001)
})


test_that("indefinite calculation agrees with one-dimensional integration", {
    conditional_weight <- 0.7
    unconditional_weight <- 1.3
    ncp <- 1.2
    shift <- 0.15
    answer <- .variance_brownian_superiority(
        conditional_weights = conditional_weight,
        conditional_noncentrality = ncp,
        unconditional_weights = unconditional_weight,
        conditional_shift = shift,
        n = 2,
        accuracy = 1e-8
    )
    reference <- stats::integrate(
        function(z) {
            stats::pchisq(
                (conditional_weight * z + shift) /
                    unconditional_weight,
                df = 1
            ) * stats::dchisq(z, df = 1, ncp = ncp)
        },
        lower = 0, upper = Inf, rel.tol = 1e-9
    )$value
    expect_equal(answer$probability, reference, tolerance = 2e-6)
})


test_that("deterministic result agrees with Monte Carlo", {
    set.seed(271828)
    conditional_weights <- c(0.2, 1.4)
    conditional_ncp <- c(0.5, 2)
    unconditional_weights <- c(0.6, 1.1)
    answer <- .variance_brownian_superiority(
        conditional_weights,
        conditional_ncp,
        unconditional_weights,
        n = 3
    )
    draws <- 200000L
    conditional <- conditional_weights[1L] *
        stats::rchisq(draws, 1, conditional_ncp[1L]) +
        conditional_weights[2L] *
        stats::rchisq(draws, 1, conditional_ncp[2L])
    unconditional <- unconditional_weights[1L] *
        stats::rchisq(draws, 1) +
        unconditional_weights[2L] *
        stats::rchisq(draws, 1)
    monte_carlo <- mean(conditional > unconditional)
    mc_se <- sqrt(monte_carlo * (1 - monte_carlo) / draws)
    expect_equal(
        answer$probability, monte_carlo,
        tolerance = 5 * mc_se + 5e-4
    )
})


test_that("degenerate and negligible spectra are handled explicitly", {
    one <- .variance_brownian_superiority(1, 0, 1, n = 1)
    expect_true(is.na(one$probability))

    zero <- .variance_brownian_superiority(
        numeric(), numeric(), numeric(), n = 2
    )
    expect_true(is.na(zero$probability))
    expect_identical(zero$method, "degenerate")
    tied_nonzero <- .variance_brownian_superiority(
        numeric(), numeric(), numeric(),
        conditional_shift = 2, unconditional_shift = 2, n = 2
    )
    expect_identical(tied_nonzero$probability, 0)

    near_zero <- .variance_brownian_superiority(
        c(1, 1e-16), c(0, 1e12), c(1), n = 3, tol = 1e-12
    )
    expect_length(near_zero$conditional_weights, 1L)
    expect_gt(near_zero$probability, 0.5)
    expect_lt(near_zero$probability, 0.501)
})

test_that("negligible noncentral components retain their deterministic limit", {
    limiting_probability <- stats::pchisq(1, df = 1)
    probabilities <- vapply(
        c(1.1e-12, 0.9e-12, 0.5e-12),
        function(weight) {
            answer <- .variance_brownian_superiority(
                conditional_weights = weight,
                conditional_noncentrality = 1 / weight,
                unconditional_weights = 1,
                n = 2,
                tol = 1e-12
            )
            answer$probability
        },
        numeric(1)
    )

    expect_equal(probabilities, rep(limiting_probability, 3), tolerance = 1e-4)
    expect_lt(diff(range(probabilities)), 1e-4)

    deterministic <- .variance_brownian_superiority(
        conditional_weights = 0.5e-12,
        conditional_noncentrality = 2 / 0.5e-12,
        unconditional_weights = numeric(),
        n = 2,
        tol = 1e-12
    )
    expect_identical(deterministic$probability, 1)
})


test_that("Imhof is a deterministic fallback and failures are diagnosed", {
    testthat::local_mocked_bindings(
        .quadratic_form_tail_davies = function(...) {
            list(
                probability = NA_real_, convergence = FALSE,
                diagnostic = "forced Davies failure"
            )
        },
        .package = "PastMoments"
    )
    fallback <- .variance_brownian_superiority(1, 0, 2, n = 2)
    expect_identical(fallback$method, "imhof")
    expect_true(fallback$convergence)

    testthat::local_mocked_bindings(
        .quadratic_form_tail_imhof = function(...) {
            list(
                probability = NA_real_, convergence = FALSE,
                diagnostic = "forced Imhof failure"
            )
        },
        .package = "PastMoments"
    )
    expect_warning(
        failed <- .variance_brownian_superiority(1, 0, 2, n = 2),
        "did not converge"
    )
    expect_false(failed$convergence)
    expect_true(is.na(failed$probability))
})


test_that("time-dependent Brownian comparison is invariant and optional", {
    tree <- ape::read.tree(
        text = "((a:1,b:0.6):1.2,(c:0.8,d:1.4):0.9);"
    )
    trait <- c(a = -1, b = 0.5, c = 2, d = 3)
    times <- c(0.5, 1, 1.5, 2)

    baseline <- reconstruct_empirical_moments(
        tree, trait, times = times, sigma2 = 0.8,
        quantities = "disparity_level"
    )
    shifted <- reconstruct_empirical_moments(
        tree, trait + 17, times = times, sigma2 = 0.8,
        quantities = "disparity_level"
    )
    scale <- 3.5
    rescaled <- reconstruct_empirical_moments(
        tree, scale * trait, times = times,
        sigma2 = scale^2 * 0.8,
        quantities = "disparity_level"
    )
    disabled <- reconstruct_empirical_moments(
        tree, trait, times = times, sigma2 = 0.8,
        quantities = "variance"
    )

    expect_true("disparity_level" %in% names(baseline$data))
    valid <- baseline$data$n_lineages >= 2L
    expect_true(all(
        baseline$data$disparity_level[valid] >= 0 &
            baseline$data$disparity_level[valid] <= 1
    ))
    expect_equal(
        shifted$data$disparity_level,
        baseline$data$disparity_level,
        tolerance = 2e-5
    )
    expect_equal(
        rescaled$data$disparity_level,
        baseline$data$disparity_level,
        tolerance = 2e-5
    )
    expect_false("disparity_level" %in% names(disabled$data))
    expect_equal(
        disabled$data,
        baseline$data[setdiff(
            names(baseline$data), "disparity_level"
        )]
    )
})


test_that("unconditional active covariance has the Brownian contrast geometry", {
    tree <- example_tree()
    trait <- example_trait()
    setup <- PastMoments:::.prepare_reml_setup(tree, trait, sigma2 = 1)
    parent_depth <- setup$depths_tree[tree$edge[, 1L]]
    child_depth <- setup$depths_tree[tree$edge[, 2L]]
    time <- max(setup$depths_tree)
    active_edges <- which(parent_depth < time & child_depth >= time)
    active <- PastMoments:::.compute_active_conditional(
        tree, setup, active_edges, time
    )

    n <- length(tree$tip.label)
    projector <- diag(n) - matrix(1 / n, n, n)
    original_covariance <- ape::vcv.phylo(tree)
    original_covariance <- original_covariance[
        tree$tip.label, tree$tip.label
    ]
    expect_equal(
        projector %*% active$unconditional_covariance0 %*% projector,
        projector %*% original_covariance %*% projector,
        tolerance = 1e-10
    )
})
