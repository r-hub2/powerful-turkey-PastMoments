test_that("Brownian reconstruction recovers empirical moments at the tips", {
    tree <- example_tree()
    trait <- example_trait()

    result <- reconstruct_empirical_moments(
        tree,
        trait,
        times = 2,
        sigma2 = 1,
        mean_probs = NULL,
        variance_probs = NULL
    )

    expect_s3_class(result, "empirical_moment_reconstruction")
    expect_s3_class(result$tree, "phylo")
    expect_equal(result$trait, trait[result$tree$tip.label])
    expect_identical(result$clades, list(all = result$tree$tip.label))
    expect_identical(result$model$name, "Brownian")
    expect_equal(result$model$sigma2, 1)

    mean_row <- result$data[result$data$quantity == "empirical_mean", ]
    variance_row <- result$data[result$data$quantity == "empirical_variance", ]

    expect_equal(mean_row$n_lineages, 4L)
    expect_equal(mean_row$expectation, mean(trait), tolerance = 1e-8)
    expect_equal(mean_row$variance, 0, tolerance = 1e-8)
    expect_equal(variance_row$expectation, stats::var(trait), tolerance = 1e-8)
    expect_equal(variance_row$variance, 0, tolerance = 1e-8)
    expect_equal(variance_row$expected_variance, stats::var(trait),
                 tolerance = 1e-8)
    expect_equal(variance_row$expected_sd, sqrt(stats::var(trait)),
                 tolerance = 1e-8)
    expect_false(any(c(
        "sqrt_expected_variance", "jensen_ratio", "jensen_gap"
    ) %in% names(variance_row)))
})


test_that("combined Brownian interface agrees with separate interfaces", {
    tree <- example_tree()
    trait <- example_trait()
    times <- c(0.5, 1.5, 2)

    both <- reconstruct_empirical_moments(
        tree, trait, times = times, sigma2 = 1,
        mean_probs = NULL, variance_probs = NULL
    )
    mean_only <- reconstruct_empirical_moments(
        tree, trait, times = times, sigma2 = 1,
        quantities = "mean", mean_probs = NULL
    )
    variance_only <- reconstruct_empirical_moments(
        tree, trait, times = times, sigma2 = 1,
        quantities = "variance", variance_probs = NULL
    )

    expect_equal(
        unname(both$data[
            both$data$quantity == "empirical_mean",
            names(mean_only$data)
        ]),
        unname(mean_only$data),
        ignore_attr = TRUE
    )
    expect_equal(
        unname(both$data[
            both$data$quantity == "empirical_variance",
            names(variance_only$data)
        ]),
        unname(variance_only$data),
        ignore_attr = TRUE
    )
})


test_that("reconstruction objects expose their table explicitly", {
    result <- reconstruct_empirical_moments(
        example_tree(), example_trait(), times = 1, sigma2 = 1,
        quantities = "mean"
    )

    expect_equal(as.data.frame(result), result$data)
    expect_output(print(result), "empirical_moment_reconstruction")
})


test_that("combined reconstructions expose disparity levels", {
    result <- reconstruct_empirical_moments(
        example_tree(), example_trait(), times = c(1, 2), sigma2 = 1,
        mean_probs = NULL, variance_probs = NULL,
        quantities = c("mean", "variance", "disparity_level")
    )

    expect_true("disparity_level" %in% names(result$data))
    expect_true(all(is.na(
        result$data$disparity_level[
            result$data$quantity == "empirical_mean"
        ]
    )))
    expect_true(all(is.finite(
        result$data$disparity_level[
            result$data$quantity == "empirical_variance"
        ]
    )))
    expect_output(print(result), "disparity level: data\\$disparity_level")
})


test_that("quantities defines the single reconstruction API", {
    api <- formals(reconstruct_empirical_moments)
    expect_true("quantities" %in% names(api))
    expect_false("disparity_level" %in% names(api))
    expect_false(any(c(
        "reconstruct_empirical_mean", "reconstruct_empirical_variance"
    ) %in% getNamespaceExports("PastMoments")))

    defaults <- reconstruct_empirical_moments(
        example_tree(), example_trait(), times = 1, sigma2 = 1,
        mean_probs = NULL, variance_probs = NULL
    )
    expect_identical(
        defaults$settings$quantities,
        c("mean", "variance", "disparity_level")
    )
    expect_true("disparity_level" %in% names(defaults$data))

    disparity_only <- reconstruct_empirical_moments(
        example_tree(), example_trait(), times = 1, sigma2 = 1,
        quantities = "disparity_level",
        mean_probs = NULL, variance_probs = NULL
    )
    expect_identical(unique(disparity_only$data$quantity),
                     "empirical_variance")
    expect_identical(
        disparity_only$settings$quantities,
        c("variance", "disparity_level")
    )

    expect_error(
        reconstruct_empirical_moments(
            example_tree(), example_trait(), quantities = character()
        ),
        "quantities must select"
    )
    expect_warning(
        supported <- reconstruct_empirical_moments(
            example_tree(), example_trait(), times = 1, sigma2 = 1,
            quantities = c("mean", "expected_sd"), mean_probs = NULL
        ),
        "'expected_sd'"
    )
    expect_identical(unique(supported$data$quantity), "empirical_mean")
    expect_warning(
        expect_error(
            reconstruct_empirical_moments(
                example_tree(), example_trait(), quantities = "median"
            ),
            "No supported quantity remains"
        ),
        "'median'"
    )
})


test_that("negligible negative covariance residuals are treated as zero", {
    mu <- c(-1, 1)
    tiny_negative <- -1e-14 * diag(2)

    expect_warning(
        moments <- PastMoments:::.empvar_moments(mu, tiny_negative, sigma2 = 1),
        NA
    )
    expect_equal(unname(moments["variance"]), 0)
    expect_warning(
        PastMoments:::.empvar_decomposition(
            mu, diag(c(-2e-12, 1)), tol = 1e-12
        ),
        NA
    )

    expect_warning(
        PastMoments:::.empvar_moments(mu, -0.01 * diag(2), sigma2 = 1),
        "Invalid Var"
    )
    expect_warning(
        PastMoments:::.empvar_decomposition(mu, diag(c(-0.01, 1))),
        "negative eigenvalue"
    )
})
