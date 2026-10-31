expected_sd_from_form <- function(lambda = numeric(0),
                                  eta = numeric(length(lambda)),
                                  shift = 0, sigma2 = 1) {
    decomposition <- list(
        n = 2L,
        lambda = lambda,
        eta = eta,
        shift = shift
    )
    PastMoments:::.empvar_expected_sd(decomposition, sigma2)
}


normal_quadrature <- function(order = 80L) {
    J <- matrix(0, order, order)
    off_diagonal <- sqrt(seq_len(order - 1L))
    J[cbind(seq_len(order - 1L), 2:order)] <- off_diagonal
    J[cbind(2:order, seq_len(order - 1L))] <- off_diagonal
    eig <- eigen(J, symmetric = TRUE)
    list(nodes = eig$values, weights = eig$vectors[1L, ]^2)
}


independent_normal_reference <- function(lambda, eta, shift = 0,
                                         sigma2 = 1, order = 80L) {
    gh <- normal_quadrature(order)
    grid <- expand.grid(seq_along(gh$nodes), seq_along(gh$nodes))
    z1 <- gh$nodes[grid[[1L]]]
    z2 <- gh$nodes[grid[[2L]]]
    values <- sqrt(
        shift + lambda[1L] * (eta[1L] + sqrt(sigma2) * z1)^2 +
            lambda[2L] * (eta[2L] + sqrt(sigma2) * z2)^2
    )
    sum(gh$weights[grid[[1L]]] * gh$weights[grid[[2L]]] * values)
}


previous_fractional_moment_reference <- function(lambda, b, shift = 0) {
    mean_q <- shift + sum(lambda * (1 + b^2))
    if (mean_q == 0) {
        return(0)
    }
    integrand <- function(u) {
        vapply(u, function(ui) {
            if (ui <= 0) {
                return(mean_q)
            }
            if (ui >= 1) {
                return(1)
            }
            t <- (ui / (1 - ui))^2
            denominator <- 1 + 2 * lambda * t
            log_laplace <- -shift * t - 0.5 * sum(log1p(2 * lambda * t)) -
                sum(lambda * b^2 * t / denominator)
            -expm1(log_laplace) / ui^2
        }, numeric(1))
    }
    stats::integrate(
        integrand, 0, 1, rel.tol = 1e-11, abs.tol = 1e-12,
        subdivisions = 1000L
    )$value / sqrt(pi)
}


test_that("expected SD handles positive and zero deterministic forms", {
    expect_identical(expected_sd_from_form(shift = 9), 3)
    expect_identical(expected_sd_from_form(), 0)

    summary <- PastMoments:::.summarize_empirical_variance(
        mu = c(1, 1), C0 = matrix(0, 2, 2), sigma2 = 1,
        probs = NULL, quantile_method = "qapprox"
    )
    expect_identical(summary$expected_sd, 0)
    expect_false(any(c(
        "sqrt_expected_variance", "jensen_gap", "jensen_ratio"
    ) %in% names(summary)))
})


test_that("central scaled chi-square forms agree with the exact moment", {
    cases <- expand.grid(nu = c(1L, 2L, 5L, 20L),
                         lambda = c(0.03, 1, 12))
    for (i in seq_len(nrow(cases))) {
        nu <- cases$nu[i]
        scale <- cases$lambda[i]
        exact <- sqrt(2 * scale) * exp(
            lgamma((nu + 1) / 2) - lgamma(nu / 2)
        )
        actual <- expected_sd_from_form(rep(scale, nu), rep(0, nu))
        expect_equal(actual, exact, tolerance = 2e-8)
    }
})


test_that("rank-one central and noncentral forms have exact Gaussian moments", {
    lambda <- 0.37
    expect_equal(
        expected_sd_from_form(lambda),
        sqrt(lambda) * sqrt(2 / pi),
        tolerance = 2e-8
    )

    for (b in c(-2.3, 0.4, 3)) {
        folded_normal_mean <- 2 * stats::dnorm(b) +
            b * (2 * stats::pnorm(b) - 1)
        expect_equal(
            expected_sd_from_form(lambda, eta = b),
            sqrt(lambda) * folded_normal_mean,
            tolerance = 2e-8
        )
    }
})


test_that("derivative integral agrees with the previous fractional identity", {
    set.seed(26071986)
    for (i in seq_len(8L)) {
        rank <- sample.int(5L, 1L)
        lambda <- exp(stats::runif(rank, log(1e-4), log(20)))
        b <- stats::rnorm(rank, sd = 1.5)
        shift <- sample(c(0, stats::runif(1L, 0.01, 3)), 1L)
        reference <- previous_fractional_moment_reference(lambda, b, shift)
        expect_equal(
            expected_sd_from_form(lambda, b, shift),
            reference,
            tolerance = 3e-8
        )
    }
})


test_that("Q convention is normalized exactly once for general sample size", {
    n <- 6L
    lambda <- c(0.2, 1.4, 3)
    eta <- c(-0.3, 1.1, 0)
    shift <- 0.7
    decomposition <- list(
        n = n, lambda = lambda, eta = eta, shift = shift, tol = 1e-12
    )
    expected_sqrt_q <- previous_fractional_moment_reference(
        lambda, eta, shift
    )
    expect_equal(
        PastMoments:::.empvar_expected_sd(decomposition, sigma2 = 1),
        expected_sqrt_q / sqrt(n - 1),
        tolerance = 3e-8
    )
})


test_that("spectral tolerance clamps roundoff but rejects negative eigenvalues", {
    negligible <- list(
        n = 2L, lambda = c(-5e-13, 1), eta = c(10, 0), shift = 0,
        tol = 1e-12
    )
    expect_equal(
        PastMoments:::.empvar_expected_sd(negligible, sigma2 = 1),
        sqrt(2 / pi),
        tolerance = 2e-8
    )

    material <- negligible
    material$lambda[1L] <- -1e-4
    expect_error(
        PastMoments:::.empvar_expected_sd(material, sigma2 = 1),
        "materially negative eigenvalue"
    )
})


test_that("Jensen inequality holds for central and noncentral forms", {
    forms <- list(
        list(lambda = c(0.2, 1.7), eta = c(0, 0), shift = 0),
        list(lambda = c(0.1, 0.4, 2), eta = c(-2, 0.5, 1), shift = 0.3),
        list(lambda = c(1e-6, 1, 1e6), eta = c(3, -1, 0.2), shift = 2),
        list(lambda = numeric(0), eta = numeric(0), shift = 4)
    )
    for (form in forms) {
        actual <- expected_sd_from_form(
            form$lambda, form$eta, form$shift
        )
        root_mean <- sqrt(
            form$shift + sum(form$lambda * (1 + form$eta^2))
        )
        expect_lte(actual, root_mean + 1e-8 * max(1, root_mean))
    }
    expect_identical(expected_sd_from_form(shift = 4), 2)
})


test_that("unequal central weights match an independent polar integral", {
    for (lambda in list(c(0.2, 1.7), c(1, 8), c(0.01, 3))) {
        angular_mean <- stats::integrate(
            function(theta) sqrt(
                lambda[1L] * cos(theta)^2 +
                    lambda[2L] * sin(theta)^2
            ),
            0, 2 * pi,
            rel.tol = 1e-11
        )$value / (2 * pi)
        reference <- sqrt(pi / 2) * angular_mean
        expect_equal(
            expected_sd_from_form(lambda), reference, tolerance = 2e-8
        )
    }
})


test_that("noncentral forms match independent Gaussian quadrature", {
    forms <- list(
        list(lambda = c(0.4, 1.6), eta = c(0.7, -1.2), shift = 0.3),
        list(lambda = c(0.05, 3), eta = c(2.5, 0.4), shift = 1)
    )
    for (form in forms) {
        reference <- independent_normal_reference(
            form$lambda, form$eta, form$shift, order = 100L
        )
        expect_equal(
            expected_sd_from_form(form$lambda, form$eta, form$shift),
            reference,
            tolerance = 2e-7
        )
    }
})


test_that("zero eigenvalues, singular forms, and wide scales are stable", {
    expect_equal(
        expected_sd_from_form(c(0, 1, 1e-30), c(100, 0, 1e10)),
        sqrt(2 / pi),
        tolerance = 2e-8
    )
    decomposition <- PastMoments:::.empvar_decomposition(
        mu = c(0, 0, 0),
        C0 = diag(c(1, 1e-30, 0)),
        tol = 1e-12
    )
    expect_length(decomposition$lambda, 1L)

    wide <- expected_sd_from_form(
        lambda = c(1e-12, 1, 1e12),
        eta = c(1, -0.5, 2),
        shift = 1e-8
    )
    expect_true(is.finite(wide))
    expect_gt(wide, 0)
    expect_lte(wide, sqrt(1e-8 + sum(
        c(1e-12, 1, 1e12) * (1 + c(1, -0.5, 2)^2)
    )))
})
