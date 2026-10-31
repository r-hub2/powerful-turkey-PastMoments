
## =============================================================================
## Empirical mean
## =============================================================================

.summarize_empirical_mean <- function(mu, C0, sigma2, probs) {
    n <- length(mu)
    if (n == 0L) {
        return(list(
            expectation = NA_real_,
            variance = NA_real_,
            quantiles = rep(NA_real_, length(probs))
        ))
    }

    expectation <- mean(mu)
    variance <- sigma2 * sum(C0) / (n * n)

    if (is.finite(variance) && variance < 0 && variance > -1e-10) {
        variance <- 0
    }
    if (!is.finite(variance) || variance < 0) {
        warning(sprintf("Invalid conditional variance of empirical mean: %g", variance))
        variance <- max(0, variance)
    }

    quantiles <- if (length(probs) == 0L) {
        numeric(0)
    } else {
        stats::qnorm(probs, mean = expectation, sd = sqrt(variance))
    }

    list(
        expectation = expectation,
        variance = variance,
        quantiles = quantiles
    )
}


## =============================================================================
## Empirical variance: analytic moments
## =============================================================================

.empvar_moments <- function(mu, C0, sigma2) {
    n <- length(mu)
    if (n <= 1L) {
        return(c(expectation = 0, variance = 0))
    }

    mu <- as.numeric(mu)
    C <- sigma2 * ((C0 + t(C0)) / 2)
    Pmu <- mu - mean(mu)

    # Double centering computes P C P directly and avoids subtracting large,
    # nearly equal quantities in tr(P C).
    row_mean <- rowMeans(C)
    col_mean <- colMeans(C)
    grand_mean <- mean(C)
    CPC <- sweep(C, 1L, row_mean, `-`)
    CPC <- sweep(CPC, 2L, col_mean, `-`)
    CPC <- CPC + grand_mean
    CPC <- (CPC + t(CPC)) / 2

    trace_PC <- sum(diag(CPC))
    expectation <- (sum(Pmu^2) + trace_PC) / (n - 1)

    # For Q = X' P X:
    # Var(Q) = 2 tr(P C P C) + 4 mu' P C P mu.
    # Since P C P is symmetric, tr(P C P C) = tr((P C P)^2).
    mean_covariance_term <- 4 * as.numeric(crossprod(Pmu, CPC %*% Pmu))
    var_Q <- 2 * sum(CPC^2) + mean_covariance_term
    variance <- var_Q / (n - 1)^2

    expectation_scale <- max(1, sum(Pmu^2), sum(abs(diag(CPC))))
    expectation_tol <- 1e-10 * expectation_scale
    variance_scale <- max(
        1,
        sum(CPC^2),
        sum(Pmu^2)^2,
        abs(mean_covariance_term)
    ) / (n - 1)^2
    variance_tol <- 1e-10 * variance_scale
    if (is.finite(expectation) && expectation < 0 &&
        expectation > -expectation_tol) {
        expectation <- 0
    }
    if (is.finite(variance) && variance < 0 && variance > -variance_tol) {
        variance <- 0
    }

    if (!is.finite(expectation) || expectation < 0) {
        warning(sprintf("Invalid E[S^2] from conditional covariance: %g", expectation))
        expectation <- max(0, expectation)
    }
    if (!is.finite(variance) || variance < 0) {
        warning(sprintf("Invalid Var(S^2) from conditional covariance: %g", variance))
        variance <- max(0, variance)
    }

    c(expectation = expectation, variance = variance)
}



## =============================================================================
## Empirical variance: quadratic-form decomposition and quantiles
## =============================================================================

.empvar_decomposition <- function(mu, C0, tol = 1e-12) {
    # Decompose (n - 1) S^2 as
    #
    #   shift + sum_i lambda_i (eta_i + sqrt(sigma2) Z_i)^2,
    #
    # where the decomposition below is independent of sigma2.
    #
    # The calculation is performed through an SVD in the contrast subspace.
    # This avoids treating the structural zero eigenvalue induced by the
    # centering projector as a small positive eigenvalue, which can otherwise
    # create very large spurious negative values of shift.
    n <- length(mu)
    if (n <= 1L) {
        return(list(
            n = n,
            lambda = numeric(0),
            eta = numeric(0),
            shift = 0,
            tol = tol
        ))
    }

    C0 <- (C0 + t(C0)) / 2
    Pmu <- as.numeric(mu - mean(mu))

    ec <- eigen(C0, symmetric = TRUE)
    tol_C <- tol * max(1, max(abs(ec$values)))

    warning_tol_C <- 100 * tol_C
    if (min(ec$values) < -warning_tol_C) {
        warning(sprintf(
            "Conditional covariance has a substantially negative eigenvalue: %g",
            min(ec$values)
        ))
    }

    keep_C <- ec$values > tol_C
    if (!any(keep_C)) {
        return(list(
            n = n,
            lambda = numeric(0),
            eta = numeric(0),
            shift = sum(Pmu^2),
            tol = tol
        ))
    }

    Uc <- ec$vectors[, keep_C, drop = FALSE]
    dc <- ec$values[keep_C]
    L <- sweep(Uc, 2L, sqrt(dc), `*`)

    # A = P L, where P is the centering projector. Subtracting the column
    # means is algebraically equivalent and avoids constructing P explicitly.
    A <- sweep(L, 2L, colMeans(L), `-`)
    sA <- svd(A, nu = min(dim(A)), nv = 0L)

    # Since lambda_i = d_i^2, apply the relative tolerance on squared singular
    # values. This reliably removes the structural null direction of P.
    tol_d <- sqrt(tol * max(1, max(sA$d)^2))
    keep_d <- sA$d > tol_d

    if (!any(keep_d)) {
        return(list(
            n = n,
            lambda = numeric(0),
            eta = numeric(0),
            shift = sum(Pmu^2),
            tol = tol
        ))
    }

    d <- sA$d[keep_d]
    U <- sA$u[, keep_d, drop = FALSE]
    projected_mean <- as.numeric(crossprod(U, Pmu))

    lambda <- d^2
    eta <- projected_mean / d
    shift <- sum(Pmu^2) - sum(projected_mean^2)
    shift_tol <- 100 * tol * max(1, sum(Pmu^2))

    if (is.finite(shift) && shift < 0 && shift > -shift_tol) {
        shift <- 0
    }
    if (!is.finite(shift) || shift < -shift_tol) {
        warning(sprintf(
            "Invalid deterministic shift in empirical variance after SVD: %g",
            shift
        ))
        shift <- max(0, shift)
    }

    list(n = n, lambda = lambda, eta = eta, shift = shift, tol = tol)
}


## =============================================================================
## Empirical standard deviation: fractional moment
## =============================================================================

.empvar_expected_sd <- function(decomposition, sigma2,
                                rel.tol = 1e-8, abs.tol = 1e-10,
                                subdivisions = 200L) {
    n <- decomposition$n
    lambda <- decomposition$lambda
    eta <- decomposition$eta
    shift <- decomposition$shift

    if (n <= 1L) {
        return(0)
    }

    # The internal convention is Q = (n - 1) S^2:
    #
    #   Q = shift + sum_i lambda_i (eta_i + sqrt(sigma2) Z_i)^2.
    #
    # Clamp only negative eigenvalues within the decomposition's spectral
    # tolerance.  Substantially negative values indicate an invalid quadratic
    # form and must not silently disappear from the fractional moment.
    spectral_tol <- decomposition$tol
    if (is.null(spectral_tol)) {
        spectral_tol <- 1e-12
    }
    lambda_scale <- max(1, abs(lambda), na.rm = TRUE)
    lambda_tol <- spectral_tol * lambda_scale
    if (any(!is.finite(lambda)) || any(!is.finite(eta))) {
        stop("Non-finite spectral parameter in empirical-standard-deviation integration.")
    }
    if (any(lambda < -lambda_tol)) {
        stop(sprintf(
            "Empirical-variance quadratic form has a materially negative eigenvalue: %g",
            min(lambda)
        ))
    }
    lambda[lambda < 0] <- 0
    # Positive eigenvalues have already been screened by the decomposition;
    # retain them here so deliberately wide valid spectra are not truncated a
    # second time relative to their largest member.
    keep <- lambda > 0
    lambda <- lambda[keep]
    eta <- eta[keep]

    # In the standard noncentral-chi-square notation, the stochastic weights
    # are lambda_i * sigma2 and b_i = eta_i / sqrt(sigma2).  Retaining
    # lambda_i * eta_i^2 avoids forming b_i when sigma2 is zero or tiny.
    c_value <- shift
    stochastic_weight <- lambda * sigma2
    noncentral_weight <- lambda * eta^2
    deterministic_value <- c_value + sum(noncentral_weight)

    if (length(lambda) == 0L || sigma2 == 0 ||
        !any(stochastic_weight > 0)) {
        return(sqrt(max(0, deterministic_value)) / sqrt(n - 1L))
    }

    representation_mean <- deterministic_value + sum(stochastic_weight)
    if (!is.finite(representation_mean) || representation_mean < 0) {
        stop("Invalid mean in empirical-standard-deviation integration.")
    }
    if (representation_mean == 0) {
        return(0)
    }

    # Scale Q to unit expectation.  This makes the quadrature tolerances
    # meaningful for forms whose eigenvalues are very small or very large.
    c_scaled <- c_value / representation_mean
    w_scaled <- stochastic_weight / representation_mean
    b_scaled <- noncentral_weight / representation_mean

    # Use E[sqrt(Q)] = integral t^(-1/2) {-L'(t)} dt / sqrt(pi), where
    # -L'(t) = L(t) H(t).  With t = x^2 and x = u / (1-u), the finite-interval
    # integrand below is L(t) H(t) / (1-u)^2.
    integrand_one <- function(u) {
        if (u <= 0) {
            # H(0) = E[Q / E(Q)] = 1 after scaling.
            return(1)
        }
        if (u >= 1) {
            return(0)
        }

        one_minus_u <- 1 - u
        x <- u / one_minus_u
        if (!is.finite(x) || x > sqrt(.Machine$double.xmax)) {
            return(0)
        }
        t <- x * x
        denominator <- 1 + 2 * w_scaled * t
        log_L <- -c_scaled * t - 0.5 * sum(log1p(2 * w_scaled * t)) -
            sum(b_scaled * t / denominator)
        H <- c_scaled + sum(
            w_scaled / denominator + b_scaled / denominator^2
        )
        if (!is.finite(H) || H <= 0) {
            return(0)
        }
        log_integrand <- log_L + log(H) - 2 * log(one_minus_u)
        exp(log_integrand)
    }
    integrand <- function(u) vapply(u, integrand_one, numeric(1))

    quadrature <- stats::integrate(
        integrand,
        lower = 0,
        upper = 1,
        rel.tol = rel.tol,
        abs.tol = abs.tol,
        subdivisions = subdivisions,
        stop.on.error = FALSE
    )
    if (!identical(quadrature$message, "OK") ||
        !is.finite(quadrature$value)) {
        stop(
            "Expected empirical standard deviation integration failed: ",
            quadrature$message,
            call. = FALSE
        )
    }

    2 * sqrt(representation_mean) * quadrature$value /
        sqrt(pi * (n - 1L))
}

.empvar_quantile_rank1 <- function(target_survival, lambda, eta, sigma2) {
    s <- sqrt(sigma2)

    survival_q <- function(q) {
        if (q <= 0) {
            return(1)
        }
        r <- sqrt(q)
        stats::pnorm((-r - eta) / s) +
            stats::pnorm((r - eta) / s, lower.tail = FALSE)
    }

    f <- function(q) survival_q(q) - target_survival
    mean_q <- eta^2 + sigma2
    var_q <- 2 * (sigma2^2 + 2 * sigma2 * eta^2)
    upper <- max(mean_q + 10 * sqrt(var_q), .Machine$double.eps)
    f_upper <- f(upper)

    while (!is.finite(f_upper) || f_upper > 0) {
        upper <- 2 * upper
        f_upper <- f(upper)
        if (!is.finite(upper) || upper > 1e100) {
            stop("Could not bracket rank-one empirical-variance quantile.")
        }
    }

    lambda * stats::uniroot(
        f,
        lower = 0,
        upper = upper,
        f.lower = 1 - target_survival,
        f.upper = f_upper,
        tol = 1e-14
    )$root
}


.empvar_quantile_qapprox <- function(target_survival, lambda, eta, sigma2,
                                     method = "SW") {
    if (!requireNamespace("Qapprox", quietly = TRUE)) {
        stop(
            "Package 'Qapprox' is required for quantile_method = 'qapprox'.",
            call. = FALSE
        )
    }

    keep <- is.finite(lambda) & is.finite(eta) & lambda > 0
    lambda <- lambda[keep]
    eta <- eta[keep]
    m <- length(lambda)

    if (m == 0L) {
        return(0)
    }
    if (m == 1L) {
        return(.empvar_quantile_rank1(
            target_survival, lambda[[1L]], eta[[1L]], sigma2
        ))
    }

    Sigma <- diag(sigma2, m)
    mean_q <- sum(lambda * (eta^2 + sigma2))
    var_q <- 2 * sum(lambda^2 * (sigma2^2 + 2 * sigma2 * eta^2))
    scale <- mean_q

    if (!is.finite(scale) || scale <= 0) {
        stop("Invalid scale in Qapprox quantile computation.")
    }

    A_scaled <- diag(lambda / scale, m)

    f <- function(z) {
        Qapprox::Qapprox_nc(
            q = z,
            mu = eta,
            Sigma = Sigma,
            A = A_scaled,
            method = method
        ) - target_survival
    }

    upper <- max(
        (mean_q + 10 * sqrt(var_q)) / scale,
        .Machine$double.eps
    )
    f_upper <- f(upper)

    while (!is.finite(f_upper) || f_upper > 0) {
        upper <- 2 * upper
        f_upper <- f(upper)
        if (!is.finite(upper) || upper > 1e100) {
            stop("Could not bracket Qapprox empirical-variance quantile.")
        }
    }

    z <- stats::uniroot(
        f,
        lower = 0,
        upper = upper,
        f.lower = 1 - target_survival,
        f.upper = f_upper,
        tol = 1e-12
    )$root

    scale * z
}


.empvar_quantile_davies <- function(target_survival, lambda, eta, sigma2) {
    if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
        stop(
            "Package 'CompQuadForm' is required for quantile_method = 'davies'.",
            call. = FALSE
        )
    }

    if (length(lambda) == 0L) {
        return(0)
    }
    if (length(lambda) == 1L) {
        return(.empvar_quantile_rank1(
            target_survival, lambda[[1L]], eta[[1L]], sigma2
        ))
    }

    delta2 <- eta^2 / sigma2
    f <- function(x) {
        if (x <= 0) {
            return(1 - target_survival)
        }
        CompQuadForm::davies(
            x / sigma2,
            lambda = lambda,
            delta = delta2
        )$Qq - target_survival
    }

    mean_q <- sum(lambda * (eta^2 + sigma2))
    var_q <- 2 * sum(lambda^2 * (sigma2^2 + 2 * sigma2 * eta^2))
    upper <- max(mean_q + 10 * sqrt(var_q), .Machine$double.eps)
    f_upper <- f(upper)

    while (!is.finite(f_upper) || f_upper > 0) {
        upper <- 2 * upper
        f_upper <- f(upper)
        if (!is.finite(upper) || upper > 1e100) {
            stop("Could not bracket Davies empirical-variance quantile.")
        }
    }

    stats::uniroot(
        f,
        lower = 0,
        upper = upper,
        f.lower = 1 - target_survival,
        f.upper = f_upper
    )$root
}


.empvar_quantile_liu <- function(target_survival, lambda, eta, sigma2) {
    if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
        stop(
            "Package 'CompQuadForm' is required for quantile_method = 'liu'.",
            call. = FALSE
        )
    }

    if (length(lambda) == 0L) {
        return(0)
    }
    if (length(lambda) == 1L) {
        return(.empvar_quantile_rank1(
            target_survival, lambda[[1L]], eta[[1L]], sigma2
        ))
    }

    delta2 <- eta^2 / sigma2
    survival <- function(x) {
        as.numeric(CompQuadForm::liu(
            x / sigma2,
            lambda = lambda,
            delta = delta2
        ))
    }
    f <- function(x) survival(x) - target_survival

    mean_q <- sum(lambda * (eta^2 + sigma2))
    var_q <- 2 * sum(lambda^2 * (sigma2^2 + 2 * sigma2 * eta^2))
    upper <- max(mean_q + 10 * sqrt(var_q), .Machine$double.eps)
    f_upper <- f(upper)

    while (!is.finite(f_upper) || f_upper > 0) {
        upper <- 2 * upper
        f_upper <- f(upper)
        if (!is.finite(upper) || upper > 1e100) {
            stop("Could not bracket Liu empirical-variance quantile.")
        }
    }

    stats::uniroot(f, lower = 0, upper = upper)$root
}


.empvar_quantiles <- function(mu, C0, sigma2, probs,
                              method = c("qapprox", "davies", "liu"),
                              tol = 1e-12, decomposition = NULL) {
    method <- match.arg(method)
    n <- length(mu)

    if (length(probs) == 0L) {
        return(numeric(0))
    }
    if (n <= 1L) {
        return(rep(0, length(probs)))
    }

    # When sigma2 is zero, the conditional distribution is deterministic.
    if (sigma2 == 0) {
        value <- sum((mu - mean(mu))^2) / (n - 1)
        return(rep(value, length(probs)))
    }

    if (is.null(decomposition)) {
        decomposition <- .empvar_decomposition(mu, C0, tol = tol)
    }

    if (length(decomposition$lambda) == 0L) {
        value <- decomposition$shift / (n - 1)
        return(rep(value, length(probs)))
    }

    stochastic_quantile <- function(prob) {
        target_survival <- 1 - prob

        switch(
            method,
            qapprox = .empvar_quantile_qapprox(
                target_survival,
                decomposition$lambda,
                decomposition$eta,
                sigma2
            ),
            davies = .empvar_quantile_davies(
                target_survival,
                decomposition$lambda,
                decomposition$eta,
                sigma2
            ),
            liu = .empvar_quantile_liu(
                target_survival,
                decomposition$lambda,
                decomposition$eta,
                sigma2
            )
        )
    }

    q_stochastic <- vapply(probs, stochastic_quantile, numeric(1))
    (decomposition$shift + q_stochastic) / (n - 1)
}


.summarize_empirical_variance <- function(mu, C0, sigma2, probs,
                                          quantile_method, tol = 1e-12) {
    moments <- .empvar_moments(mu, C0, sigma2)
    decomposition <- .empvar_decomposition(mu, C0, tol = tol)
    expected_sd <- .empvar_expected_sd(decomposition, sigma2)
    quantiles <- .empvar_quantiles(
        mu, C0, sigma2, probs,
        method = quantile_method,
        tol = tol,
        decomposition = decomposition
    )

    expected_variance <- unname(moments["expectation"])

    list(
        expectation = expected_variance,
        variance = unname(moments["variance"]),
        expected_variance = expected_variance,
        expected_sd = expected_sd,
        quantiles = quantiles,
        decomposition = decomposition
    )
}
