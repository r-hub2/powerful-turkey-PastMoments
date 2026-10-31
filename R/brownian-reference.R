## =============================================================================
## Brownian reference comparison for empirical variance
## =============================================================================

.clean_empvar_spectrum <- function(weights, noncentrality, tol, label) {
    weights <- as.numeric(weights)
    noncentrality <- as.numeric(noncentrality)
    if (length(noncentrality) == 0L && length(weights) > 0L) {
        noncentrality <- rep(0, length(weights))
    }
    if (length(weights) != length(noncentrality) ||
        any(!is.finite(weights)) || any(!is.finite(noncentrality)) ||
        any(noncentrality < 0)) {
        stop(label, " weights and noncentralities must be finite, aligned, ",
             "and have nonnegative noncentralities.", call. = FALSE)
    }
    scale <- max(1, abs(weights))
    threshold <- tol * scale
    if (any(weights < -threshold)) {
        stop(label, " spectrum contains a materially negative eigenvalue.",
             call. = FALSE)
    }
    weights[weights < 0] <- 0
    keep <- weights > threshold
    discarded <- !keep
    # For w * chi-square_1(delta), both the mean w * (1 + delta) and
    # variance 2 * w^2 * (1 + 2 * delta) vanish with w when delta is
    # bounded. If delta grows like 1 / w, however, the random term still
    # converges to the finite deterministic limit w * delta. Preserve that
    # limit instead of dropping the entire component at the spectral cutoff.
    discarded_shift <- sum(
        weights[discarded] * (1 + noncentrality[discarded])
    )
    if (!is.finite(discarded_shift)) {
        stop(
            label,
            " negligible spectral components have a non-finite deterministic ",
            "limit.",
            call. = FALSE
        )
    }
    list(
        weights = weights[keep],
        noncentrality = noncentrality[keep],
        discarded_shift = discarded_shift
    )
}


.quadratic_form_tail_davies <- function(q, weights, df, ncp, accuracy) {
    if (length(weights) == 1L && weights[[1L]] > 0) {
        probability <- if (q <= 0) {
            1
        } else {
            stats::pchisq(
                q / weights[[1L]],
                df = df[[1L]],
                ncp = ncp[[1L]],
                lower.tail = FALSE
            )
        }
        return(list(
            probability = as.numeric(probability),
            convergence = is.finite(probability),
            diagnostic = "analytic scaled chi-squared tail"
        ))
    }
    answer <- suppressWarnings(
        CompQuadForm::davies(
            q = q, lambda = weights, h = df, delta = ncp,
            lim = 10000L, acc = accuracy
        )
    )
    list(
        probability = as.numeric(answer$Qq),
        convergence = identical(as.integer(answer$ifault), 0L),
        diagnostic = paste0("ifault=", answer$ifault)
    )
}


.quadratic_form_tail_imhof <- function(q, weights, df, ncp, accuracy) {
    answer <- CompQuadForm::imhof(
        q = q, lambda = weights, h = df, delta = ncp,
        epsabs = accuracy, epsrel = accuracy, limit = 10000L
    )
    list(
        probability = as.numeric(answer$Qq),
        convergence = is.finite(answer$Qq) && is.finite(answer$abserr) &&
            answer$abserr <= max(10 * accuracy, 1e-5),
        diagnostic = paste0("abserr=", format(answer$abserr, digits = 4L))
    )
}


.indefinite_quadratic_probability <- function(
    conditional_weights,
    conditional_noncentrality,
    unconditional_weights,
    conditional_shift = 0,
    unconditional_shift = 0,
    tol = 1e-12,
    accuracy = 1e-6,
    warn = TRUE
) {
    if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
        stop(
            "Package 'CompQuadForm' is required for the exact Brownian ",
            "variance comparison.",
            call. = FALSE
        )
    }
    scalars <- c(
        conditional_shift = conditional_shift,
        unconditional_shift = unconditional_shift,
        tol = tol,
        accuracy = accuracy
    )
    if (any(!is.finite(scalars)) || tol <= 0 || accuracy <= 0) {
        stop("Shifts must be finite and tolerances strictly positive.",
             call. = FALSE)
    }

    conditional <- .clean_empvar_spectrum(
        conditional_weights, conditional_noncentrality, tol, "Conditional"
    )
    unconditional <- .clean_empvar_spectrum(
        unconditional_weights,
        rep(0, length(unconditional_weights)),
        tol,
        "Unconditional"
    )
    weights <- c(conditional$weights, -unconditional$weights)
    ncp <- c(
        conditional$noncentrality,
        rep(0, length(unconditional$weights))
    )
    df <- rep(1, length(weights))
    conditional_total_shift <- conditional_shift +
        conditional$discarded_shift
    unconditional_total_shift <- unconditional_shift +
        unconditional$discarded_shift
    shift <- conditional_total_shift - unconditional_total_shift
    scale <- max(1, abs(shift), abs(weights))
    zero_tol <- tol * scale

    if (length(weights) == 0L) {
        both_zero <- abs(conditional_total_shift) <= zero_tol &&
            abs(unconditional_total_shift) <= zero_tol
        probability <- if (both_zero) {
            NA_real_
        } else if (shift > zero_tol) {
            1
        } else {
            0
        }
        return(list(
            probability = probability,
            weights = weights,
            df = df,
            ncp = ncp,
            shift = shift,
            method = "degenerate",
            convergence = is.finite(probability),
            diagnostic = if (both_zero) {
                "both distributions are degenerate at zero"
            } else {
                "deterministic comparison"
            }
        ))
    }

    # CompQuadForm::davies() and imhof() return Qq = Pr(form > q), not the
    # lower tail. Q = shift + form, hence Pr(Q > 0) = Pr(form > -shift).
    q <- -shift
    reverse_tail <- length(weights) > 0L && all(weights < 0)
    evaluation_weights <- if (reverse_tail) -weights else weights
    evaluation_q <- if (reverse_tail) -q else q
    orient_probability <- function(probability) {
        if (reverse_tail && is.finite(probability)) {
            1 - probability
        } else {
            probability
        }
    }
    primary <- tryCatch(
        .quadratic_form_tail_davies(
            evaluation_q, evaluation_weights, df, ncp, accuracy
        ),
        error = function(e) list(
            probability = NA_real_, convergence = FALSE,
            diagnostic = conditionMessage(e)
        )
    )
    primary$probability <- orient_probability(primary$probability)
    result <- primary
    method <- "davies"

    valid_probability <- function(z) {
        is.finite(z) && z >= -10 * accuracy && z <= 1 + 10 * accuracy
    }
    if (!isTRUE(primary$convergence) ||
        !valid_probability(primary$probability)) {
        fallback <- tryCatch(
            .quadratic_form_tail_imhof(
                evaluation_q, evaluation_weights, df, ncp, accuracy
            ),
            error = function(e) list(
                probability = NA_real_, convergence = FALSE,
                diagnostic = conditionMessage(e)
            )
        )
        fallback$probability <- orient_probability(fallback$probability)
        result <- fallback
        method <- "imhof"
        result$diagnostic <- paste0(
            "Davies: ", primary$diagnostic, "; Imhof: ",
            fallback$diagnostic
        )
    }

    if (!isTRUE(result$convergence) ||
        !valid_probability(result$probability)) {
        if (isTRUE(warn)) {
            warning(
                "Brownian variance comparison did not converge reliably (",
                result$diagnostic, ").",
                call. = FALSE
            )
        }
        probability <- NA_real_
        convergence <- FALSE
    } else {
        probability <- min(1, max(0, result$probability))
        convergence <- TRUE
    }

    list(
        probability = probability,
        weights = weights,
        df = df,
        ncp = ncp,
        shift = shift,
        method = method,
        convergence = convergence,
        diagnostic = result$diagnostic
    )
}


.empvar_spectral_representation <- function(decomposition, sigma2) {
    lambda <- as.numeric(decomposition$lambda)
    eta <- as.numeric(decomposition$eta)
    if (length(lambda) != length(eta)) {
        stop("Invalid empirical-variance spectral decomposition.", call. = FALSE)
    }
    if (sigma2 > 0) {
        list(
            weights = lambda * sigma2,
            noncentrality = eta^2 / sigma2,
            shift = decomposition$shift
        )
    } else {
        list(
            weights = numeric(0),
            noncentrality = numeric(0),
            shift = decomposition$shift + sum(lambda * eta^2)
        )
    }
}


## Internal convenience wrapper around the indefinite quadratic-form
## probability calculation. It accepts the two spectral representations of
## (n - 1) S^2 and returns cleaned spectra together with numerical diagnostics.
.variance_brownian_superiority <- function(
    conditional_weights,
    conditional_noncentrality = rep(0, length(conditional_weights)),
    unconditional_weights,
    conditional_shift = 0,
    unconditional_shift = 0,
    n = max(
        length(conditional_weights),
        length(unconditional_weights)
    ) + 1L,
    tol = 1e-12,
    accuracy = 1e-6
) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0) {
        stop("n must be a single nonnegative finite number.", call. = FALSE)
    }
    if (n < 2) {
        return(list(
            probability = NA_real_,
            conditional_weights = numeric(0),
            conditional_noncentrality = numeric(0),
            unconditional_weights = numeric(0),
            method = "degenerate",
            convergence = FALSE,
            diagnostic = "fewer than two points"
        ))
    }

    answer <- .indefinite_quadratic_probability(
        conditional_weights = conditional_weights,
        conditional_noncentrality = conditional_noncentrality,
        unconditional_weights = unconditional_weights,
        conditional_shift = conditional_shift,
        unconditional_shift = unconditional_shift,
        tol = tol,
        accuracy = accuracy
    )
    n_conditional <- sum(answer$weights > 0)
    conditional_index <- if (n_conditional > 0L) {
        seq_len(n_conditional)
    } else {
        integer(0)
    }
    unconditional_index <- if (n_conditional < length(answer$weights)) {
        seq.int(n_conditional + 1L, length(answer$weights))
    } else {
        integer(0)
    }
    list(
        probability = answer$probability,
        conditional_weights = answer$weights[conditional_index],
        conditional_noncentrality = answer$ncp[conditional_index],
        unconditional_weights = -answer$weights[unconditional_index],
        conditional_shift = conditional_shift,
        unconditional_shift = unconditional_shift,
        method = answer$method,
        convergence = answer$convergence,
        diagnostic = answer$diagnostic
    )
}
