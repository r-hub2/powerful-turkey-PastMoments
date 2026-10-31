## =============================================================================
## Shared reconstruction engine
## =============================================================================

.build_reconstruction_table <- function(quantity, times, clade_names,
                                        n_lineages, expectation, variance,
                                        quantiles, probs) {
    q_names <- .probability_column_names(probs)
    out <- vector("list", length(clade_names))

    for (j in seq_along(clade_names)) {
        d <- data.frame(
            time = times,
            clade = clade_names[[j]],
            n_lineages = n_lineages[, j],
            quantity = quantity,
            expectation = expectation[, j],
            variance = variance[, j],
            stringsAsFactors = FALSE
        )

        if (length(probs) > 0L) {
            for (h in seq_along(probs)) {
                d[[q_names[h]]] <- quantiles[, j, h]
            }
        }

        out[[j]] <- d
    }

    out <- do.call(rbind.data.frame, out)
    rownames(out) <- NULL
    out
}


.reconstruct_empirical <- function(tree, trait, clades = NULL,
                                   times = NULL, n_times = 2000L,
                                   sigma2 = NULL,
                                   compute_mean = TRUE,
                                   compute_variance = TRUE,
                                   disparity_level = FALSE,
                                   mean_probs = c(0.25, 0.50, 0.75),
                                   variance_probs = NULL,
                                   quantile_method = c("qapprox", "davies", "liu"),
                                   verbose = FALSE,
                                   check_psd = FALSE,
                                   psd_tol = 1e-10,
                                   eigen_tol = 1e-12) {
    if (!is.logical(compute_mean) || length(compute_mean) != 1L || is.na(compute_mean) ||
        !is.logical(compute_variance) || length(compute_variance) != 1L ||
        is.na(compute_variance) || (!compute_mean && !compute_variance)) {
        stop("At least one of compute_mean and compute_variance must be TRUE.")
    }
    if (!is.logical(disparity_level) ||
        length(disparity_level) != 1L || is.na(disparity_level)) {
        stop("disparity_level must be TRUE or FALSE.")
    }
    if (disparity_level && !compute_variance) {
        stop("disparity_level requires compute_variance = TRUE.")
    }
    if (disparity_level &&
        !requireNamespace("CompQuadForm", quietly = TRUE)) {
        stop(
            "Package 'CompQuadForm' is required when quantities includes ",
            "'disparity_level'.",
            call. = FALSE
        )
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE.")
    }
    if (!is.logical(check_psd) || length(check_psd) != 1L || is.na(check_psd)) {
        stop("check_psd must be TRUE or FALSE.")
    }
    if (!is.numeric(psd_tol) || length(psd_tol) != 1L ||
        !is.finite(psd_tol) || psd_tol <= 0) {
        stop("psd_tol must be a strictly positive finite scalar.")
    }
    if (!is.numeric(eigen_tol) || length(eigen_tol) != 1L ||
        !is.finite(eigen_tol) || eigen_tol <= 0) {
        stop("eigen_tol must be a strictly positive finite scalar.")
    }

    mean_probs <- .validate_probs(mean_probs, "mean_probs", allow_null = TRUE)
    variance_probs <- .validate_probs(
        variance_probs, "variance_probs", allow_null = TRUE
    )
    quantile_method <- match.arg(quantile_method)

    setup <- .prepare_reml_setup(tree, trait, sigma2 = sigma2)
    times <- .prepare_times(setup$total_depth, times, n_times)
    prep <- .prepare_clade_membership(tree, clades, setup$depths_tree)

    clade_names <- names(prep$clades)
    nT <- length(times)
    nC <- length(clade_names)

    n_lineages <- matrix(
        0L, nrow = nT, ncol = nC,
        dimnames = list(NULL, clade_names)
    )

    mean_expectation <- mean_variance <- NULL
    mean_quantiles <- NULL
    if (compute_mean) {
        mean_expectation <- matrix(NA_real_, nT, nC)
        mean_variance <- matrix(NA_real_, nT, nC)
        mean_quantiles <- array(
            NA_real_, dim = c(nT, nC, length(mean_probs))
        )
    }

    empvar_expectation <- empvar_variance <- NULL
    empvar_expected_sd <- NULL
    empvar_superiority <- NULL
    empvar_quantiles <- NULL
    if (compute_variance) {
        empvar_expectation <- matrix(NA_real_, nT, nC)
        empvar_variance <- matrix(NA_real_, nT, nC)
        empvar_expected_sd <- matrix(NA_real_, nT, nC)
        if (disparity_level) {
            empvar_superiority <- matrix(NA_real_, nT, nC)
        }
        empvar_quantiles <- array(
            NA_real_, dim = c(nT, nC, length(variance_probs))
        )
    }

    parent_depth <- setup$depths_tree[tree$edge[, 1]]
    child_depth <- setup$depths_tree[tree$edge[, 2]]

    pb <- NULL
    if (isTRUE(verbose)) {
        pb <- utils::txtProgressBar(min = 0, max = nT, style = 3)
        on.exit(close(pb), add = TRUE)
    }

    for (k in seq_along(times)) {
        time <- times[k]
        # Preserve the parent-open/child-closed convention at node times, while
        # allowing for roundoff when equivalent root-to-tip paths are summed.
        time_tol <- 64 * .Machine$double.eps * time
        active_edges <- which(
            parent_depth < time - time_tol & child_depth >= time - time_tol
        )

        active <- .compute_active_conditional(
            tree,
            setup,
            active_edges,
            time,
            check_psd = check_psd,
            psd_tol = psd_tol
        )
        unconditional_spectrum_cache <- if (disparity_level) {
            new.env(parent = emptyenv())
        } else {
            NULL
        }

        for (j in seq_along(clade_names)) {
            if (length(prep$tip_indices[[j]]) == 0L ||
                !is.finite(prep$time_min[j]) ||
                time < prep$time_min[j] - time_tol) {
                next
            }

            idx <- which(prep$flags[active_edges, j])
            m <- length(idx)
            if (m == 0L) {
                # This can happen after the last fossil occurrence of a focal set.
                next
            }

            n_lineages[k, j] <- m
            mu_sub <- active$mean[idx]
            C0_sub <- active$covariance0[idx, idx, drop = FALSE]

            if (compute_mean) {
                sm <- .summarize_empirical_mean(
                    mu_sub, C0_sub, setup$sigma2, mean_probs
                )
                mean_expectation[k, j] <- sm$expectation
                mean_variance[k, j] <- sm$variance
                if (length(mean_probs) > 0L) {
                    mean_quantiles[k, j, ] <- sm$quantiles
                }
            }

            if (compute_variance) {
                sv <- tryCatch(
                    .summarize_empirical_variance(
                        mu_sub,
                        C0_sub,
                        setup$sigma2,
                        variance_probs,
                        quantile_method,
                        tol = eigen_tol,
                        covariance_scale = max(abs(
                            active$unconditional_covariance0[idx, idx, drop = FALSE]
                        ))
                    ),
                    error = function(e) {
                        stop(
                            sprintf(
                                "Empirical-variance summary failed at time %g for clade '%s': %s",
                                time, clade_names[[j]], conditionMessage(e)
                            ),
                            call. = FALSE
                        )
                    }
                )
                empvar_expectation[k, j] <- sv$expectation
                empvar_variance[k, j] <- sv$variance
                empvar_expected_sd[k, j] <- sv$expected_sd
                if (disparity_level && m >= 2L) {
                    cache_key <- paste(idx, collapse = ",")
                    if (exists(
                        cache_key, envir = unconditional_spectrum_cache,
                        inherits = FALSE
                    )) {
                        unconditional_spectrum <- get(
                            cache_key, envir = unconditional_spectrum_cache,
                            inherits = FALSE
                        )
                    } else {
                        C0_unconditional <-
                            active$unconditional_covariance0[
                                idx, idx, drop = FALSE
                            ]
                        unconditional_decomposition <- .empvar_decomposition(
                            rep(0, m), C0_unconditional, tol = eigen_tol
                        )
                        unconditional_spectrum <-
                            .empvar_spectral_representation(
                                unconditional_decomposition, setup$sigma2
                            )
                        assign(
                            cache_key, unconditional_spectrum,
                            envir = unconditional_spectrum_cache
                        )
                    }
                    conditional_spectrum <- .empvar_spectral_representation(
                        sv$decomposition, setup$sigma2
                    )
                    comparison <- .indefinite_quadratic_probability(
                        conditional_weights = conditional_spectrum$weights,
                        conditional_noncentrality =
                            conditional_spectrum$noncentrality,
                        unconditional_weights =
                            unconditional_spectrum$weights,
                        conditional_shift = conditional_spectrum$shift,
                        unconditional_shift = unconditional_spectrum$shift,
                        tol = eigen_tol
                    )
                    empvar_superiority[k, j] <- comparison$probability
                }
                if (length(variance_probs) > 0L) {
                    empvar_quantiles[k, j, ] <- sv$quantiles
                }
            }
        }

        if (!is.null(pb)) {
            utils::setTxtProgressBar(pb, k)
        }
    }

    result <- list()

    if (compute_mean) {
        result$mean <- .build_reconstruction_table(
            quantity = "empirical_mean",
            times = times,
            clade_names = clade_names,
            n_lineages = n_lineages,
            expectation = mean_expectation,
            variance = mean_variance,
            quantiles = mean_quantiles,
            probs = mean_probs
        )
    }

    if (compute_variance) {
        result$variance <- .build_reconstruction_table(
            quantity = "empirical_variance",
            times = times,
            clade_names = clade_names,
            n_lineages = n_lineages,
            expectation = empvar_expectation,
            variance = empvar_variance,
            quantiles = empvar_quantiles,
            probs = variance_probs
        )
        result$variance$expected_variance <- result$variance$expectation
        result$variance$expected_sd <- as.vector(empvar_expected_sd)
        if (disparity_level) {
            result$variance$disparity_level <- as.vector(
                empvar_superiority
            )
        }
    }

    result$model <- "Brownian"
    result$sigma2 <- setup$sigma2
    result$times <- times
    result$clades <- prep$clades
    result$mean_probs <- mean_probs
    result$variance_probs <- variance_probs
    result$quantile_method <- quantile_method
    result$quantities <- c(
        if (compute_mean) "mean",
        if (compute_variance) "variance",
        if (disparity_level) "disparity_level"
    )
    result$tree <- tree
    result$trait <- setup$trait
    result
}



.new_reconstruction <- function(data, engine) {
    structure(
        list(
            data = data,
            tree = engine$tree,
            trait = engine$trait,
            clades = engine$clades,
            model = list(
                name = engine$model,
                sigma2 = engine$sigma2
            ),
            settings = list(
                times = engine$times,
                mean_probs = engine$mean_probs,
                variance_probs = engine$variance_probs,
                quantile_method = engine$quantile_method,
                quantities = engine$quantities
            )
        ),
        class = "empirical_moment_reconstruction"
    )
}



.bind_reconstruction_tables <- function(x, y) {
    # Bind two reconstruction tables even if they use different quantile sets.
    all_names <- union(names(x), names(y))

    add_missing <- function(d) {
        missing <- setdiff(all_names, names(d))
        for (nm in missing) {
            d[[nm]] <- NA_real_
        }
        d[all_names]
    }

    out <- rbind(add_missing(x), add_missing(y))
    rownames(out) <- NULL
    out
}
