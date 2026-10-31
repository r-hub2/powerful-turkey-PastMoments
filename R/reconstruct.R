## =============================================================================
## Public interface
## =============================================================================

#' Reconstruct empirical moments through phylogenetic time
#'
#' Computes any combination of the conditional empirical mean, empirical
#' variance, and disparity level during one traversal of the time grid.
#'
#' @param tree A rooted phylogenetic tree of class `"phylo"` with branch lengths.
#' @param trait Numeric tip values, preferably named by `tree$tip.label`.
#' @param clades `NULL` for all lineages, a character vector for one focal set,
#'   or a named list of character vectors for several focal sets.
#' @param times Positive times measured from the original root. If `NULL`, a
#'   regular grid of `n_times` values is used.
#' @param n_times Number of default grid points when `times` is `NULL`.
#' @param sigma2 Brownian variance parameter. If `NULL`, it is estimated by REML.
#' @param quantities Character vector selecting `"mean"`, `"variance"`, and/or
#'   `"disparity_level"`. All three are computed by default. Selecting
#'   `"disparity_level"` also produces empirical-variance rows because the
#'   disparity level is defined from that distribution. Unsupported values are
#'   ignored with a warning; an error is raised only when none remain.
#' @param mean_probs Quantile probabilities for the empirical mean, or `NULL`.
#' @param variance_probs Quantile probabilities for the empirical variance, or
#'   `NULL` to skip their numerical calculation.
#' @param quantile_method Numerical method used for variance quantiles:
#'   `"qapprox"`, `"davies"`, or `"liu"`.
#' @param verbose Display a progress bar.
#' @param check_psd Check conditional covariance matrices numerically.
#' @param psd_tol Relative tolerance used by the PSD diagnostic.
#' @param eigen_tol Relative eigenvalue tolerance in quadratic decompositions.
#'
#' @return An `empirical_moment_reconstruction` object. Its `data` component has
#'   one row per time and clade for each requested moment. Variance rows contain
#'   `expected_variance` and `expected_sd`; when requested, they also contain
#'   `disparity_level`. Quantile columns are named from their probabilities, for
#'   example `q25`, `q50`, and `q75`.
#'
#' @references Didier, G. (draft). *Phylogenetic reconstruction of trait
#'   summary statistics and disparity dynamics*.
#'
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
#' trait <- c(a = -1, b = 0.5, c = 2, d = 3)
#' result <- reconstruct_empirical_moments(
#'     tree, trait, times = c(0.5, 1, 2), sigma2 = 1
#' )
#' result$data
#'
#' @export
reconstruct_empirical_moments <- function(
    tree,
    trait,
    clades = NULL,
    times = NULL,
    n_times = 2000L,
    sigma2 = NULL,
    quantities = c("mean", "variance", "disparity_level"),
    mean_probs = c(0.25, 0.50, 0.75),
    variance_probs = NULL,
    quantile_method = c("qapprox", "davies", "liu"),
    verbose = FALSE,
    check_psd = FALSE,
    psd_tol = 1e-10,
    eigen_tol = 1e-12
) {
    valid_quantities <- c("mean", "variance", "disparity_level")
    if (!is.character(quantities) || length(quantities) == 0L ||
        anyNA(quantities) || any(!nzchar(quantities))) {
        stop(
            "quantities must select at least one of 'mean', 'variance', ",
            "or 'disparity_level'.",
            call. = FALSE
        )
    }
    quantities <- unique(quantities)
    unsupported <- setdiff(quantities, valid_quantities)
    if (length(unsupported) > 0L) {
        warning(
            "Unsupported quantities were ignored and not computed: ",
            paste(sprintf("'%s'", unsupported), collapse = ", "), ".",
            call. = FALSE
        )
        quantities <- intersect(quantities, valid_quantities)
    }
    if (length(quantities) == 0L) {
        stop(
            "No supported quantity remains to compute. Select at least one ",
            "of 'mean', 'variance', or 'disparity_level'.",
            call. = FALSE
        )
    }
    compute_mean <- "mean" %in% quantities
    disparity_level <- "disparity_level" %in% quantities
    compute_variance <- "variance" %in% quantities || disparity_level

    engine <- .reconstruct_empirical(
        tree = tree,
        trait = trait,
        clades = clades,
        times = times,
        n_times = n_times,
        sigma2 = sigma2,
        compute_mean = compute_mean,
        compute_variance = compute_variance,
        disparity_level = disparity_level,
        mean_probs = mean_probs,
        variance_probs = variance_probs,
        quantile_method = quantile_method,
        verbose = verbose,
        check_psd = check_psd,
        psd_tol = psd_tol,
        eigen_tol = eigen_tol
    )

    out <- if (compute_mean && compute_variance) {
        .bind_reconstruction_tables(engine$mean, engine$variance)
    } else if (compute_mean) {
        engine$mean
    } else {
        engine$variance
    }
    .new_reconstruction(out, engine)
}
