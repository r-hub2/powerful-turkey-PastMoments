#' Coerce a reconstruction to a data frame
#'
#' @param x An empirical-moment reconstruction.
#' @param row.names Optional row names.
#' @param optional Logical; forwarded to [base::as.data.frame()].
#' @param ... Additional arguments forwarded to [base::as.data.frame()].
#'
#' @return The reconstruction table stored in `x$data`.
#' @export
as.data.frame.empirical_moment_reconstruction <- function(
    x, row.names = NULL, optional = FALSE, ...
) {
    .validate_reconstruction_object(x)
    as.data.frame(
        x$data, row.names = row.names, optional = optional, ...
    )
}


#' @export
print.empirical_moment_reconstruction <- function(x, ...) {
    .validate_reconstruction_object(x)
    quantities <- unique(as.character(x$data$quantity))

    cat("<empirical_moment_reconstruction>\n")
    cat("  model: ", x$model$name, "\n", sep = "")
    cat("  sigma2: ", format(x$model$sigma2), "\n", sep = "")
    cat("  tips: ", length(x$tree$tip.label), "\n", sep = "")
    cat("  clades: ", paste(names(x$clades), collapse = ", "), "\n", sep = "")
    cat("  quantities: ", paste(quantities, collapse = ", "), "\n", sep = "")
    if ("disparity_level" %in% names(x$data)) {
        cat("  disparity level: data$disparity_level\n")
    }
    cat("  rows: ", nrow(x$data), "\n", sep = "")
    invisible(x)
}


.validate_reconstruction_object <- function(x) {
    if (!inherits(x, "empirical_moment_reconstruction") || !is.list(x)) {
        stop("x must be an empirical_moment_reconstruction object.")
    }

    required <- c("data", "tree", "trait", "clades", "model", "settings")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0L) {
        stop(
            "The reconstruction object is missing: ",
            paste(missing, collapse = ", "), "."
        )
    }

    .validate_reconstruction_plot_data(x$data)
    .validate_tree(x$tree)
    x$trait <- .align_trait_to_tree(x$tree, x$trait, "x$trait")
    normalized_clades <- .normalize_clades(x$tree, x$clades)
    if (!identical(names(normalized_clades), names(x$clades))) {
        stop("x$clades has invalid or unstable clade names.")
    }

    result_clades <- unique(as.character(x$data$clade))
    missing_results <- setdiff(names(x$clades), result_clades)
    if (length(missing_results) > 0L) {
        stop(
            "x$data contains no reconstruction rows for clades: ",
            paste(missing_results, collapse = ", "), "."
        )
    }

    invisible(TRUE)
}
