# Empirical moments of Brownian traits conditional on tip values
# -----------------------------------------------------------------------------
# This script reorganizes the Brownian reconstruction code around a single
# computational engine and three public functions:
#
#   reconstruct_empirical_moments()
#
# The empirical mean is conditionally Gaussian, so its variance and quantiles
# are inexpensive once the conditional Gaussian distribution of the active
# lineages has been computed.
#
# For the empirical variance, its expectation and variance are available from
# analytic quadratic-form formulas. Its expected square root is evaluated by a
# deterministic fractional-moment integral. Quantiles are optional because they
# require numerical inversion of a non-central Gaussian quadratic-form
# distribution.
#
# Time is measured from the original root of the supplied tree. The root itself
# (time 0) is excluded from the default grid because the current edge-based
# convention represents lineages immediately after strictly positive times.
#
# The unknown root value is integrated out by re-rooting at the last tip.
# No package is attached globally. Qapprox or CompQuadForm is required only when
# empirical-variance quantiles are requested with the corresponding method.


## =============================================================================
## Basic linear algebra
## =============================================================================

.chol_spd <- function(M, name = "matrix") {
    M <- (M + t(M)) / 2
    tryCatch(
        chol(M),
        error = function(e) {
            stop(
                name, " is not numerically positive definite: ", e$message,
                call. = FALSE
            )
        }
    )
}


.chol_solve <- function(U, b) {
    # Solve M x = b from the upper Cholesky factor U, with M = t(U) %*% U.
    backsolve(U, forwardsolve(t(U), b))
}


.compute_variance_chol <- function(y, mean, UYY) {
    # The re-rooting tip is fixed. The contrast dimension is therefore length(y).
    x <- as.numeric(y - mean)
    as.numeric(crossprod(x, .chol_solve(UYY, x))) / length(x)
}


.check_psd <- function(M, name = "matrix", tol = 1e-10) {
    if (length(M) == 0L) {
        return(invisible(NA_real_))
    }

    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    scale <- max(1, max(abs(ev)))

    if (min(ev) < -tol * scale) {
        warning(sprintf(
            "%s is not numerically positive semidefinite: minimum eigenvalue = %g",
            name, min(ev)
        ))
    }

    invisible(min(ev))
}


## =============================================================================
## Tree utilities
## =============================================================================

.validate_tree <- function(tree) {
    if (!inherits(tree, "phylo")) {
        stop("tree must inherit from class 'phylo'.")
    }
    if (is.null(tree$edge) || !is.matrix(tree$edge) || ncol(tree$edge) != 2L) {
        stop("tree$edge must be a two-column parent-child matrix.")
    }
    if (is.null(tree$edge.length)) {
        stop("tree$edge.length is required.")
    }
    if (length(tree$edge.length) != nrow(tree$edge)) {
        stop("tree$edge.length must have one value per edge.")
    }
    if (any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
        stop("tree$edge.length must contain finite nonnegative values.")
    }
    if (is.null(tree$tip.label) || length(tree$tip.label) < 2L) {
        stop("tree must contain at least two labelled tips.")
    }
    if (is.null(tree$Nnode) || length(tree$Nnode) != 1L || tree$Nnode < 1L) {
        stop("tree$Nnode must be a positive scalar.")
    }
    if (anyDuplicated(tree$tip.label)) {
        stop("tree$tip.label contains duplicated labels.")
    }

    # The re-rooting code preserves edge rows and assumes the standard ape
    # numbering: tips are 1,...,Ntip and internal nodes follow consecutively.
    parents <- tree$edge[, 1]
    children <- tree$edge[, 2]
    tips <- setdiff(children, parents)
    internals <- setdiff(parents, tips)
    expected_tips <- seq_along(tree$tip.label)
    expected_internals <- length(tree$tip.label) + seq_len(tree$Nnode)

    if (!setequal(tips, expected_tips)) {
        stop("tree tips must be numbered 1,...,Ntip.")
    }
    if (!setequal(internals, expected_internals)) {
        stop("internal nodes must be numbered Ntip+1,...,Ntip+Nnode.")
    }
    if (max(tree$edge) != length(tree$tip.label) + tree$Nnode) {
        stop("tree node numbering is not consecutive.")
    }

    root <- setdiff(parents, children)
    if (length(root) != 1L) {
        stop("Could not identify a unique root.")
    }

    invisible(TRUE)
}


.compute_node_depths <- function(tree) {
    n_nodes <- length(tree$tip.label) + tree$Nnode
    if (n_nodes <= 0L) {
        stop("The tree contains no node.")
    }

    outgoing <- vector("list", n_nodes)
    for (e in seq_len(nrow(tree$edge))) {
        parent <- tree$edge[e, 1]
        outgoing[[parent]] <- c(outgoing[[parent]], e)
    }

    root <- setdiff(tree$edge[, 1], tree$edge[, 2])
    if (length(root) != 1L) {
        stop("Could not identify a unique root.")
    }
    root <- root[[1L]]

    depths <- numeric(n_nodes)
    stack <- integer(n_nodes)
    top <- 1L
    stack[top] <- root

    while (top > 0L) {
        node <- stack[top]
        top <- top - 1L

        for (e in outgoing[[node]]) {
            child <- tree$edge[e, 2]
            depths[child] <- depths[node] + tree$edge.length[e]
            top <- top + 1L
            stack[top] <- child
        }
    }

    depths
}


#' Compute phylogenetic node depths
#'
#' Compute the distance from the root to every tip and internal node of a
#' rooted phylogenetic tree.
#'
#' @param tree A rooted phylogenetic tree of class `"phylo"` with branch
#'   lengths and standard `"phylo"` node numbering.
#'
#' @return A numeric vector containing one depth per node. Tips occupy positions
#'   `1, ..., Ntip`, followed by the internal nodes in their usual `"phylo"`
#'   numbering. The root depth is zero.
#'
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
#' compute_node_depths(tree)
#'
#' @export
compute_node_depths <- function(tree) {
    .validate_tree(tree)
    .compute_node_depths(tree)
}


#' Make selected tree tips contemporaneous
#'
#' Adjust the terminal branches of selected tips so that they all end at the
#' same time, while leaving every internal and unselected branch unchanged.
#'
#' @param tree A rooted phylogenetic tree of class `"phylo"` with branch
#'   lengths.
#' @param tips Tip labels or one-based tip indices. A list is flattened before
#'   matching, which is convenient for reusing lists of focal tips.
#' @param target_time The common distance from the root at which the selected
#'   tips should end. When `NULL`, the greatest current depth among the selected
#'   tips is used, so branches are only extended.
#' @param tol Nonnegative numerical tolerance. A target up to `tol` before a
#'   selected tip's parent is treated as coinciding with that parent.
#'
#' @return The modified `"phylo"` object. Its `contemporaneous_tips` attribute
#'   is a data frame recording each selected tip's old and new times and the
#'   change in terminal branch length.
#'
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:0.5):1,(c:1,d:0.5):1);")
#' contemporaneous <- make_tips_contemporaneous(tree, c("b", "d"))
#' attr(contemporaneous, "contemporaneous_tips")
#'
#' @export
make_tips_contemporaneous <- function(tree, tips, target_time = NULL,
                                      tol = 1e-10) {
    .validate_tree(tree)

    if (is.list(tips)) {
        tips <- unlist(tips, recursive = TRUE, use.names = FALSE)
    }

    if (is.character(tips)) {
        if (length(tips) == 0L || anyNA(tips) || any(!nzchar(tips))) {
            stop("tips must contain at least one non-empty tip label.")
        }
        missing_tips <- setdiff(unique(tips), tree$tip.label)
        if (length(missing_tips) > 0L) {
            stop(
                "Unknown tip labels: ",
                paste(head(missing_tips, 10L), collapse = ", "),
                if (length(missing_tips) > 10L) " ..." else ""
            )
        }
        tip_indices <- match(unique(tips), tree$tip.label)
    } else if (is.numeric(tips)) {
        if (length(tips) == 0L || anyNA(tips) || any(!is.finite(tips)) ||
            any(tips != as.integer(tips))) {
            stop("Numeric tips must be finite integer tip indices.")
        }
        tip_indices <- unique(as.integer(tips))
        if (any(tip_indices < 1L | tip_indices > length(tree$tip.label))) {
            stop("Numeric tips must lie between 1 and the number of tips.")
        }
    } else {
        stop("tips must be a character vector, an integer vector, or a list.")
    }

    tol <- as.numeric(tol)
    if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
        stop("tol must be a finite nonnegative scalar.")
    }

    depths <- .compute_node_depths(tree)
    old_depths <- depths[tip_indices]

    if (is.null(target_time)) {
        target_time <- max(old_depths)
    } else {
        target_time <- as.numeric(target_time)
        if (length(target_time) != 1L || !is.finite(target_time) ||
            target_time < 0) {
            stop("target_time must be a finite nonnegative scalar.")
        }
    }

    incoming_edges <- match(tip_indices, tree$edge[, 2])
    if (anyNA(incoming_edges)) {
        stop("Could not find the terminal edge of every selected tip.")
    }

    parent_depths <- depths[tree$edge[incoming_edges, 1]]
    old_lengths <- tree$edge.length[incoming_edges]
    new_lengths <- target_time - parent_depths

    if (any(new_lengths < -tol)) {
        offending <- tree$tip.label[tip_indices[new_lengths < -tol]]
        stop(
            "target_time precedes the parent node of: ",
            paste(offending, collapse = ", ")
        )
    }

    new_lengths[new_lengths < 0] <- 0
    tree$edge.length[incoming_edges] <- new_lengths

    attr(tree, "contemporaneous_tips") <- data.frame(
        tip = tree$tip.label[tip_indices],
        old_time = old_depths,
        new_time = rep(target_time, length(tip_indices)),
        terminal_edge_change = new_lengths - old_lengths,
        stringsAsFactors = FALSE
    )

    tree
}


.get_mrca_function <- function(tree) {
    edge <- tree$edge
    n_nodes <- max(edge)

    parent <- rep(NA_integer_, n_nodes)
    children <- vector("list", n_nodes)

    for (e in seq_len(nrow(edge))) {
        p <- edge[e, 1]
        child <- edge[e, 2]
        parent[child] <- p
        children[[p]] <- c(children[[p]], child)
    }

    root <- setdiff(edge[, 1], edge[, 2])
    if (length(root) != 1L) {
        stop("Could not identify a unique root.")
    }
    root <- root[[1L]]

    level <- rep(NA_integer_, n_nodes)
    level[root] <- 0L
    stack <- root

    while (length(stack) > 0L) {
        u <- stack[[length(stack)]]
        stack <- stack[-length(stack)]

        for (v in children[[u]]) {
            level[v] <- level[u] + 1L
            stack <- c(stack, v)
        }
    }

    function(a, b) {
        while (level[a] > level[b]) {
            a <- parent[a]
        }
        while (level[b] > level[a]) {
            b <- parent[b]
        }
        while (a != b) {
            a <- parent[a]
            b <- parent[b]
        }
        a
    }
}


.compute_mrca_matrix <- function(tree) {
    get_mrca <- .get_mrca_function(tree)
    n_nodes <- max(tree$edge)
    mrca <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
    diag(mrca) <- seq_len(n_nodes)

    if (n_nodes > 1L) {
        for (i in seq_len(n_nodes - 1L)) {
            for (j in (i + 1L):n_nodes) {
                ancestor <- get_mrca(i, j)
                mrca[i, j] <- ancestor
                mrca[j, i] <- ancestor
            }
        }
    }

    mrca
}


.re_root_tree_at_last_tip <- function(tree) {
    # The original algorithm re-roots at the last tip and treats its observed
    # value as fixed. Keeping the last tip avoids a complete ape reindexing.
    node <- length(tree$tip.label)
    edge <- which(tree$edge[, 2] == node)

    if (length(edge) != 1L) {
        stop("The last tip must have exactly one incoming edge.")
    }

    prec <- node
    cur <- tree$edge[edge, 1]
    tree$edge[edge, 1] <- tree$edge[edge, 2]
    tree$edge[edge, 2] <- cur

    while (length(edge <- which(tree$edge[, 2] == cur)) > 0L) {
        if (length(edge) > 2L) {
            stop("Unexpected tree indexing while re-rooting.")
        }

        edge <- edge[tree$edge[edge, 1] != prec]
        if (length(edge) == 0L) {
            break
        }

        prec <- cur
        cur <- tree$edge[edge, 1]
        tree$edge[edge, 1] <- tree$edge[edge, 2]
        tree$edge[edge, 2] <- cur
    }

    tree$tip.label <- tree$tip.label[-node]
    tree$Nnode <- max(tree$edge) - length(tree$tip.label)
    tree
}


## =============================================================================
## Input normalization
## =============================================================================

#' Match a phylogenetic tree to trait data
#'
#' Extract one trait from a vector or two-dimensional table, remove tree tips
#' without a finite value, and return the trait vector in the exact order of
#' the pruned tree tips.
#'
#' @param tree A rooted phylogenetic tree of class `"phylo"` with branch
#'   lengths.
#' @param traits A trait vector, or a data frame, matrix, or two-dimensional
#'   array containing taxa and trait values. A named vector is matched to tree
#'   tips by name. An unnamed vector must have one value per tree tip and is
#'   assumed to follow `tree$tip.label` order.
#' @param trait_col Name or one-based index of the trait column. When
#'   `taxa_col` is supplied, this refers to the table after removal of the taxa
#'   column.
#' @param taxa_col `NULL` to obtain taxa from row names, or the name or
#'   one-based index of the column containing taxon labels.
#'
#' @return A list with `tree`, the pruned `"phylo"` object, and `trait`, a named
#'   numeric vector ordered as `tree$tip.label`.
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
#' traits <- data.frame(taxon = c("d", "b", "a", "c"),
#'                      size = c(4, NA, 1, 3))
#' prepared <- prepare_tree_trait(
#'     tree, traits, trait_col = "size", taxa_col = "taxon"
#' )
#' prepared$tree$tip.label
#' prepared$trait
#'
#' @export
prepare_tree_trait <- function(tree, traits, trait_col = 1L,
                               taxa_col = NULL) {
    .validate_tree(tree)

    is_trait_vector <- is.atomic(traits) && is.null(dim(traits))
    if (!is_trait_vector && !is.data.frame(traits) && !is.matrix(traits) &&
        !(is.array(traits) && length(dim(traits)) == 2L)) {
        stop(paste(
            "traits must be a vector, data.frame, matrix, or",
            "two-dimensional array."
        ))
    }

    if (is_trait_vector) {
        if (!is.null(taxa_col)) {
            stop("taxa_col cannot be used when traits is a vector.")
        }
        taxa <- names(traits)
        if (is.null(taxa)) {
            if (length(traits) != length(tree$tip.label)) {
                stop(
                    "An unnamed traits vector has length ", length(traits),
                    ", whereas the tree has ", length(tree$tip.label),
                    " tips."
                )
            }
            taxa <- tree$tip.label
        }
        traits_table <- data.frame(
            trait = unname(traits),
            stringsAsFactors = FALSE
        )
    } else {
        traits_table <- as.data.frame(traits, stringsAsFactors = FALSE)
    }

    resolve_column <- function(column, table, name) {
        if (is.character(column)) {
            if (length(column) != 1L || is.na(column) || !nzchar(column)) {
                stop(name, " must be one non-empty column name or index.")
            }
            matches <- which(names(table) == column)
            if (length(matches) != 1L) {
                stop(name, " must identify exactly one column in traits.")
            }
            return(matches)
        }
        if (is.numeric(column) && length(column) == 1L &&
            is.finite(column) && column == as.integer(column)) {
            column <- as.integer(column)
            if (column >= 1L && column <= ncol(table)) {
                return(column)
            }
        }
        stop(name, " must be a valid column name or one-based index.")
    }

    if (is_trait_vector) {
        # Taxa were obtained from vector names (or tree order) above.
    } else if (is.null(taxa_col)) {
        taxa <- rownames(traits_table)
        if (is.null(taxa)) {
            stop("taxa_col is NULL but traits has no row names.")
        }
    } else {
        taxa_index <- resolve_column(taxa_col, traits_table, "taxa_col")
        taxa <- traits_table[[taxa_index]]
        traits_table <- traits_table[, -taxa_index, drop = FALSE]
    }

    taxa <- as.character(taxa)
    if (anyNA(taxa) || any(!nzchar(taxa))) {
        stop("Taxon labels must be non-missing and non-empty.")
    }
    if (anyDuplicated(taxa)) {
        stop("Taxon labels must be unique.")
    }

    trait_index <- resolve_column(trait_col, traits_table, "trait_col")
    raw_values <- traits_table[[trait_index]]
    if (is.factor(raw_values)) {
        raw_values <- as.character(raw_values)
    }
    values <- suppressWarnings(as.numeric(raw_values))
    introduced_na <- is.na(values) & !is.na(raw_values)
    if (any(introduced_na)) {
        stop("trait_col contains values that cannot be converted to numeric.")
    }
    names(values) <- taxa

    keep_taxa <- intersect(names(values)[is.finite(values)], tree$tip.label)
    if (length(keep_taxa) < 2L) {
        stop("Fewer than two tree taxa have a finite trait value.")
    }

    pruned_tree <- ape::drop.tip(
        tree,
        setdiff(tree$tip.label, keep_taxa)
    )
    trait <- values[pruned_tree$tip.label]
    names(trait) <- pruned_tree$tip.label

    list(tree = pruned_tree, trait = trait)
}


.align_trait_to_tree <- function(tree, trait, name = "trait") {
    # Return a numeric vector ordered as tree$tip.label.
    if (is.data.frame(trait) || is.matrix(trait)) {
        if (ncol(trait) != 1L) {
            stop(name, " must be a vector or a one-column matrix/data.frame.")
        }
        rn <- rownames(trait)
        trait <- trait[, 1]
        if (is.null(names(trait)) && !is.null(rn)) {
            names(trait) <- rn
        }
    }

    trait_names <- names(trait)
    trait <- as.numeric(trait)
    if (!is.null(trait_names)) {
        names(trait) <- trait_names
    }

    if (is.null(names(trait))) {
        if (length(trait) != length(tree$tip.label)) {
            stop(
                name, " is unnamed and has length ", length(trait),
                ", whereas the tree has ", length(tree$tip.label), " tips."
            )
        }
        names(trait) <- tree$tip.label
        return(trait)
    }

    missing_tips <- setdiff(tree$tip.label, names(trait))
    if (length(missing_tips) > 0L) {
        stop(
            name, " is missing values for tree tips: ",
            paste(head(missing_tips, 10L), collapse = ", "),
            if (length(missing_tips) > 10L) " ..." else ""
        )
    }

    trait <- trait[tree$tip.label]
    if (anyNA(trait)) {
        stop(name, " contains NA after matching to tree$tip.label.")
    }

    names(trait) <- tree$tip.label
    trait
}


.normalize_clades <- function(tree, clades) {
    # clades = NULL means that all active lineages are treated as one group.
    if (is.null(clades)) {
        return(list(all = tree$tip.label))
    }

    # A character vector is accepted as a convenient single-clade shorthand.
    if (is.character(clades)) {
        clades <- list(clade = clades)
    }

    if (!is.list(clades) || length(clades) == 0L) {
        stop("clades must be NULL, a character vector, or a non-empty list.")
    }

    nm <- names(clades)
    if (is.null(nm)) {
        nm <- paste0("clade_", seq_along(clades))
    } else {
        empty <- is.na(nm) | !nzchar(nm)
        nm[empty] <- paste0("clade_", which(empty))
    }
    if (anyDuplicated(nm)) {
        stop("clade names must be unique.")
    }

    clades <- lapply(clades, function(z) unique(as.character(z)))
    names(clades) <- nm
    clades
}


.validate_probs <- function(probs, name = "probs", allow_null = TRUE) {
    if (is.null(probs)) {
        if (allow_null) {
            return(numeric(0))
        }
        stop(name, " cannot be NULL.")
    }

    probs <- as.numeric(probs)
    if (length(probs) == 0L || any(!is.finite(probs)) ||
        any(probs <= 0) || any(probs >= 1)) {
        stop(name, " must contain finite probabilities strictly between 0 and 1.")
    }
    if (anyDuplicated(probs)) {
        stop(name, " contains duplicated probabilities.")
    }

    probs
}


.probability_column_names <- function(probs) {
    if (length(probs) == 0L) {
        return(character(0))
    }

    pct <- format(100 * probs, scientific = FALSE, trim = TRUE, digits = 12)
    pct <- gsub("\\.", "_", pct)
    make.unique(paste0("q", pct), sep = "_")
}


.prepare_times <- function(total_depth, times = NULL, n_times = 2000L) {
    if (!is.finite(total_depth) || total_depth <= 0) {
        stop("The maximum tree depth must be strictly positive.")
    }

    if (is.null(times)) {
        if (!is.numeric(n_times) || length(n_times) != 1L ||
            !is.finite(n_times) || n_times < 1 || n_times != as.integer(n_times)) {
            stop("n_times must be a strictly positive integer.")
        }
        n_times <- as.integer(n_times)
        return(seq(total_depth / n_times, total_depth, length.out = n_times))
    }

    times <- as.numeric(times)
    if (length(times) == 0L || any(!is.finite(times))) {
        stop("times must contain finite values.")
    }

    tol <- 1e-12 * max(1, total_depth)
    if (any(times <= 0)) {
        stop("times must be strictly positive; time 0 is not represented.")
    }
    if (any(times > total_depth + tol)) {
        stop("times cannot exceed the maximum tip depth.")
    }

    times[times > total_depth] <- total_depth
    sort(unique(times))
}


## =============================================================================
