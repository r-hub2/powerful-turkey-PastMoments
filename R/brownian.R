## REML setup and conditional Gaussian distribution
## =============================================================================

.prepare_reml_setup <- function(tree, trait, sigma2 = NULL) {
    .validate_tree(tree)
    trait <- .align_trait_to_tree(tree, trait)

    n_tips <- length(tree$tip.label)
    reroot <- n_tips
    retree <- .re_root_tree_at_last_tip(tree)

    depths_tree <- .compute_node_depths(tree)
    depths_retree <- .compute_node_depths(retree)
    n_tips_retree <- n_tips - 1L
    mrca_retree <- .compute_mrca_matrix(retree)

    # Brownian covariance matrix of the non-fixed tips on the re-rooted tree.
    tip_mrca <- mrca_retree[
        seq_len(n_tips_retree), seq_len(n_tips_retree), drop = FALSE
    ]
    MYY <- matrix(
        depths_retree[tip_mrca],
        nrow = n_tips_retree,
        ncol = n_tips_retree
    )

    UYY <- .chol_spd(MYY, "MYY")
    y_centered <- as.numeric(trait[-reroot] - trait[reroot])
    alpha_y <- .chol_solve(UYY, y_centered)

    if (is.null(sigma2)) {
        sigma2 <- .compute_variance_chol(
            trait[-reroot], trait[reroot], UYY
        )
    } else {
        sigma2 <- as.numeric(sigma2)
    }

    if (length(sigma2) != 1L || !is.finite(sigma2) || sigma2 < 0) {
        stop("sigma2 must be a finite nonnegative scalar.")
    }

    list(
        trait = trait,
        sigma2 = sigma2,
        reroot = reroot,
        retree = retree,
        depths_tree = depths_tree,
        depths_retree = depths_retree,
        n_tips_retree = n_tips_retree,
        total_depth = max(depths_tree),
        mrca_retree = mrca_retree,
        MYY = MYY,
        UYY = UYY,
        alpha_y = alpha_y
    )
}


.prepare_clade_membership <- function(tree, clades, depths_tree) {
    clades <- .normalize_clades(tree, clades)

    n_edges <- nrow(tree$edge)
    n_clades <- length(clades)
    tree_size <- max(tree$edge)

    # edge_parent[node] is the index of the edge ending at node.
    edge_parent <- rep(NA_integer_, tree_size)
    for (e in seq_len(n_edges)) {
        edge_parent[tree$edge[e, 2]] <- e
    }

    flags <- matrix(FALSE, nrow = n_edges, ncol = n_clades)
    colnames(flags) <- names(clades)
    time_min <- rep(NA_real_, n_clades)
    tip_indices <- vector("list", n_clades)
    names(tip_indices) <- names(clades)

    for (j in seq_along(clades)) {
        nm <- names(clades)[j]
        idx <- match(clades[[j]], tree$tip.label)
        missing_tips <- clades[[j]][is.na(idx)]

        if (length(missing_tips) > 0L) {
            warning(
                "Tips absent from tree for clade ", nm, ": ",
                paste(missing_tips, collapse = ", ")
            )
        }

        idx <- idx[!is.na(idx)]
        tip_indices[[j]] <- idx
        if (length(idx) == 0L) {
            next
        }

        clade_eff <- numeric(n_edges)

        # Mark every edge lying on a path from a selected terminal taxon to the
        # original root. At a given time, subsetting the active edges by this
        # flag gives the lineages contributing to the focal clade.
        for (tip in idx) {
            cur <- tip
            while (!is.na(edge_parent[cur]) && !flags[edge_parent[cur], j]) {
                flags[edge_parent[cur], j] <- TRUE
                cur <- tree$edge[edge_parent[cur], 1]
            }
        }

        # Count selected descendant tips along each ancestral path in order to
        # identify the crown time of the selected set.
        for (tip in idx) {
            cur <- tip
            while (!is.na(edge_parent[cur])) {
                clade_eff[edge_parent[cur]] <- clade_eff[edge_parent[cur]] + 1
                cur <- tree$edge[edge_parent[cur], 1]
            }
        }

        cur <- idx[[1L]]
        while (!is.na(edge_parent[cur]) &&
               clade_eff[edge_parent[cur]] < length(idx)) {
            cur <- tree$edge[edge_parent[cur], 1]
        }
        time_min[j] <- depths_tree[cur]
    }

    list(
        clades = clades,
        flags = flags,
        time_min = time_min,
        tip_indices = tip_indices
    )
}


.compute_active_conditional <- function(tree, setup, active_edges, time,
                                        check_psd = FALSE, psd_tol = 1e-10) {
    n <- length(active_edges)
    if (n == 0L) {
        return(list(
            mean = numeric(0),
            covariance0 = matrix(0, 0, 0),
            unconditional_covariance0 = matrix(0, 0, 0),
            active_edges = integer(0)
        ))
    }

    retree <- setup$retree
    depths_tree <- setup$depths_tree
    depths_retree <- setup$depths_retree
    mrca_retree <- setup$mrca_retree

    x_depth <- numeric(n)

    for (i in seq_along(active_edges)) {
        e <- active_edges[i]

        # Edge rows are preserved by the custom re-rooting operation. The test
        # below determines whether the orientation of edge e was reversed.
        if (retree$edge[e, 2] == tree$edge[e, 2]) {
            x_depth[i] <- depths_retree[retree$edge[e, 1]] +
                time - depths_tree[tree$edge[e, 1]]
        } else {
            x_depth[i] <- depths_retree[retree$edge[e, 1]] -
                time + depths_tree[tree$edge[e, 2]]
        }
    }

    active_children <- retree$edge[active_edges, 2]
    mrca_xx <- mrca_retree[
        active_children, active_children, drop = FALSE
    ]
    MXX <- matrix(depths_retree[mrca_xx], nrow = n, ncol = n)
    child_by_row <- matrix(active_children, nrow = n, ncol = n)
    depth_by_row <- matrix(x_depth, nrow = n, ncol = n)
    row_is_ancestor <- mrca_xx == child_by_row
    col_is_ancestor <- mrca_xx == t(child_by_row)
    MXX[row_is_ancestor] <- depth_by_row[row_is_ancestor]
    MXX[col_is_ancestor] <- t(depth_by_row)[col_is_ancestor]

    tip_indices <- seq_len(setup$n_tips_retree)
    mrca_xy <- mrca_retree[active_children, tip_indices, drop = FALSE]
    MXY <- matrix(
        depths_retree[mrca_xy],
        nrow = n,
        ncol = setup$n_tips_retree
    )
    child_is_ancestor <- mrca_xy == matrix(
        active_children,
        nrow = n,
        ncol = setup$n_tips_retree
    )
    x_depth_matrix <- matrix(
        x_depth,
        nrow = n,
        ncol = setup$n_tips_retree
    )
    MXY[child_is_ancestor] <- x_depth_matrix[child_is_ancestor]

    mean_vec <- setup$trait[setup$reroot] +
        as.vector(MXY %*% setup$alpha_y)

    # Schur complement computed from the Cholesky factor of MYY.
    Z <- forwardsolve(t(setup$UYY), t(MXY))
    C0 <- MXX - crossprod(Z)
    C0 <- (C0 + t(C0)) / 2

    if (isTRUE(check_psd)) {
        .check_psd(
            C0, sprintf("conditional covariance at time %g", time), psd_tol,
            covariance_scale = max(abs(MXX))
        )
    }

    list(
        mean = mean_vec,
        covariance0 = C0,
        # Brownian covariance at the same active points before conditioning on
        # the observed tips. A fixed root value is sufficient because only
        # centered contrasts enter empirical variance.
        unconditional_covariance0 = (MXX + t(MXX)) / 2,
        active_edges = active_edges
    )
}
