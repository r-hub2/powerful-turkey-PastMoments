test_that("compute_node_depths exposes root distances for every node", {
    tree <- ape::read.tree(text = "((a:1,b:2):3,c:4);")

    depths <- compute_node_depths(tree)

    expect_type(depths, "double")
    expect_length(depths, length(tree$tip.label) + tree$Nnode)
    expect_equal(depths, c(4, 5, 4, 0, 3))
})


test_that("compute_node_depths is independent of edge order and validates input", {
    tree <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
    order <- rev(seq_len(nrow(tree$edge)))
    tree$edge <- tree$edge[order, , drop = FALSE]
    tree$edge.length <- tree$edge.length[order]

    expect_equal(compute_node_depths(tree), c(4, 5, 4, 0, 3))
    expect_error(compute_node_depths(list()), "inherit from class 'phylo'")

    tree_without_lengths <- tree
    tree_without_lengths$edge.length <- NULL
    expect_error(compute_node_depths(tree_without_lengths),
                 "edge.length is required")
})


test_that("make_tips_contemporaneous extends selected terminal branches", {
    tree <- ape::read.tree(text = "((a:1,b:0.5):1,(c:0.25,d:1):1);")
    original_lengths <- tree$edge.length

    result <- make_tips_contemporaneous(tree, c("b", "c", "d"))
    report <- attr(result, "contemporaneous_tips")

    expect_s3_class(result, "phylo")
    expect_equal(report$tip, c("b", "c", "d"))
    expect_equal(report$old_time, c(1.5, 1.25, 2))
    expect_equal(report$new_time, rep(2, 3))
    expect_equal(report$terminal_edge_change, c(0.5, 0.75, 0))

    selected_edges <- match(match(c("b", "c", "d"), tree$tip.label),
                            tree$edge[, 2])
    expect_equal(result$edge.length[-selected_edges],
                 original_lengths[-selected_edges])
})


test_that("make_tips_contemporaneous accepts indices and an explicit target", {
    tree <- ape::read.tree(text = "((a:1,b:0.5):1,(c:0.25,d:1):1);")

    result <- make_tips_contemporaneous(tree, c(1L, 2L), target_time = 1.75)
    report <- attr(result, "contemporaneous_tips")

    expect_equal(report$tip, c("a", "b"))
    expect_equal(report$new_time, c(1.75, 1.75))
    expect_equal(result$edge.length[match(1:2, result$edge[, 2])],
                 c(0.75, 0.75))
})


test_that("make_tips_contemporaneous validates tips and target time", {
    tree <- example_tree()

    expect_error(make_tips_contemporaneous(tree, "missing"),
                 "Unknown tip labels")
    expect_error(make_tips_contemporaneous(tree, 0L),
                 "between 1 and the number of tips")
    expect_error(make_tips_contemporaneous(tree, "a", target_time = 0.5),
                 "precedes the parent node")
})


test_that("MRCA matrix caches all pairwise ancestors", {
    tree <- example_tree()
    get_mrca <- PastMoments:::.get_mrca_function(tree)
    mrca <- PastMoments:::.compute_mrca_matrix(tree)
    n_nodes <- max(tree$edge)

    expect_equal(dim(mrca), c(n_nodes, n_nodes))
    expect_equal(diag(mrca), seq_len(n_nodes))
    expect_equal(mrca, t(mrca))

    for (i in seq_len(n_nodes)) {
        for (j in seq_len(n_nodes)) {
            expect_equal(mrca[i, j], get_mrca(i, j))
        }
    }
})


test_that("prepare_tree_trait extracts, prunes, and aligns trait values", {
    tree <- example_tree()
    traits <- data.frame(
        taxon = c("d", "outside", "b", "a", "c"),
        mass = c(4, 10, NA, 1, 3),
        length = 5:1
    )

    prepared <- prepare_tree_trait(
        tree,
        traits,
        trait_col = "mass",
        taxa_col = "taxon"
    )

    expect_s3_class(prepared$tree, "phylo")
    expect_setequal(prepared$tree$tip.label, c("a", "c", "d"))
    expect_equal(names(prepared$trait), prepared$tree$tip.label)
    expect_equal(unname(prepared$trait),
                 traits$mass[match(prepared$tree$tip.label, traits$taxon)])
})


test_that("prepare_tree_trait supports row names and validates columns", {
    tree <- example_tree()
    traits <- data.frame(mass = c(4, 2, 1, 3))
    rownames(traits) <- c("d", "b", "a", "c")

    prepared <- prepare_tree_trait(tree, traits)
    expect_equal(prepared$trait, c(a = 1, b = 2, c = 3, d = 4))

    expect_error(prepare_tree_trait(tree, traits, trait_col = "missing"),
                 "exactly one column")
    expect_error(
        prepare_tree_trait(
            tree,
            data.frame(taxon = c("a", "a"), mass = c(1, 2)),
            taxa_col = "taxon"
        ),
        "unique"
    )
    expect_error(
        prepare_tree_trait(
            tree,
            data.frame(taxon = c("a", "b"), mass = c("one", "two")),
            taxa_col = "taxon"
        ),
        "converted to numeric"
    )
})


test_that("prepare_tree_trait supports named and unnamed vectors", {
    tree <- example_tree()

    named <- c(d = 4, outside = 10, b = NA, a = 1, c = 3)
    prepared <- prepare_tree_trait(tree, named)
    expect_setequal(prepared$tree$tip.label, c("a", "c", "d"))
    expect_equal(prepared$trait, named[prepared$tree$tip.label])

    unnamed <- c(1, 2, 3, 4)
    prepared <- prepare_tree_trait(tree, unnamed)
    expect_equal(prepared$trait, setNames(unnamed, tree$tip.label))

    expect_error(
        prepare_tree_trait(tree, unnamed[-1]),
        "whereas the tree has 4 tips"
    )
    expect_error(
        prepare_tree_trait(tree, named, taxa_col = 1),
        "cannot be used"
    )
})
