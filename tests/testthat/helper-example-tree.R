example_tree <- function() {
    ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
}

example_trait <- function() {
    c(a = 1, b = 2, c = 3, d = 4)
}
