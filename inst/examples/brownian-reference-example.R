## Exact Brownian empirical-variance reference and optional figure panel.
## Run from the repository root:
## Rscript PastMoments/inst/examples/brownian-reference-example.R

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) > 0L) {
    args[[1L]]
} else {
    "PastMoments/brownian-reference-example.png"
}

tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
trait <- c(a = 1, b = 2, c = 3, d = 4)
reconstruction <- PastMoments::reconstruct_empirical_moments(
    tree,
    trait,
    times = seq(0.1, 2, length.out = 100),
    sigma2 = 1,
    variance_probs = NULL,
    quantities = c("mean", "variance", "disparity_level")
)

grDevices::png(output_file, width = 1200, height = 1400, res = 140)
plot(
    reconstruction,
    variance_envelope_style = "none",
    panels = c("tree", "mean", "variance", "disparity_level")
)
grDevices::dev.off()
