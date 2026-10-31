## Reproducible combined figure with a 265 Ma geological time axis.
## Run from the repository root after installing PastMoments and palaeoverse:
## Rscript PastMoments/inst/examples/geological-scale-example.R

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) > 0L) {
    args[[1L]]
} else {
    "PastMoments/geological-scale-example.png"
}

tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
tree$edge.length <- tree$edge.length * 132.5
trait <- c(a = 1, b = 2, c = 3, d = 4)

reconstruction <- PastMoments::reconstruct_empirical_moments(
    tree,
    trait,
    times = seq(5, 265, length.out = 100),
    sigma2 = 1,
    variance_probs = c(0.1, 0.5, 0.9)
)

grDevices::png(output_file, width = 1200, height = 1000, res = 120)
plot(
    reconstruction,
    time_orientation = "before_present",
    geological_scale = TRUE,
    geological_rank = c("stage", "period", "era"),
    geological_axis_args = list(abbr = TRUE, lab_size = 0.8)
)
grDevices::dev.off()
