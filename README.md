# PastMoments

`PastMoments` reconstructs empirical means and variances of continuous traits
through phylogenetic time under Brownian models, conditionally on observed tip
values. It includes clade-specific graphical summaries.

## Installation

From the package directory:

```r
install.packages(".", repos = NULL, type = "source")
```

or from the source archive:

```r
install.packages(
    "PastMoments_0.1.tar.gz",
    repos = NULL,
    type = "source"
)
```

## Minimal example

```r
library(PastMoments)

tree <- ape::read.tree(
    text = "((a:1,b:1):1,(c:1,d:1):1);"
)
trait <- c(a = 1, b = 2, c = 3, d = 4)

result <- reconstruct_empirical_moments(
    tree,
    trait,
    times = seq(0.1, 2, length.out = 50),
    sigma2 = 1
)

print(result)
plot(result)
```

`reconstruct_empirical_moments()` is the single reconstruction entry point.
Its `quantities` argument accepts `"mean"`, `"variance"`, and
`"disparity_level"`; all three are selected by default. For example, a fast
mean-only reconstruction uses `quantities = "mean"`, while
`quantities = c("mean", "variance")` skips the disparity comparison.
Unknown entries are ignored with a warning as long as at least one supported
quantity remains.

## Expected empirical standard deviation

For an empirical variance `V(t)`, variance reconstruction rows include:

- `expected_variance`: `E[V(t) | data]`;
- `expected_sd`: the expected empirical standard deviation
  `E[sqrt(V(t)) | data]`.

`expected_sd` is evaluated without simulation. Internally,
`Q = (n - 1) V` is represented as a weighted sum of independent noncentral
chi-square variables plus a nonnegative deterministic term. With `L(t)` its
Laplace transform, the computation uses

```text
E[sqrt(Q)] = 1/sqrt(pi) integral_0^infinity t^(-1/2) {-L'(t)} dt,
L(t) = exp(-c t) product_j (1 + 2 lambda_j t)^(-1/2)
       exp{-lambda_j b_j^2 t / (1 + 2 lambda_j t)}.
```

This fractional moment is evaluated using deterministic one-dimensional
adaptive integration. Select this quantity as a mean-panel envelope with:

```r
plot(
    result,
    panels = "mean",
    mean_envelope = "expected_sd",
    mean_envelope_legend_label = "Expected empirical SD"
)
```

This envelope describes expected empirical spread among contemporaneous
lineages. It is not an uncertainty interval for the empirical mean.
Use `mean_envelope = "none"` to suppress it. The central curve is controlled
independently with `mean_display = "expectation"` or `"none"`. Likewise,
`variance_display` selects `"both"`, `"expectation"`, `"median"`, or `"none"`.

All reconstruction graphics use the standard `plot()` method. Reconstruction
objects retain their tree, aligned tip trait and normalized clades, so plots
need no separate biological inputs:

```r
plot(result, panels = "variance")
plot(
    result,
    panels = c("tree", "mean", "variance")
)
```

Common moment-panel aesthetics are supplied through `common_style`; values in
`mean_style` and `variance_style` override them for one panel. Tree-specific
settings belong in `tree_style`.

## Unconditional Brownian variance reference

The Brownian reference compares two independent variables at exactly
the same active lineages:

1. the reconstructed empirical variance conditional on the observed tips;
2. the unconditional empirical variance under Brownian motion with the same
   `sigma2`.

It is computed by default. It can also be selected explicitly before requesting
its separate plot panel:

```r
result_reference <- reconstruct_empirical_moments(
    tree,
    trait,
    times = seq(0.1, 2, length.out = 50),
    sigma2 = 1,
    quantities = c("mean", "variance", "disparity_level")
)

plot(
    result_reference,
    panels = c("tree", "mean", "variance", "disparity_level")
)
```

Variance rows then contain `disparity_level`, equal to
`P(reconstructed variance > Brownian variance)` for independent draws from the
conditional and unconditional distributions. The calculation uses Davies'
algorithm with Imhof's method as a deterministic fallback. It is a descriptive
probability of superiority, not a p-value, and its null distribution is not
generally uniform.

Selected non-contemporaneous tips can be aligned before reconstruction by
changing only their terminal branches:

```r
non_contemporaneous_tree <- ape::read.tree(
    text = "((a:1,b:0.5):1,(c:0.25,d:1):1);"
)
tree <- make_tips_contemporaneous(
    non_contemporaneous_tree,
    c("a", "b", "c", "d")
)
```

Variance quantiles are optional. They require `Qapprox` for
`quantile_method = "qapprox"`, or `CompQuadForm` for the Davies and Liu methods.
The Brownian reference requires `CompQuadForm`.

## Before a CRAN submission

Run:

```r
devtools::document()
devtools::test()
devtools::check()
```
