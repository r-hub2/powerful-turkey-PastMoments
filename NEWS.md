# PastMoments 0.1

* Initial package structure.
* Brownian reconstruction of empirical means and variances.
* Added `make_tips_contemporaneous()` to align selected tips by changing only
  their terminal branches.
* Optional quadratic-form quantiles and clade-aware plotting functions.
* Added deterministic reconstruction of the expected empirical standard
  deviation, Jensen diagnostics, and selectable empirical-spread envelopes.
* Consolidated reconstruction graphics into the standard `plot()` method,
  with selectable tree, mean, and variance panels and shared or panel-specific
  style lists. Removed the redundant exported plotting functions.
* Reconstruction results are now self-contained list objects holding their
  data, tree, aligned trait, clades, model, and calculation settings. `plot()`
  reads these inputs directly from the result.
* Mean plotting now separates the central curve (`mean_display`), envelope
  quantity (`mean_envelope`), and envelope rendering
  (`mean_envelope_style`). Variance central curves are selected with
  `variance_display`; the former interval and median options were replaced by
  the corresponding envelope and display options.
* Renamed the plot argument `time_scale` to `time_orientation`.
* Added the exact unconditional Brownian empirical-variance reference,
  selected with `quantities = "disparity_level"`, optional time-dependent
  `disparity_level` output, and an optional `"disparity_level"` entry in
  `plot()`'s `panels` argument.
* Consolidated reconstruction calculations into
  `reconstruct_empirical_moments()`. Its `quantities` argument selects
  `"mean"`, `"variance"`, and/or `"disparity_level"`, with all three enabled
  by default. Removed the separate mean and variance reconstruction functions
  and the former logical `disparity_level` argument.
* Geological scales can now select and stack any GTS2020 ranks with
  `geological_rank` (for example `c("stage", "period", "era")`), and pass
  label and rendering options to `palaeoverse::axis_geo()` through
  `geological_axis_args`. The first requested rank supplies the transparent
  geological backgrounds.
