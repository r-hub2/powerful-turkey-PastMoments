\name{make_clade_colors}
\alias{make_clade_colors}
\title{Assign Colours to Clades}
\description{Assign a distinct plotting colour to each reconstructed clade.}
\usage{make_clade_colors(clades, colors = NULL, palette = "Dark 3", alpha = 1)}
\arguments{
  \item{clades}{Character vector of unique clade names.}
  \item{colors}{Optional vector of colours, unnamed or named by clade.}
  \item{palette}{HCL palette name used when \code{colors} is \code{NULL}.}
  \item{alpha}{Transparency passed to \code{grDevices::hcl.colors}.}
}
\value{A named character vector of colours.}
