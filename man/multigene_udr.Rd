% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multigene_udr.R
\name{multigene_udr}
\alias{multigene_udr}
\title{Perform udr Analysis on Multigene Data}
\usage{
multigene_udr(
  combined_susie_list,
  coverage,
  independent_variant_list,
  n_random,
  n_null,
  seed,
  exclude_condition = NULL
)
}
\arguments{
\item{combined_susie_list}{A list containing the combined SuSiE and summary statistics results.}

\item{coverage}{A numeric vector representing the coverage values.}

\item{independent_variant_list}{A list of independent variants to be used as a filter.}

\item{n_random}{An integer specifying the number of random samples to generate.}

\item{n_null}{An integer specifying the number of null samples to generate.}

\item{seed}{An integer specifying the seed for random number generation.}

\item{exclude_condition}{A character vector specifying conditions to be excluded from the analysis. Defaults to NULL.}
}
\value{
A list containing the results of the prior, or NULL if conditions are not met.
}
\description{
Perform udr Analysis on Multigene Data
}
