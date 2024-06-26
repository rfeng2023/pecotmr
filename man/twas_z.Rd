% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twas.R
\name{twas_z}
\alias{twas_z}
\title{Calculate TWAS z-score and p-value}
\usage{
twas_z(weights, z, R = NULL, X = NULL)
}
\arguments{
\item{weights}{A numeric vector of weights.}

\item{z}{A numeric vector of z-scores.}

\item{R}{An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.}

\item{X}{An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.}
}
\value{
A list containing the following elements:
\itemize{
  \item z: The TWAS z-score.
  \item pval: The corresponding p-value.
}
}
\description{
This function calculates the TWAS z-score and p-value given the weights, z-scores,
and optionally the correlation matrix (R) or the genotype matrix (X).
}
