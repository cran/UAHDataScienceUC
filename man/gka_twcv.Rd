% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkmeans.R
\name{gka_twcv}
\alias{gka_twcv}
\title{Total Within Cluster Variation (TWCV) computation}
\usage{
gka_twcv(data, k, population, centers)
}
\arguments{
\item{data}{a set of observations, presented as a matrix-like object where
every row is a new observation. The matrix is of size \code{n} by \code{m}.}

\item{k}{the number of clusters.}

\item{population}{a matrix of size \code{p} by \code{n} with the cluster
assignments for each observation.}

\item{centers}{a 3D array of size \code{p} by \code{k} by \code{m} with the
cluster centers for each individual in the population.}
}
\value{
a vector of size \code{p} with the total within cluster variation of
each individual in the population.
}
\description{
Total Within Cluster Variation (TWCV) computation
}
\author{
Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@edu.uah.es}
}
