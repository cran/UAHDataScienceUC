% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkmeans.R
\name{gka_crossover}
\alias{gka_crossover}
\title{Crossover method i.e. K-Means Operator}
\usage{
gka_crossover(data, centers)
}
\arguments{
\item{data}{a set of observations, presented as a matrix-like object where
every row is a new observation. The matrix is of size \code{n} by \code{m}.}

\item{centers}{a matrix of size \code{k} by \code{m} with the cluster centers
for a specific individual in the population.}
}
\value{
a vector of size \code{n} with the cluster assignments for each
observation i.e. a new chromosome.
}
\description{
K-Means Operator (KMO) which replaces the crossover operator in
the Genetic K-Means algorithm (GKA).
}
\author{
Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@edu.uah.es}
}
