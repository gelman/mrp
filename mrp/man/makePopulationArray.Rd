\name{makePopulationArray}
\alias{makePopulationArray}
\title{Make an array suitable for poststratification}
\usage{
  makePopulationArray(population, pop.weights,
    population.varnames)
}
\arguments{
  \item{population}{A \code{data.frame} containing
  population (e.g. census) containing numeric weight and
  some factor variables for cross-classification}

  \item{pop.weights}{character. The column of the
  \code{population} data that contains frequencies or
  proportions (of the entire population) for the cells
  defined by population.varnames}

  \item{population.varnames}{character vector of names of
  columns describing the cells wanted in the array}
}
\value{
  NWayData (array with some extra stuff)
}
\description{
  Make an array suitable for poststratification
}

