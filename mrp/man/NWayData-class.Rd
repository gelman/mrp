\docType{class}
\name{NWayData-class}
\alias{NWayData-class}
\alias{getDesignEffect}
\alias{getDesignEffect,NWayData-method}
\alias{getN}
\alias{getN,NWayData-method}
\alias{getNEffective}
\alias{getNEffective,NWayData-method}
\alias{getYbarWeighted}
\alias{getYbarWeighted,NWayData-method}
\alias{is.NWayData}
\alias{jagsNWayData-class}
\title{N-Dimensional Arrays for Multilevel Regression and Poststratification}
\description{
  Arrays used in multilevel regression and
  poststratification are low-dimensional (usually 2 or 3)
  slices of survey respondents or population frequencies.
  For the multilevel regression step, dimensions of such
  arrays are then summarized using methods here taking into
  account survey weights and associated design effects.
}
\section{Slots}{
  \describe{ \item{.Data}{An array.} \item{type}{Character,
  either \sQuote{poll} or \sQuote{pop} indicating the
  origin of the data contained in the array.}
  \item{levels}{When the array is collapsed back to a
  two-dimensional form, array dimension labels become
  character and any ordering of the original factor levels
  is lost. To restore level names and order, those
  attributes (in a \code{list}) are stored here.} }
}
\seealso{
  \code{\link{mrp-class}}
}
\keyword{classes}

