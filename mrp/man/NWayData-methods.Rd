\docType{methods}
\name{getNumberWays}
\alias{getNumberWays}
\alias{getNumberWays,NWayData-method}
\alias{getNumberWays,mrp-method}
\alias{getNumberWays-method}
\title{Get Number of Ways for MRP analysis}
\description{
  A cross-classified dataset for multilevel regression and
  poststratification is an \eqn{N}-dimensional array. Each
  dimension in the array is one of the \dQuote{ways} by
  which it can later be poststratified. For various reasons
  it may be useful to query a \code{\link{mrp-class}}
  object or an \code{\link{NWayData-class}} object for
  these dimensions and their names.
}
\section{Methods}{
  \describe{ \item{list("signature(object = \"mrp\")")}{
  Returns a vector of length 2 with names \sQuote{poll} and
  \sQuote{pop}. In the special case where no population
  data has been assigned (and an array of 1s is used to
  return aggregated fitted values from the multilevel
  model) the value returned for \sQuote{pop} is 0.  }

  \item{list("signature(object = \"NWayData\")")}{ Returns
  the number of ways and an attribute \sQuote{ways}
  containing the character vector of names of the
  \dQuote{ways} by which it is constructed.  } }
}
\keyword{methods}

