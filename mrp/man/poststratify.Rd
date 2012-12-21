\name{poststratify}
\alias{poststratify}
\alias{poststratify,mrp-method}
\alias{poststratify,NWayData-method}
\alias{poststratify,jagsNWayData-method}

\title{Poststratification method}
\description{
  Poststratify multilevel regression model by an arbitrary number 
  of strata or \dQuote{ways}. By default this method returns a 
  single poststratified predicted value.
}
\usage{
  \S4method{poststratify}{mrp}(object, formula)
  \S4method{poststratify}{NWayData}(object, formula, fun, population)
  \S4method{poststratify}{jagsNWayData}(object, formula = NULL, fun = mean, population =
                 Census.NWay)
}
\arguments{
  \item{object}{An \code{mrp}, \code{NWayData}, or \code{jagsNWayData}
    object.}
  \item{formula}{
    A formula representation of the desired poststratification. 
    The formula is \code{NULL} on the left-hand side and 
    right-hand side variable names corresponding to the 
    \dQuote{ways} in the population data by which to poststratify. 
    The right-hand side can also be a character vector of 
    such names or a logical vector of length \dQuote{ways}.
    
    See example in \code{\link{mrp}}.}
  \item{fun}{The function (default=\emph{mean}) to summarize the
  collapsed dimensions.}
  \item{population}{An array or \code{NWayData} with dimensions matching
    \code{object}, used to produce population-weighted estimates from
    \code{jagsNWayData.}}
  }
  
\seealso{\code{\link{mrp-class}} for an example.
	\code{\link{mrp-class}} for other methods on the objects produced by \code{mrp()};
    \code{\link{plotmrp}} for how to plot poststratified results onto maps.}