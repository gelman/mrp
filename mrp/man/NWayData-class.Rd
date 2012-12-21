\name{NWayData-class}
\docType{class}
\alias{NWayData-class}

\alias{getDesignEffect,NWayData-method}
\alias{getN,NWayData-method}
\alias{getNEffective,NWayData-method}
\alias{getYbarWeighted,NWayData-method}
\alias{is.NWayData}

\alias{getDesignEffect}
\alias{getN}
\alias{getNEffective}
\alias{getYbarWeighted}  

\title{N-Dimensional Arrays for Multilevel Regression and Poststratification}
\description{
  Arrays used in multilevel regression and poststratification are
  low-dimensional (usually 2 or 3) slices of survey respondents or
  population frequencies. For the multilevel regression step, dimensions
  of such arrays are then summarized using methods here taking
  into account survey weights and associated design effects.
}

\section{Slots}{
  \describe{
    \item{.Data}{An array.}
    \item{type}{Character, either \sQuote{poll} or \sQuote{pop}
      indicating the origin of the data contained in the array.}
    \item{levels}{When the array is collapsed back to a two-dimensional
      form, array dimension labels become character and any ordering of the
      original factor levels is lost. To restore level names and order,
      those attributes (in a \code{list}) are stored here.}
  }
}
\section{Methods}{
  \describe{
    \item{getNumberWays}{}
    \item{getYbarWeighted}{}
    \item{getDesignEffect}{}
    \item{getN}{}
    \item{getNEffective}{}
    \item{getDesignEffectByCell}{}
    }}
\seealso{\code{\link{mrp-class}}}
\author{Andrew Gelman <gelman@stat.columbia.edu>,
  Daniel Lee <bearlee@alum.mit.edu>,
  Yu-Sung Su <ys463@columbia.edu>,
  Michael Malecki <malecki@wustl.edu>
}
% \examples{
% }
\keyword{classes}
