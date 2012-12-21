\name{mrp-class}
\docType{class}
\alias{mrp-class}

\alias{getData,mrp-method}
\alias{getResponse,mrp-method}
\alias{getThetaHat,mrp-method}
\alias{getPopulation,mrp-method}
\alias{getPopulation}
\alias{setPopulation,mrp-method}
\alias{setPopulationOnes,mrp-method}
\alias{getFormula,mrp-method}
\alias{setFormula,mrp-method}
\alias{getFormula}
\alias{getResponse}
\alias{getThetaHat}
\alias{getEstimates}
\alias{getEstimates,mrp-method}
\alias{setFormula}
\alias{setPopulationOnes}
\alias{setPopulation}


\title{Multilevel regression and poststratification}
\description{
  The \code{mrp} class holds N-dimensional cross-classified arrays for a survey
  and a population, a 2-dimensional summary of survey data used to fit a
  multilevel model, and a fitted model. Methods described here operate
  on \code{mrp} objects which are typically created by calling the
  \code{\link{mrp}} function.
}

\section{Slots}{
  \describe{
    \item{\code{data}:}{ the
      \code{data.frame} used in the multilevel regression step. This is
      created by summarizing the NWayData array in poll, and optionally
      joining additional columns of other predictors. \bold{Note:} the row
      order of the data object, which will be the row order of the fitted
      values, \bold{must be preserved} because it corresponds to the NWayData arrays used
      in poststratification. To add new columns onto the data.frame in
  the \code{mrp@data} slot, take care to preserve this
  ordering. Base \code{merge()} will almost invariably permute it in
  unpredictable ways: \bold{use \code{\link[plyr]{join}}} from plyr instead.  
      This is done automatically when data.frames are
      joined by the \code{add} argument to \code{\link{mrp}}.}
    \item{\code{formula}:}{ The formula used in the multilevel
      regression. The left-hand side is always \sQuote{response}.}
    \item{\code{multilevelModel}:}{ The multilevel regression model 
      (class \code{mer} created by \code{\link[lme4]{glmer}})}
    \item{\code{poll}:}{ An N-way array (\code{\link{NWayData-class}})
      constructed from a survey.}
    \item{\code{population}:}{ The population distribution, also of
      \code{\link{NWayData-class}}, with dimensions matching those of
      \code{poll}. This is intended to be a probability 
      mass table (all entries sum to 1), but the package will 
      also work with unnormalized population distributions. 
      The population is not used for the multilevel regression step. 
      It is used only in poststratification.}
  }
}

\section{Methods}{
  \describe{
    \item{mr}{\code{signature(object = "mrp")}: update the formula and
      re-run the multilevel regression in the object.}
    \item{poststratify}{\code{signature(object = "mrp")}: perform
      poststratification using the fitted values of the multilevel model
      slot and the population frequencies contained in the
      population slot. When used with a
      \code{poststratification.specification} (formula, character, or
      logical; see above), gives values in each of the given cells of the
      N-way array.}
    \item{getThetaHat}{\code{signature(object = "mrp")}: estimates from the multilevel model 
	  in an N-way array corresponding to the dimension of the
      NWayData contained in slot \code{poll}.}
    \item{getEstimates}{alias for \code{getThetaHat}}
    \item{setPopulation}{\code{signature(object = "mrp")}: set
      population frequency data used for poststratification. This method
      checks that the dimension and names of the array correspond to those
      of slot \code{poll}.}
    \item{setPopulationOnes}{\code{signature(object = "mrp")}: To display model
      results without adjusting to a population frequency table, all
      dimensions can be set to 1.}
    \item{getNumberWays}{\code{signature(object = "mrp")}: returns the
      dimension of both \code{poll} and \code{population}. When a
      \dQuote{real} population array is not set (all=1), value for pop
      shows 0 here for.}
    \item{getFormula}{\code{signature(object = "mrp")}: returns the
      formula used to fit the multilevel regression.}
    \item{setFormula}{\code{signature(object = "mrp")}: set the formula; but
      you really should use the \code{mr} method instead.}
    \item{getPopulation}{\code{signature(object = "mrp")}: returns the
      \code{NWayData} population array.}
    \item{getData}{\code{signature(object = "mrp")}: return the
      \code{data.frame} used in the multilevel regression step.}
    \item{getResponse}{\code{signature(object = "mrp")}: returns a
      matrix containing the two-column \dQuote{yes/no} response used to
      fit the multilevel regression.}
  }
}
\seealso{\code{\link{mrp}}, which produces \code{mrp} objects.}

% \examples{
% }
\keyword{classes}
