\docType{class}
\name{mrp-class}
\alias{mrp-class}
\title{the mrp object}
\description{
  the mrp object
}
\section{Slots}{
  \describe{ \item{list("data")}{ the \code{data.frame}
  used in the multilevel regression step. This is created
  by summarizing the NWayData array in poll, and optionally
  joining additional columns of other predictors.
  \bold{Note:} the row order of the data object, which will
  be the row order of the fitted values, \bold{must be
  preserved} because it corresponds to the NWayData arrays
  used in poststratification. To add new columns onto the
  data.frame in the \code{data} slot, take care to preserve
  this ordering. Base \code{merge()} will almost invariably
  permute it in unpredictable ways: \bold{use
  \code{\link[plyr]{join}}} from plyr instead.  This is
  done automatically when data.frames are joined by the
  \code{add} argument to \code{\link{mrp}}.}\item{:}{ the
  \code{data.frame} used in the multilevel regression step.
  This is created by summarizing the NWayData array in
  poll, and optionally joining additional columns of other
  predictors. \bold{Note:} the row order of the data
  object, which will be the row order of the fitted values,
  \bold{must be preserved} because it corresponds to the
  NWayData arrays used in poststratification. To add new
  columns onto the data.frame in the \code{data} slot, take
  care to preserve this ordering. Base \code{merge()} will
  almost invariably permute it in unpredictable ways:
  \bold{use \code{\link[plyr]{join}}} from plyr instead.
  This is done automatically when data.frames are joined by
  the \code{add} argument to \code{\link{mrp}}.}
  \item{list("formula")}{ The formula used in the
  multilevel regression. The left-hand side is always
  \sQuote{response}.}\item{:}{ The formula used in the
  multilevel regression. The left-hand side is always
  \sQuote{response}.} \item{list("multilevelModel")}{ The
  multilevel regression model (class \code{mer} created by
  \code{\link[lme4]{glmer}})}\item{:}{ The multilevel
  regression model (class \code{mer} created by
  \code{\link[lme4]{glmer}})}
  \item{list("outcome")}{character. The name of the outcome
  variable.}\item{:}{character. The name of the outcome
  variable.} \item{list("poll")}{ An N-way array
  (\code{\link{NWayData-class}}) constructed from a
  survey.}\item{:}{ An N-way array
  (\code{\link{NWayData-class}}) constructed from a
  survey.} \item{list("population")}{ The population
  distribution, also of \code{\link{NWayData-class}}, with
  dimensions matching those of \code{poll}. This is
  intended to be a probability mass table (all entries sum
  to 1), but the package will also work with unnormalized
  population distributions.  The population is not used for
  the multilevel regression step.  It is used only in
  poststratification.}\item{:}{ The population
  distribution, also of \code{\link{NWayData-class}}, with
  dimensions matching those of \code{poll}. This is
  intended to be a probability mass table (all entries sum
  to 1), but the package will also work with unnormalized
  population distributions.  The population is not used for
  the multilevel regression step.  It is used only in
  poststratification.} }
}
\seealso{
  \code{\link{mrp}}, which produces \code{mrp} objects.
}
\keyword{classes}

