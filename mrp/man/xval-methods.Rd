\docType{methods}
\name{xval}
\alias{xval}
\alias{xval,mrp-method}
\title{Run k-fold cross validations on mrp model}
\arguments{
  \item{object}{A \code{mrp} object.}

  \item{formula}{A formula specification for the multilevel
  model on which the k-folds cross validation is run. The
  default is the formula of the \code{mrp} object. The
  left-hand side should always be \sQuote{\code{response}}.
  For convenience, the formula is handled by
  \code{update.formula} so that \code{.} indicates the
  current formula contents on either side of the \code{~},
  e.g., \code{.~.+newVar}. The initial default formula is
  constructed as just an intercept term for each of the
  variables in the main formula specification
  (\code{(1|way1)+(1|way2)} etc.)}

  \item{folds}{Number of folds for cross valiation. Default
  is 4.}

  \item{loss.type}{Type of loss measure used in cross
  validation. Currently \dQuote{log} is supported.}

  \item{\dots}{Additional arguments to be passed to the
  multilevel regression step, which uses
  \code{\link[blme]{bglmer}} by default.}
}
\value{
  a list with named elements: MulRes, MulFormula, ClaRes,
  ClaFormula, and LB WEI-- DOC THIS RETURN PLEASE
}
\description{
  Run a (binomial) multilevel regression in survey data for
  later poststratification.
}

