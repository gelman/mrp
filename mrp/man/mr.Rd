\name{mr}
\alias{mr}
\alias{mr,mrp-method}

\title{Run Multilevel Regression step of MRP Analysis}
\description{
  Run a (binomial) multilevel regression in survey data for later
  poststratification.
}

\usage{
  \S4method{mr}{mrp}(object, mr.formula, \dots)
}
\arguments{
    \item{object}{A \code{mrp} object.}  
    \item{mr.formula}{A formula specification for the multilevel model
      to run in the prepared data. The left-hand side should always be
	  \sQuote{\code{response}}. For convenience, the formula is handled by
      \code{update.formula} so that \code{.} indicates the 
      current formula contents on either side of the \code{~}, e.g.,
      \code{.~.+newVar}. The initial default formula is constructed as
      just an intercept term for each of the variables in the main
      formula specification (\code{(1|way1)+(1|way2)} etc.)}
    \item{\dots}{Additional arguments to be passed to the multilevel
      regression step, which uses \code{\link[lme4]{glmer}} by default. }
  }
\seealso{
	\code{\link{mrp-class}} for an example.
	\code{\link{mrp-class}} for other methods on the objects produced by \code{mrp()};
    \code{\link{plotmrp}} for how to plot poststratified results  onto maps.}