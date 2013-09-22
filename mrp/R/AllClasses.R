##' N-Dimensional Arrays for Multilevel Regression and Poststratification
##'
##' Arrays used in multilevel regression and poststratification are
##' low-dimensional (usually 2 or 3) slices of survey respondents or population
##' frequencies. For the multilevel regression step, dimensions of such arrays
##' are then summarized using methods here taking into account survey weights
##' and associated design effects.
##'
##'
##' @name NWayData-class
##' @aliases jagsNWayData-class
##' @aliases NWayData-class getDesignEffect,NWayData-method
##' getN,NWayData-method getNEffective,NWayData-method
##' getYbarWeighted,NWayData-method is.NWayData getDesignEffect getN
##' getNEffective getYbarWeighted
##' @docType class
##' @section Slots: \describe{ \item{.Data}{An array.} \item{type}{Character,
##' either \sQuote{poll} or \sQuote{pop} indicating the origin of the data
##' contained in the array.} \item{levels}{When the array is collapsed back to
##' a two-dimensional form, array dimension labels become character and any
##' ordering of the original factor levels is lost. To restore level names and
##' order, those attributes (in a \code{list}) are stored here.} }
##' @seealso \code{\link{mrp-class}}
##' @import grid
##' @import lattice
##' @keywords classes
setClass("NWayData",representation(type="character",levels="list"),contains="array")
setClass("jagsNWayData",representation(type="character",levels="list"),contains="array")


##' the mrp object
##'
##' @name mrp-class
##'
##' @aliases mrp-class
##' @docType class
##' @section Slots: \describe{ \item{list("data")}{ the \code{data.frame} used
##' in the multilevel regression step. This is created by summarizing the
##' NWayData array in poll, and optionally joining additional columns of other
##' predictors. \bold{Note:} the row order of the data object, which will be
##' the row order of the fitted values, \bold{must be preserved} because it
##' corresponds to the NWayData arrays used in poststratification. To add new
##' columns onto the data.frame in the \code{data} slot, take care to
##' preserve this ordering. Base \code{merge()} will almost invariably permute
##' it in unpredictable ways: \bold{use \code{\link[plyr]{join}}} from plyr
##' instead.  This is done automatically when data.frames are joined by the
##' \code{add} argument to \code{\link{mrp}}.}\item{:}{ the \code{data.frame}
##' used in the multilevel regression step. This is created by summarizing the
##' NWayData array in poll, and optionally joining additional columns of other
##' predictors. \bold{Note:} the row order of the data object, which will be
##' the row order of the fitted values, \bold{must be preserved} because it
##' corresponds to the NWayData arrays used in poststratification. To add new
##' columns onto the data.frame in the \code{data} slot, take care to
##' preserve this ordering. Base \code{merge()} will almost invariably permute
##' it in unpredictable ways: \bold{use \code{\link[plyr]{join}}} from plyr
##' instead.  This is done automatically when data.frames are joined by the
##' \code{add} argument to \code{\link{mrp}}.} \item{list("formula")}{ The
##' formula used in the multilevel regression. The left-hand side is always
##' \sQuote{response}.}\item{:}{ The formula used in the multilevel regression.
##' The left-hand side is always \sQuote{response}.}
##' \item{list("multilevelModel")}{ The multilevel regression model (class
##' \code{mer} created by \code{\link[lme4]{glmer}})}\item{:}{ The multilevel
##' regression model (class \code{mer} created by \code{\link[lme4]{glmer}})}
##' \item{list("outcome")}{character. The name of the outcome
##' variable.}\item{:}{character. The name of the outcome variable.}
##' \item{list("poll")}{ An N-way array (\code{\link{NWayData-class}})
##' constructed from a survey.}\item{:}{ An N-way array
##' (\code{\link{NWayData-class}}) constructed from a survey.}
##' \item{list("population")}{ The population distribution, also of
##' \code{\link{NWayData-class}}, with dimensions matching those of
##' \code{poll}. This is intended to be a probability mass table (all entries
##' sum to 1), but the package will also work with unnormalized population
##' distributions.  The population is not used for the multilevel regression
##' step.  It is used only in poststratification.}\item{:}{ The population
##' distribution, also of \code{\link{NWayData-class}}, with dimensions
##' matching those of \code{poll}. This is intended to be a probability mass
##' table (all entries sum to 1), but the package will also work with
##' unnormalized population distributions.  The population is not used for the
##' multilevel regression step.  It is used only in poststratification.} }
##' @seealso \code{\link{mrp}}, which produces \code{mrp} objects.
##' @keywords classes
setClass(Class="mrp",
         representation=representation(
         poll = "NWayData",
         data = "data.frame",
         formula = "formula",
         multilevelModel = "bmerMod",
         population = "NWayData",
         outcome="character"),

         validity=function (object) {
             if (is.null (object@data)) {
        stop ("[mrp: validation] flattened data is missing")
      }
      if (is.null (object@poll)) {
        stop("[mrp: validation] poll NWayData is missing.")
      }
      if (is.null (object@population)) {
        stop("[mrp: validation] pop NWayData is missing.")
      }
      if (is.null (object@formula)) {
        stop("[mrp: validation] formula is missing.")
      }
      return(TRUE)
    }
)
