.mr <- function(object,mr.formula.update,...) {
    dots <- list(...)
    if(is.null(dots$family)) {
        family <- binomial(link="logit")
    }
    if(missing(mr.formula.update)) {
        fm <- object@formula
    } else {
        fm <- update.formula(object@formula, mr.formula.update)
        object@formula <- fm
    }
    response <- as.matrix(getResponse(object))
    object@multilevelModel <- bglmer(fm,
                                     data=object@data, family=family, ...)
    return (object)
}
##'
##' ##' Run Multilevel Regression step of MRP Analysis
##'
##' Run a (binomial) multilevel regression in survey data for later
##' poststratification.
##'
##' @name mr
##' @docType methods
##' @aliases mr mr,mrp-method
##' @param object A \code{mrp} object.
##' @param mr.formula A formula specification for the multilevel model to run
##' in the prepared data. The left-hand side should always be
##' \sQuote{\code{response}}. For convenience, the formula is handled by
##' \code{update.formula} so that \code{.} indicates the current formula
##' contents on either side of the \code{~}, e.g., \code{.~.+newVar}. The
##' initial default formula is constructed as just an intercept term for each
##' of the variables in the main formula specification
##' (\code{(1|way1)+(1|way2)} etc.)
##' @param \dots Additional arguments to be passed to the multilevel regression
##' step, which uses \code{\link[blme]{bglmer}}.
##' @seealso \code{\link{mrp-class}} for an example.  \code{\link{mrp-class}}
##' for other methods on the objects produced by \code{mrp()};
##' \code{\link{plotmrp}} for how to plot poststratified results onto maps.
##' @export
setGeneric ("mr", function (object,mr.formula.update,...) { standardGeneric ("mr")})
setMethod (f="mr",
    signature=signature(object="mrp"),
    definition=.mr )
