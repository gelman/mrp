
.poststratify <- function (object, formula=NULL, population=NULL) {
    spec <- formula
    if(is.null(object@population)) {
        warning("Object does not contain population data;\nestimates returned instead.")
        return(getThetaHat(object));
    }
    if(is.null(spec)){
        spec <- rep(FALSE,getNumberWays(object@poll))
    }
    if(is.formula(spec)){
        spec <- attr(terms(spec),"term.labels")
    }
    if (is.null(population)) {
      population <- object@population
    }
    stopifnot (population != numeric(0))

    poststratified <- getThetaHat(object) * population

    if(!is.logical(spec)){
        groups <- match(spec,
                        attr(getNumberWays(object@poll),"ways"))
    } else {
        groups <- which (spec == TRUE)
    }
    if (length(groups) == 0) {
        return (sum (poststratified, na.rm=TRUE) / sum (population))
    } else {
        ans <- (apply (poststratified, groups, sum, na.rm=TRUE) /
                apply(population, groups, sum))
        ans[is.nan(ans)] <- NA
        return(ans)
    }
}

##' Poststratification method
##'
##' Poststratify multilevel regression model by an arbitrary number of strata
##' or \dQuote{ways}. By default this method returns a single poststratified
##' predicted value.
##'
##' @docType methods
##' @rdname poststratify
##' @aliases poststratify poststratify-methods poststratify,mrp-method poststratify,NWayData-method
##' @aliases poststratify,jagsNWayData-method
##' poststratify,jagsNWayData-method
##' @param object An \code{mrp}, \code{NWayData}, or \code{jagsNWayData}
##' object.
##' @param formula A formula representation of the desired poststratification.
##' The formula is \code{NULL} on the left-hand side and right-hand side
##' variable names corresponding to the \dQuote{ways} in the population data by
##' which to poststratify.  The right-hand side can also be a character vector
##' of such names or a logical vector of length \dQuote{ways}.
##' @param population An optional replacement population array
##'
##' See example in \code{\link{mrp}}.
##' @param fun The function (default=\emph{mean}) to summarize the collapsed
##' dimensions.
##' @param population An array or \code{NWayData} with dimensions matching
##' \code{object}, used to produce population-weighted estimates from
##' \code{jagsNWayData.}
##' @seealso \code{\link{mrp-class}} for an example.  \code{\link{mrp-class}}
##' for other methods on the objects produced by \code{mrp()};
##' \code{\link{plotmrp}} for how to plot poststratified results onto maps.
##' @export
setGeneric ("poststratify", function (object, formula=NULL, population=NULL, ...) { standardGeneric ("poststratify")})
setMethod (f="poststratify",
           signature=signature(object="mrp"),
           definition=.poststratify)

setMethod (f="poststratify",
           signature=signature(object="NWayData"),
           definition=function (object, formula=NULL, fun=mean,
           population=Census.NWay#, sort=NULL
           ) {
               if(object@type=="jagsNWayData"){
                   spec <- formula

                   if(is.null(spec)){
                       spec <- rep(FALSE,length(dimnames(object)))
                   }
                   if(is.formula(spec)){
                       spec <- attr(terms(spec),"term.labels")
                   }
                   stopifnot (population != numeric(0))

                   poststratified <- object * population

                   if(!is.logical(spec)){
                       groups <- match(spec,
                                       names(dimnames(object)) )
                   } else {
                       groups <- which (spec == TRUE)
                   }
                   if (length(groups) == 0) {
                       ans <- sum (poststratified, na.rm=TRUE) / sum (population)
                   } else {
                       ans <- (apply (poststratified, groups, sum, na.rm=TRUE) /
                               apply(population, groups, sum, na.rm=TRUE))
                       ans[is.nan(ans)] <- NA
                       ans <- new("NWayData",
                                  ans, type="poststratified",
                                  levels=object@levels[groups])
                   }
                   ## if (!is.null(sort)) {
                   ##   sortme <- rep(TRUE,length(dimnames(ans)))
                   ##   sortme[which(dimnames(ans) == sort)] <- FALSE
                   ##   ans <- apply(ans,sortme,sort)
                   ## }
                   return(ans)
               } else {
                   spec <- formula
                   if(is.null(spec)){
                       spec <- rep(FALSE, length(dim(object)) )
                   }
                   if(is.formula(spec)){
                       spec <- attr(terms(spec),"term.labels")
                   }
                   if(!is.logical(spec)){
                       groups <- match(spec,
                                       names(dimnames(object)) )
                   } else {
                       groups <- which (spec == TRUE)
                   }
                   if (length(groups) == 0) {
                       ans <- do.call(fun, args=list(object, na.rm=TRUE))
                   } else {
                       ans <- (apply (object, groups, fun, na.rm=TRUE))
                       ans[is.nan(ans)] <- NA
                       ans <- new("NWayData",
                                  ans, type="poststratified",
                                  levels=object@levels[groups])

                   }
                   ## if (!is.null(sort)) {
                   ##   browser()
                   ##   sortme <- rep(TRUE,length(dimnames(ans)))
                   ##   sortme[which(dimnames(ans) == sort)] <- FALSE
                   ##   #ans <- apply(ans,sortme,sort)
                   ## }
                   return(ans)
               }

           })

setMethod (f="poststratify",
           signature=signature(object="jagsNWayData"),
           definition=function (object, formula=NULL, fun=mean, population=Census.NWay) {
               spec <- formula

               if(is.null(spec)){
                   spec <- rep(FALSE,length(dimnames(object)))
               }
               if(is.formula(spec)){
                   spec <- attr(terms(spec),"term.labels")
               }
               stopifnot (population != numeric(0))

               poststratified <- object * population

               if(!is.logical(spec)){
                   groups <- match(spec,
                                   names(dimnames(object)) )
               } else {
                   groups <- which (spec == TRUE)
               }
               if (length(groups) == 0) {
                   return (sum (poststratified, na.rm=TRUE) / sum (population))
               } else {
                   ans <- (apply (poststratified, groups, sum, na.rm=TRUE) /
                           apply(population, groups, sum, na.rm=TRUE))
                   ans[is.nan(ans)] <- NA
                   return(ans)
               }
           })
