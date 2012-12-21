
setGeneric ("poststratify", function (object, formula=NULL, ...) { standardGeneric ("poststratify")})
setMethod (f="poststratify",
    signature=signature(object="mrp"),
    definition=function (object, formula=NULL) {
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
      stopifnot (object@population != numeric(0))

      poststratified <- getThetaHat(object) * object@population

      if(!is.logical(spec)){
        groups <- match(spec,
            attr(getNumberWays(object@poll),"ways"))
      } else {
        groups <- which (spec == TRUE)
      }
      if (length(groups) == 0) {
        return (sum (poststratified, na.rm=TRUE) / sum (object@population))
      } else {
        ans <- (apply (poststratified, groups, sum, na.rm=TRUE) /
              apply(object@population, groups, sum))
        ans[is.nan(ans)] <- NA
        return(ans)
      }
    })


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
