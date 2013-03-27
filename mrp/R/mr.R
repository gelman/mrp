setGeneric ("mr", function (object,mr.formula.update,...) { standardGeneric ("mr")})
#setGeneric ("multilevelRegression", function (object) { standardGeneric ("multilevelRegression")})
setMethod (f="mr",
    signature=signature(object="mrp"),
    definition=function(object,mr.formula.update,...) {
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
    } )
