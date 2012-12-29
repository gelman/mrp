setGeneric ("mr", function (object,mr.formula,...) { standardGeneric ("mr")})
#setGeneric ("multilevelRegression", function (object) { standardGeneric ("multilevelRegression")})
setMethod (f="mr",
    signature=signature(object="mrp"),
    definition=function(object,mr.formula,...) {
      if(missing(mr.formula)) {
        fm <- object@formula
      } else {
        fm <- update.formula(object@formula, mr.formula)
        object@formula <- fm
      }
      response <- as.matrix(getResponse(object))
      object@multilevelModel <- bglmer(fm,
          data=object@data,
          family=binomial(link="logit"),...)
      return (object)
    } )
