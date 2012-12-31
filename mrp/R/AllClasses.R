##' @exportPattern "."
##' @import grid
##' @import lattice
NULL

setClass("NWayData",representation(type="character",levels="list"),contains="array")
setClass("jagsNWayData",representation(type="character",levels="list"),contains="array")

setClass(Class="mrp",
    representation=representation(
        poll = "NWayData",
        data = "data.frame",
        formula = "formula",
        multilevelModel = "mer",
        population = "NWayData"),
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
