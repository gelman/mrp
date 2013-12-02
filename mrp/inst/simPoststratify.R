## Uncertainty using sim() on lme4-fitted models
## Mostly by Jeffrey R. Lax (@jrllrj) and Michael Malecki @malecki
## Use at your own risk. This needs review and revision before
## it can move into mrp itself.
###################################

poststratify.sim <- function(object, formula, n.sim=2) {

  n.sim <- max(2,n.sim)
  pop <- getPopulation(object)
    dimnametemp <- dimnames(pop)
    dimnametemp$simnum <- c(1:n.sim)
    theta.hat.sim <- array (fitted(sim(getModel(object), n.sim),getModel(object) ),
                        dim=c(dim(pop),n.sim), dimnames=dimnametemp)
   # return(theta.hat.sim)
  pop.sim <- pop
    pop.sim@levels$simnum <- c(1:n.sim)
    pop.sim@.Data <-   array (rep(pop.sim@.Data,2), dim = c(dim(pop),n.sim), dimnames = dimnametemp )


  c( dim(pop.sim@.Data), "simnum" = n.sim)

       poststratified.sim <-  theta.hat.sim * pop.sim


     spec <- formula
      if(is.null(object@population)) {
        warning("Object does not contain population data;\nestimates returned instead.")
        return(theta.hat.sim);
      }
      if(is.null(spec)){
        spec <- rep(FALSE,getNumberWays(object@poll))
      }
      if(is.formula(spec)){
        spec <- attr(terms(spec),"term.labels")
      }
      stopifnot (object@population != numeric(0))

      if(!is.logical(spec)){
        groups <- c(match(spec,
            c(attr(getNumberWays(object@poll),"ways"),"simnum" )  ) ,
                     getNumberWays(object@poll)[1] +1)
      } else {
        groups <- c(which (spec == TRUE),  getNumberWays(mrp.simple@poll)[1] +1)
      }
      if (length(groups) == 0) {
        return (sum (poststratified.sim, na.rm=TRUE) / sum (pop.sim))
      } else {
        ans <- (apply (poststratified.sim, groups, sum, na.rm=TRUE) /
              apply(pop.sim, groups, sum))
        ans[is.nan(ans)] <- NA
        return(ans)
      }
}
