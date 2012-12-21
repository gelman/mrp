setClass("NWayData",representation(type="character",levels="list"),contains="array")
setClass("jagsNWayData",representation(type="character",levels="list"),contains="array")

## save the levels on the original data for when it is plyd back
## in poststratification
## match on names of ways
saveNWayLevels <- function(df, variables=TRUE){
  fac <- sapply(df[variables],is.factor)
  lev <- lapply(df[,fac],attributes)
  return(lev[variables])
}
restoreNWayLevels <- function(df=df,nway=nway){
  pos <- na.omit(names(nway@levels)[match(names(df),names(nway@levels))])
  df[,pos] <- as.data.frame(lapply(pos, function(col) {
                df[,col] <- factor(df[,col],
                levels=nway@levels[[col]]$levels,
                ordered=nway@levels[[col]]$class[1]=="ordered")
            }))
  return(df)
}

## Returns the number of ways of the analysis and an attribute "ways"
## a character vector of names of 'ways' variables.
setGeneric ("getNumberWays", function (object) { standardGeneric ("getNumberWays") })
setMethod (f="getNumberWays",
    signature=signature(object="NWayData"),
    definition=function(object){
      ## Poll has an extra dimension containing
      ## computed values (ybar, N, design.effect.cell)
      is.poll <- ifelse(attr(object,"type")=="poll",
          TRUE,FALSE)
      w <- length(attr(object,"dim"))-is.poll
      attr(w,"ways") <- names(dimnames(object))[1:w]
      if(object@type=="ones") { w <- 0 }
      return(w)
    })

setGeneric ("getYbarWeighted", function (object) { standardGeneric ("getYbarWeighted")})
setMethod (f="getYbarWeighted",
    signature=signature(object="NWayData"),
    definition=function (object) {
      ybar.w <- do.call("[",c(x=quote(object),
              as.list(rep(TRUE,getNumberWays(object))),
              "ybar.w"))
      return(ybar.w)
    })

setGeneric ("getN", function (object) { standardGeneric ("getN")})
setMethod (f="getN",
    signature=signature(object="NWayData"),
    definition=function(object){
        N <- do.call("[",c(x=quote(object),
              as.list(rep(TRUE,getNumberWays(object))),
              "N"))
        N[is.na(N)] <- 0
        return(N)
    }  )

setGeneric ("getDesignEffect", function (object) { standardGeneric ("getDesignEffect")})
setMethod (f="getDesignEffect",
    signature=signature(object="NWayData"),
    definition=function (object) {
      D <- getDesignEffectByCell(object)
      N <- getN(object)
      #subsY <- {length(subsN)*2+1}:length(array)
      design.effect <- weighted.mean(D,N,na.rm=TRUE)
      return(c("design.effect"=design.effect))
    })


setGeneric ("getDesignEffectByCell", function (object) { standardGeneric ("getDesignEffectByCell")})
setMethod (f="getDesignEffectByCell",
    signature=signature(object="NWayData"),
    definition=function (object) {
      D <- do.call("[",c(x=quote(object),
              as.list(rep(TRUE,getNumberWays(object))),
              "design.effect.cell"))
      return(D)
    })


## Returns the number of observations (unweighted observations) that was used to create the data set.
## TODO: remove this method -- not useful. Perhaps a weighted number of observations would be more useful
## setGeneric ("getDataLength", function (object) { standardGeneric ("getDataLength")})
## setMethod (f="getDataLength",
##         signature=signature(object="NWayData"),
##         definition=function (object) {
##             return (object@dataLength)
##         })

setGeneric ("getNEffective", function (object) { standardGeneric ("getNEffective")})
setMethod (f="getNEffective",
    signature=signature(object="NWayData"),
    definition=function(object){
      getN(object) / getDesignEffect(object)
    })

setGeneric ("getDesignEffect", function (object) { standardGeneric ("getDesignEffect")})
setMethod (f="getDesignEffect",
    signature=signature(object="NWayData"),
    definition=function (object) {
      return (weighted.mean (getDesignEffectByCell (object), getN(object), na.rm=TRUE))
    })

setGeneric ("getData", function (object) { standardGeneric ("getData")})
setMethod (f="getData",
    signature=signature(object="NWayData"),
    definition=function (object) {
      return (object@data)
    })

### makeNway is meant to be called on a sliced (subsetted) data.frame
### made by the call to plyr:::daply (dataframe-to-array)
### args: response - column name of binary response
###       weights - column name of survey weight var
###       pop - logical, for population data (just sums)
### when used for population data, must supply 'weights'
setGeneric ("makeNWay", function (cell,response,weights,pop) { standardGeneric ("makeNWay")})
setMethod (f="makeNWay",
    signature=signature(cell="data.frame"),
    definition=function(cell, response="response",
        weights=1, pop=FALSE) {
      if(pop==TRUE){
        if(length(weights)!=1) {
          stop(paste("When supplying a data.frame as ",sQuote("population")," you must also indicate ",sQuote("use"), ", which column of the data.frame to use.\n")) }
        prop <- sum(cell[,weights])
        return(prop)
      }

      N <- nrow(cell)
      y <- cell[,response]
      ## quietly allow easy noweight
      if(weights==1) {
        cell$weight <- rep(1,nrow(cell))
        weights <- "weight"
      }
      w <- cell[,weights]
      ## do weighted mean
			ybar.w <- weighted.mean(y, w)

      if( N > 1 & all(w==1) ) {
				# This is just for speed. Can use the other loop without
				# an if statement and it will produce the same result.
				design.effect.cell <- 1
      } else {
        design.effect.cell <- 1+ var(w/mean(w))
      }
      ybar.w[is.nan(ybar.w)] <- 0

      ans <- c(N=N,
          design.effect.cell=design.effect.cell,
          ybar.w=ybar.w)
      #print(str(ans))
      return(ans)
    })

### turns a NWayData array back into a data.frame for lmer call.
### args: v, a vector of the plyr:::adply -sliced NWayData array;
###       design.effect, the averaged design effect
setGeneric ("flattenNWay", function (v,design.effect) { standardGeneric ("flattenNWay")})
setMethod (f="flattenNWay",
    definition=function(v,design.effect){
      if(is.na(v["N"])) {
        v["N"] <- 0
        v["design.effect.cell"] <- NA
        v["ybar.w"] <- .5
        response.yes <- response.no <- 0
      } else {
        ## do n.eff
        N.eff <- v["N"] / design.effect

        ybar.w <- v["ybar.w"]
        ## do ybar.w with cases
        response.yes <- ybar.w*N.eff
        response.no <- (1-ybar.w)*N.eff
      }
      ans <- c(response.yes,response.no,v)
      names(ans)[1:2] <- c("response.yes","response.no")
      return(ans)
    })

setGeneric ("makeOnesNWay", function (object) { standardGeneric ("makeOnesNWay")})
setMethod (f="makeOnesNWay",
    signature=signature(object="NWayData"),
    definition=function(object) {
      pop.nway <-  array (1,
          dim(getYbarWeighted(object)),
          dimnames=dimnames(getYbarWeighted(object)))
      pop.nway <- new("NWayData",pop.nway,type="ones",
          levels=object@levels)
      return(pop.nway)
    } )

## Convenience
is.NWayData <- function(object) {
  inherits(object,"NWayData")
}

NWayData <- function (df, variables, response, weights, type="poll", reference.poll) {
  if (type=="poll") {
    nway <- daply(df, .variables=variables, pop=FALSE,
        .fun=makeNWay, .progress="text",
        response=response, weights=weights)
  } else if (type == "population") {
    nway <- daply(df, .variables=variables, pop=TRUE,
        .fun=makeNWay, .progress="text",
        weights=weights)
    nway <- array(rep(nway,
            length.out=length(reference.poll)),
                  ## this dim/dimnames should be "the last dimension"
            dim(getNEffective(reference.poll)), dimnames(getNEffective(reference.poll)))
  }
  nway <- new("NWayData", nway, type=type,
      levels=saveNWayLevels(df, variables))

  return (nway)
}

NWayData2df <- function (nway) {
  data <- adply(nway, .margins=1:getNumberWays(nway),
      flattenNWay,
      design.effect=getDesignEffect(nway),
                .progress="text")
  data <- restoreNWayLevels(df=data,nway=nway)

  return (data)
}

## Transform for fitting in Stan, etc.

"write.stan" <- function(data, x, y, fileprefix="stan") {
  data <- as.data.frame(sapply(data, function(col) {
    if(is.factor(col) & length(levels(col)) == 2) {
      col <- as.integer(col) - 1
    }
    if(is.factor(col) & length(levels(col)) > 2) {
      col <- as.integer(col)
    }
    return(col)
  }))
  y <- data[,y]
  x <- data[,x]
  x$intercept <- 1
  write.table(x, file=paste(fileprefix,"x.dat",sep=""),
              row.names=FALSE, col.names=FALSE)
  write.table(y,file=paste(fileprefix,"y.dat",sep=""),
              row.names=FALSE, col.names=FALSE)
  data <- cbind(y,x)
  invisible(data)
}

## just print the array when asked interactively
setMethod(show, "NWayData",
          definition=function(object) show(object@.Data))
#setOldClass("array")

## set up 'sweep' on poststratified NWayData
setGeneric ("sweep")
setMethod(sweep, "NWayData",
          definition=function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
            ans <- base::sweep(x, MARGIN, STATS, FUN, check.margin = check.margin, ...)
            ans <- new("NWayData", ans, type="sweep", levels=x@levels)
            return(ans)
          })



jags2NWay <- function (df, variables, response) {
  nway <- daply(df, .variables=variables,
                .fun=makeJagsNWay, .progress="text",
                response=list(party=c("R","I","D")))
  d <- dimnames(nway)
  names(d)[length(d)] <- names(response)
  dimnames(nway) <- d
  d <- list(levels=response[[1]], class="factor")
  ## TODO remove `[variables]` when package reinstalled!
  l <- c(saveNWayLevels(df,variables)[variables],new=list(d))
  names(l) <- c(names(l)[1:length(variables)],names(response))

  nway <- new("NWayData", nway, type="population",
              levels=l)
  return (nway)
}

setGeneric ("makeJagsNWay", function (cell,response) {
  standardGeneric ("makeJagsNWay")
})
setMethod (f="makeJagsNWay",
           signature=signature(cell="data.frame"),
           definition=function(cell, response=response) {
                                        #N <- nrow(cell)
             ## apply over new dims (such as est party)
             if(nrow(cell)==1) {
               ans <- cell[, response[[1]] ]
             } else {
               ans <- sapply(response[[1]], function(d) {
                 mean(cell[,d])
               })
             }

           names(ans) <- response[[1]]
           return(ans)
           })

### these are convenient constructor for arrays, that won't break on an already-nway
setGeneric("newNWayData", function(object=NULL, type=NULL, levels=NULL) { standardGeneric("newNWayData")} )
setMethod(f="newNWayData",
          signature=signature(object="array"),
          definition=function(object, type="generic", levels=NULL){
            if(missing(levels)) {
              levels <- dimnames(object)
            }
            new("NWayData", object, type=type, levels=levels)
          })
setMethod(f="newNWayData",
          signature=signature(object="NWayData"),
          definition=function(object,type=object@type,levels=object@levels){
            new("NWayData", object, type=type, levels=levels)
          })
