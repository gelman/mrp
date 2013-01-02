

mrp <- function(formula,
                data, poll.weights=1,
                population=NULL,
                pop.weights=NULL,
                pop.margin=NULL,
                population.formula=formula,
                add=NULL, mr.formula=NULL,
                ...) {
    poll <- data ## 'data' later becomes the binomial form
    ## geographic-demographic data.frame, but
    ## is the natural argument name in the function.
    pop <- population
    mrp.formula <- as.formula(formula)
    mrp.terms <- terms(mrp.formula)
    mrp.varnames <- attr(mrp.terms,"term.labels")
    population.formula <- update(mrp.formula, population.formula)
    population.terms <- terms(population.formula)
    population.varnames <- attr(terms(population.formula),"term.labels")
    population.varnames <- reorder.popterms(mrp.varnames, population.varnames)

    response <- poll[, as.character (formula[[2]])]
    response <- checkResponse(response)

    allvars <- all.vars(mrp.formula)
    if(poll.weights!=1){ allvars <- c(allvars,poll.weights) }

    ## for complete, imputed datasets don't drop any columns.
    ## Useful especially when continuous imputed rather than rescaled
    ## below as an 'add'.
    if(any(is.na(poll))) {
        poll <- na.omit(poll[,allvars])
    }
    poll <- checkPoll(poll)

    ## Set up and store poll NWayData
    cat("\nExpanding data to array:\n")
    if (sum(mrp.varnames %in% names(poll)) != length(mrp.varnames) ) {
        stop(paste("\nVariable ",sQuote(mrp.varnames[!(mrp.varnames %in% names(poll))]),
                   " not found in poll data."))
    }

    ## 1. make a base data array with levels contained in the data
    ## 2. make a population array with the full set of levels
    ## 3. for levels in population NOT in poll,
    ##    fill [yes,no,N.eff] with 0.
    ## 4. flatten and then do joins and expressions.
    poll.array <- NWayData(df=poll, variables=mrp.varnames,
                           response=as.character(mrp.formula[[2]]),
                           weights=poll.weights, type="poll")

    if(is.data.frame(pop)) {
        cat("\nMatching poll data to population cells.\n")
        checkPopulationData(pop, population.varnames)
        pop.array <- makePopulationArray(pop, pop.weights, population.varnames,
                                         pop.margin=pop.margin)


        pop.subscripts <- lapply(population.varnames$inpop,
                                       findBsubscriptsInA,
                                       A=pop.array, B=poll.array)
        pop.dimnames <- dimnames(pop.array)

        if(population.formula != formula) {
            addTheseSubscripts <- lapply(population.varnames$notinpop,
                                         addSubscriptsForPollControls,
                                         poll.array=poll.array)
            addTheseDimnames <- lapply(population.varnames$notinpop,
                                       addDimnamesForPollControls,
                                       poll.array=poll.array)
            names(addTheseDimnames) <- population.varnames$notinpop
            pop.subscripts <- c(pop.subscripts, addTheseSubscripts)
            pop.dimnames <- c(pop.dimnames, addTheseDimnames)
        }
        pop.subscripts <- na.omit(as.matrix(expand.grid(pop.subscripts)))
        colnames(pop.subscripts) <- names(pop.dimnames)
        pop.array <- array(pop.array,
                           dim=lapply(pop.dimnames, length),
                           dimnames=pop.dimnames)
        poll.array <- expandPollArrayToMatchPopulation(poll.array, pop.array,
                                                       pop.subscripts)
        pop.array <- new("NWayData", pop.array, type="population",
                         levels=dimnames(pop.array))
    } else { ## No population supplied
        pop.array <- makeOnesNWay(poll.array)
    }



    cat("\nCondensing full data array to matrix for modeling:\n")
    data <- NWayData2df (poll.array)
    data.expressions <- as.expression(add[sapply(add, is.expression)])
    data.merges <- add[sapply(add, is.data.frame)]
    data$finalrow <- 1:nrow(data)

    if(length(data.expressions)>0){
        data <- within(data, sapply(data.expressions, eval.parent, n=2))
    }
    ## Attempt merges. ##
    if(length(data.merges)>0){
        for(d in 1:length(data.merges)){
            data <- join(data,data.merges[[d]], type="left")
        }
    }

    ## build the default formula unless one has been supplied
    mr.f <- formula(paste("response ~",
                          paste(paste("(1|",
                                      mrp.varnames,")"),
                                collapse="+"))
                    )
    if (!missing(mr.formula)){
        mr.f <- update.formula(mr.f, mr.formula)
    }
    mrp <- new("mrp",
               poll=poll.array,
               data=data,
               formula=mr.f,
               population=pop.array,
               outcome=as.character(formula[[2]])
               )

    cat("\nRunning Multilevel Regression step.\n")
    response <- as.matrix(getResponse(mrp))
    try(mrp <- mr(mrp,
                  ## blmer options here, possibly moved to blmer defaults
                  ...))
    return(mrp)
}

checkResponse <- function(response, varname) {
    if(is.ordered(response) && length(levels(response))==2){
        warning("Assuming ordered factor 2 levels represent 1=FALSE, 2=TRUE\n")
        response <- as.integer(response)-1
    }
    if (!is.numeric (response)) {
        stop (paste0(sQuote(varname),
                     " must be integer values of 0 / 1 or logical"))
    }
    if (length (unique(na.exclude(response))) != 2) {
        stop (paste0(sQuote(varname),
                     " must have two values"))
    }
    if (all (c(0, 1) != sort(unique(na.exclude(response))))) {
        stop (paste0(sQuote(varname), " has values of ",
                     sort(unique(na.exclude(response)))))
    }
    if(is.logical(response)) {
        response <- as.integer(response)
    }
    response
}
checkPoll <- function(poll){
    return(poll)
}
makePopulationArray <- function(pop, pop.weights, population.varnames,
                                pop.margin=pop.margin) {

    if (is.null(pop.weights)) {
        warning("Warning: pop.weights unspecified. Assuming all population cells of equal weight.") }
    main.pop.formula <- paste0(pop.weights, "~",
                                            paste(population.varnames$inpop,
                                                  collapse="+"))
    ## xtabs are arrays formed using formula interface
    ## prop.table with no 'margin'
    pop.array <- prop.table(xtabs(main.pop.formula, data=pop),
                            margin=pop.margin)
    pop.array
}

  ## For population array, if there are "ways" present in poll but constant
## in population, move those terms to the end. Will become constant (1s).
reorder.popterms <- function(poll, pop){
    inpop <- poll[poll %in% pop]
    notinpop <- poll[!{poll %in% pop}]

    return(list(inpop=inpop, notinpop=notinpop))
}

findBsubscriptsInA <- function(dim, A, B) {
    match(dimnames(B)[[dim]], dimnames(A)[[dim]])
}

checkPopulationData <- function(population.varnames, pop) {
    if (sum(population.varnames$inpop %in% names(pop)) !=
                length(population.varnames$inpop) ) {
                stop(paste("\nVariable ",
                           sQuote(
                           population.varnames$inpop[!(population.varnames$inpop
                                                       %in% names(pop))]),
                           " not found in population."))
            }
}
addSubscriptsForPollControls <- function(var, poll.array){
    1:length(dimnames(poll.array)[[var]])
}
addDimnamesForPollControls <- function(var, poll.array){
    dimnames(poll.array)[[var]]
}
expandPollArrayToMatchPopulation <- function(poll.array, pop.array,
                                             populationSubscripts){
    out.dims <- c(3, dim(pop.array))
    out.dimnames <- c(list(cellSummary=c("N", "design.effect.cell", "ybar.w")),
                         dimnames(pop.array))
    warnAboutMissingCells(dimnames(poll.array), dimnames(pop.array))

    ## fill all cells as though empty
    out <- array(c(0,1,.5), dim=out.dims, dimnames=out.dimnames)
    ## put in expected order
    out <- aperm(out, c(2:length(dim(out)),1))
    poll.array <- aperm(poll.array, c(length(dim(poll.array)),
                                      seq_len(length(dim(pop.array)))))
    poll.matrix <- matrix(poll.array, nrow=dim(populationSubscripts),
                          ncol=3, byrow=TRUE)
    ## should be able to fix this in makeNWay which otherwise is fine
    ## but empty cells should be (0,1,.5)
    poll.matrix <- t(apply(poll.matrix,1, fillNAs))
    colnames(poll.matrix) <- c("N", "design.effect.cell", "ybar.w")
    ## fill with poll data where it exists
    for(i in 1:3) {
        indToInsertFromPoll <- cbind(populationSubscripts, i)
        out[indToInsertFromPoll] <- poll.matrix[,i]
    }
    out <- new("NWayData", out, type="poll",
               levels=dimnames(pop.array))
    out
}
warnAboutMissingCells <- function(poll.dims, pop.dims) {
    for(d in names(pop.dims)) {
        if(length(poll.dims[[d]]) < length(pop.dims[[d]])) {
           warning("No data in ", sQuote(d),
                   " : ", serialPaste(setdiff(pop.dims[[d]],
                                            poll.dims[[d]])))
       }
    }
}

##' serial paste
##'
##' Function to paste together a list of items, separated by commas
##' (if more than 2), and with the last one having the collapse string.
##'
##' @param x vector or list
##' @param collapse default="and"
##' @export
serialPaste <- function (x, collapse="and") {
	##
	if (length(x)>1) x[length(x)] <- paste(collapse, x[length(x)])
	return(ifelse(length(x)>2, paste(x, collapse=", "),
                      paste(x, collapse=" ")))
}

fillNAs <- function(row) {
   if(all(is.na(row))){
       c(0,1,.5)
   } else {
       row
   }
}

## Definining Methods
## Getters and Setters
setGeneric("getData", function(object) {standardGeneric("getData")})
setMethod (f="getData",
    signature=signature(object="mrp"),
    definition=function(object) {
      return (object@data)
    })

setGeneric("getResponse", function(object) {standardGeneric("getResponse")})
setMethod( f="getResponse",
    signature=signature(object="mrp"),
    definition=function(object) {
      return(as.matrix(object@data[,c("response.yes","response.no")]))
    })

setGeneric ("getPopulation", function (object) { standardGeneric ("getPopulation")})
setMethod (f="getPopulation",
		signature=signature(object="mrp"),
		definition=function(object) {
			stopifnot (class (object) == "mrp")

			return (object@population)
		})

setGeneric ("setPopulation", function (object, population) { standardGeneric ("setPopulation")})
setMethod (f="setPopulation",
    signature=signature(object="mrp"),
    definition=function(object, population) {
      stopifnot (class (population) == "array")
      #stopifnot (dim (population) == dim (object@theta.hat))

      ## Need more checks here.
      object@population <- replace (population, is.na (population), 0)
      return (object)
    })

setGeneric ("setPopulationOnes", function (object) { standardGeneric ("setPopulationOnes")})
setMethod (f="setPopulationOnes",
    signature=signature(object="mrp"),
    definition=function(object) {
      if(object@population@type=="population") {
        warning("Population appears to contain real data. Replacing with ones!")
      }
      object@population <- makeOnesNWay(object@poll)
      return(object)

    } )

setGeneric ("setFormula", function (object, formula) { standardGeneric ("setFormula")})
setMethod (f="setFormula",
    signature=signature(object="mrp"),
    definition=function (object, formula) {
      object@formula <- formula
      return (object)
    })

setGeneric ("getFormula", function (object) { standardGeneric ("getFormula")})
setMethod (f="getFormula",
    signature=signature(object="mrp"),
    definition=function (object) {
      return (object@formula)
    })

.getThetaHat <- function(object) {
    pop <- getPopulation(object)
    theta.hat <- array (fitted(getModel(object)),
                        dim=dim(pop), dimnames=dimnames(pop))
    return(theta.hat)
}
setGeneric ("getThetaHat", function (object) { standardGeneric ("getThetaHat")})
setMethod(f="getThetaHat",signature(object="mrp"),
    definition=.getThetaHat)

setGeneric ("getEstimates", function (object) { standardGeneric ("getEstimates")})
setMethod(f="getEstimates",signature(object="mrp"),
    definition=function(object) {
      return(getThetaHat(object))
    })
setGeneric ("getModel", function (object) { standardGeneric ("getModel")})
setMethod(f="getModel",signature(object="mrp"),
    definition=function(object) {
      return(object@multilevelModel)
    })
setGeneric ("getData", function (object) { standardGeneric ("getData")})
setMethod(f="getData",signature(object="mrp"),
    definition=function(object) {
      return(object@data)
    })
setGeneric ("getOutcome", function (object) { standardGeneric ("getOutcome")})
setMethod(f="getOutcome",signature(object="mrp"),
    definition=function(object) {
      return(object@outcome)
    })
setGeneric ("getAdded", function (object) { standardGeneric ("getAdded")})
setMethod(f="getAdded",signature(object="mrp"),
    definition=function(object) {
      return(object@added)
    })





##### NEW SHIFT FUNCTION for state vote total.
##### only one margin now.

## shift(mrp,
## turnoutData (vector or char name of data col),
## ~shiftvar)

## p1 = collapse across population by shiftvar
## delta = apply (p1, shiftfun, turnoutdata)
##  do this optimization and get a bunch of
##  shiftresults
## p2 = sweep delta+population array
## return a full-dimension shifted array



.poststratify <- function (object, formula=NULL) {
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
    }
setGeneric ("poststratify", function (object, formula=NULL, ...) { standardGeneric ("poststratify")})
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


##### INTERCEPT SHIFT BY STATE FOR KNOWN TURNOUT
##    from Yair. See paper, p 12.
##    this will go in .fun arg of aaply()
################################################
##
## logit <- function (a) log(a/(1-a))
## find.delta <- function(delta, a, w, x0)
## abs(x0-sum(invlogit(logit(a) + delta)*w))
## correct.weighted <- function(a, w, x0) {
##   delta <- optimize(find.delta, interval=c(-5,5), a, w, x0)$minimum
##   corrected <- invlogit(logit(a) + delta)
##   return(list(delta=delta, corrected=corrected))
## }
## )


## newMrp <- function (response, vars, population, weight=rep(1, length(response))) {
##     # check inputs
##     if ("data.frame" != class(vars)) {
##         stop ("vars must be a data.frame.")
##     }
##     if (length(response) != nrow(vars)) {
##         stop ("response must have the same length as the vars.")
##     }
##     if (nlevels(response) != 2) {
##         stop (paste ("response must have 2 levels, found:", nlevels(response)))
##     }
##     # make inputs factors
##     response <- factor (response)
##     for (i in 1:length (vars)) {
## 		if (is.factor (vars[,i]) == FALSE) {
## 			vars[,i] <- factor(vars[,i])
## 		}
##     }

##     data <- newNWayData (ncol(vars), data.frame (response, vars, weight))
##     if (missing(population)) {
##         population <- array (1, dim(data@ybarWeighted))
##     }

##     if (all (dim(data@ybarWeighted) == dim (population)) == FALSE) {
##       stop (paste ("dim (population) must match dim (data@ybarWeighted).\n\tExpected:", dim (data@ybarWeighted)," but found: ", dim (population)))
##     }

##     formula <- paste ("cbind (response.yes, response.no) ~ 1 +",
##                       paste ("(1 | ", names (vars), ")", sep="", collapse=" + "))
##     return (new(Class="mrp", data=data, population=population, formula=formula))
##   }

##########
## TODO : vignette note: don't make dumb interactions because you're not using fucking stata
## doc note: first formula needs to be clear that these are ALL the RANDOM INTERCEPTS
## "defined at the individual level in the megapoll"
## ANY and ALL group-level terms, predictors, interactions, etc:
## are added via the add list.
## consider aliasing or renaming "add" to "grouplevel" or "group.vars"
##
## p.formula, = base cells
## intercept.formula, = p.forumla plus some stuff before response frame is made
## mr.formula = intercepts plus any from grouplevel

# intercept.formula = age + edu + income + poll, # in megapoll at indiv level.
# p.formula = age + edu + income,
## mr.formula = (1|age) + (1|edu) + (1|income) + (1|poll) + (1|region) + (1|age.edu)

## provide helper function to make a poll.data data.frame
## several responses that we want to deal with.
