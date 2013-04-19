##' Multilevel regression and poststratification
##' 
##' Set up survey and population data, and a multilevel regression model used
##' for poststratifying by an arbitrary number of strata or \dQuote{ways}.
##' 
##' 
##' @param formula.cell A formula representation of the binary outcome variable
##' and the desired eventual poststratification; i.e., the \dQuote{ways} by
##' which to break down the poll and population data, which should always be
##' given as factor variables with matching names and levels in both
##' \code{data} and \code{population}. By default, this formula will also be
##' used to construct the multilevel model, with an intercept estimated for
##' each named stratum. See \code{formula} below to easily fit more complex
##' models.
##' @param data A \code{data.frame} representing a survey, containing (at
##' least) the named variables in \code{formula.cell}.  Those variables should
##' be \code{\link[base]{factor}}s.  The LHS response is expected to be
##' dichotomous, and will be coerced to binary-logical (if factor, 1 for
##' \sQuote{yes}, 0 for \sQuote{no}).
##' @param population A \code{data.frame} containing population (e.g. census)
##' data with variable names and factor levels matching those in \code{data}
##' and specified in \code{cell.formula}.  \emph{Note: } As in \code{data}, the
##' cell formula variables should be of type \code{\link[base]{factor}}.
##' @param pop.weights character. The column of the \code{population} data that
##' contains frequencies or proportions (of the entire population) for the
##' cells defined in \code{formula}.
##' @param grouplevel.data.frames A \code{list} of \code{data.frame}s to be
##' left-joined onto the data via \code{\link[plyr]{join}}.  An example is
##' mrp.regions, which contains two columns, \sQuote{state} (the matched key)
##' and \sQuote{region}, a column being added. Multiple keys (e.g., age and
##' education) are supported.
##' @param grouplevel.expressions A \code{list} of \code{expression}s to be
##' evaluated in the data (with the grouplevel.data.frames already joined). In
##' the example below we construct two types: interaction terms and a linear
##' prior mean for a factor: \code{expression(interaction1and2 <-
##' interaction(term1,term2))} \code{expression(z.age <- rescale(age))}
##' 

##' @param formula.model.update A formula specification for the multilevel
##' model to run in the prepared data. The left-hand side should always be
##' \sQuote{response}. For convenience, the formula is handled by
##' \code{update.formula} so that \code{.} indicates the current formula
##' contents on either side of the \code{~}, e.g., \code{.~.+newVar}. The
##' initial default formula is constructed as just an intercept term for each
##' of the variables in the main formula specification
##' (\code{(1|way1)+(1|way2)} etc.)
##' @param formula.pop.update Any modifications to be made to the
##' \code{formula.cell} above. In the example below, we control for poll, but
##' the population is the same for all values of \emph{poll}.  If used, should
##' be of the \code{\link[stats]{update.formula}} template form,
##' \code{.~.-var}.
##' 
##' This replicates the population data for all levels of the variable
##' \sQuote{excluded} in this fashion.  \emph{Note:} This argument \bold{will
##' change}, to \code{constant.population.dimensions="chr"}.
##' @param poll.weights Name of variable of the survey weights for respondents
##' in the poll. This is used to compute the effective \eqn{N}, weighted
##' \eqn{\bar{Y}}, and Design Effect. Default is to make all weights equal.
##' 
##' Ideally, the dimensions specified by \code{formula.cell} would account for
##' all of the variation of survey weights in respondents in all cells.
##' Sometimes, survey researchers design samples that leave some variation in
##' the design weights of poststratification cells. Using \code{poll.weights}
##' inflates or deflates the effective \eqn{N} of respondents in
##' poststratification cells based on the average variance of design weights in
##' all cells and each cell's deviation from that overall design effect.
##' 
##' If multiple polls are included and contain poll.weights, they must be
##' normalized within each poll before \emph{mrp} attempts to normalize the
##' weights across all polls for all cells.
##' @param pop.margin Margin of population data on which to sum other cells.
##' Used in modeling PartyID, but implementation is not stable.
##' @param \dots Additional arguments to be passed to the multilevel regression
##' \code{\link[blme]{bglmer}} step.
##' @seealso \code{\link{mrp-class}} for other methods on the objects produced
##' by \code{mrp()}; \code{\link{plotmrp}} for how to plot poststratified
##' results onto maps.
##' @examples
##' 
##' \donttest{
##' library(mrpdata)
##' library(mrp)
##' 
##' ## Load example data.
##' data(CCES.complete)
##' 
##' ## Helper datasets for other US applications of MRP:
##' data(spmap.states) # projected US state map
##' data(mrp.census)   # census with common demo strata
##' data(mrp.regions)  # regions data.frame with DC separate
##' 
##' ## To ensure matching of strata between poll and population,
##' ## both should be factors with identical names and levels.
##' CCES.complete <- within (CCES.complete, {
##'   education <- factor(education,exclude=NA)
##'   female <- factor(sex=="Female", labels=c("Male","Female"), exclude=NA)
##'   race <- factor(race,exclude=NA) 
##'   f.race <- interaction(female,race)
##' })
##' 
##' ## Poll has four levels of education, so we need to combine
##' ## the top two levels in the census data. We'll also go ahead
##' ## and trim it down to just the variables used here.
##' 
##' mrp.census <- within(mrp.census,{
##'     age <- factor(age,exclude=NA,labels=c("18-29","30-44","45-64","65+"))
##'     education[education=="postgraduate"] <- "college graduate"
##'     education <- factor(education,exclude=NA)
##'     edu <- factor(education,exclude=NA,labels=c("< High School",
##'                                          "High School",
##'                                          "Some College",
##'                                          "Graduated College"))
##'     state <- factor(state,exclude=NA)
##'     race <- factor(race,exclude=NA)
##'     f.race <- interaction(sex,race)
##' })
##' mrp.census <- na.omit(mrp.census)
##' 
##' ## Ready to run simple mrp with poll and population:
##' mrp.simple <- mrp(ban.gaymarr ~ state+age+education+race, 
##'                   data=CCES.complete,
##'                   population=mrp.census,
##'                   pop.weights="weighted2004")
##' print(100*poststratify(mrp.simple, ~ education+age), digits=2)
##' \dontrun{
##' ## Fit a fuller model, adding state-level predictors:
##' ## This model is also used in the not-run example
##' ## for plotting maps.
##' mrp.statelevel <- mrp(ban.gaymarr~
##'                       state+f.race+age+education,
##'                       data=CCES.complete,
##'                       population=mrp.census, pop.weights="weighted2008",
##'                       population.formula.update= .~.-age,
##'                       grouplevel.data.frames=list(Statelevel,
##'                         mrp.regions),
##'                       grouplevel.expressions=list(
##'                         expression(age.edu <- interaction(age,education)),
##'                         ## an ordered factor, we use a normalized
##'                         ## continuous z.age as the prior mean for
##'                         ## varying intercepts of the 'age' groups.
##'                         ## That is, the prior mean for
##'                         ## age cat 1 of 4 (18-29) becomes (-.58)
##'                         expression(z.age <- rescale(age)))
##'                       )
##' ## Note: the formula is expanded from the condensed version in "formula" to
##' ##  an expanded version.
##' getFormula(mrp.statelevel)
##' 
##' ## Update the model.formula on already-prepared mrp object and re-fit:
##' mrp.statelevel <- mr(mrp.statelevel, .~.+(1|region)+ (1|age.edu)+
##'                      z.age+p.relig.full+p.kerry.full)
##' 
##' ## Fine plot control is shown with this example in plotmrp documentation!
##' }
##' }
##' 
mrp <- function(formula.cell,
                data,
                population=NULL,
                pop.weights=NULL,
                grouplevel.data.frames=NULL,
                grouplevel.expressions=NULL,
                formula.model.update=NULL,
                poll.weights=1,
                formula.pop.update=formula.cell,
                pop.margin=NULL,
                ...) {
    poll <- data ## 'data' later becomes the binomial form
    ## geographic-demographic data.frame, but
    ## is the natural argument name in the function.
    pop <- population
    formula.cell <- as.formula(formula.cell)
    mrp.terms <- terms(formula.cell)
    mrp.varnames <- attr(mrp.terms,"term.labels")
    population.formula <- update(formula.cell, formula.pop.update)
    population.terms <- terms(population.formula)
    population.varnames <- attr(terms(population.formula),"term.labels")
    population.varnames <- reorder.popterms(mrp.varnames, population.varnames)
    response.varname <- as.character (formula.cell[[2]])
    response <- poll[, response.varname]
    response <- checkResponse(response, response.varname)

    allvars <- all.vars(formula.cell)
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
                           response=as.character(formula.cell[[2]]),
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

        if(population.formula != formula.cell) {
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
    data$finalrow <- 1:nrow(data)

    ## Attempt merges. ##
    if(length(grouplevel.data.frames)>0){
        for(d in 1:length(grouplevel.data.frames)){
            data <- join(data,grouplevel.data.frames[[d]], type="left")
        }
    }
    if(length(grouplevel.expressions)>0){
        data <- within(data, sapply(grouplevel.expressions, eval.parent, n=2))
    }

    ## build the default formula unless one has been supplied
    mr.f <- formula(paste("response ~",
                          paste(paste("(1|",
                                      mrp.varnames,")"),
                                collapse="+"))
                    )
    if (!missing(formula.model.update)){
        mr.f <- update.formula(mr.f, formula.model.update)
    }
    mrp <- new("mrp",
               poll=poll.array,
               data=data,
               formula=mr.f,
               population=pop.array,
               outcome=as.character(formula.cell[[2]])
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
        ## include the factor levels in this warning.
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
##' Function to paste together a list of items, separated by commas (if more
##' than 2), and with the last one having the collapse string.
##' 
##' 
##' @param x vector or list
##' @param collapse default="and"
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
