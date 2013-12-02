## Uncertainty using sim() on lme4-fitted models
## Mostly by Jeffrey R. Lax (@jrllrj) and Michael Malecki @malecki
## Use at your own risk. This needs review and revision before
## it can move into mrp itself.
###################################


##random set up stuff you might need to run first
library(devtools)
install_github("mrp", "malecki", sub="mrpdata")
install_github("mrp", "malecki", sub="mrp")
install_url(url = "http://cran.r-project.org/src/contrib/Archive/blme/blme_0.01-6.tar.gz")
install_url(url = "http://cran.r-project.org/src/contrib/Archive/lme4/lme4_0.999999-2.tar.gz")
install_url(url = "http://cran.r-project.org/src/contrib/Archive/arm/arm_1.6-07.01.tar.gz")


rm(list=ls(all=TRUE))


#load packages
library(boot)
library(arm)
library(foreign)
library(car)
library(blme)
library(plyr)
library(mrp)
library(mrpdata)




#loads some data I might use
data(dpsp)
data(mrp.regions)

#set signif digits to 3 and push against using scientific notation
options(digits = 3)
options(scipen = 100)





##############  vignette

## Load example data.
data(CCES.complete)

## Helper datasets for other US applications of MRP:
data(spmap.states) # projected US state map
data(mrp.census)   # census with common demo strata
data(mrp.regions)  # regions data.frame with DC separate

## To ensure matching of strata between poll and population,
## both should be factors with identical names and levels.
CCES.complete <- within (CCES.complete, {
  education <- factor(education,exclude=NA)
  female <- factor(sex=="Female", labels=c("Male","Female"), exclude=NA)
  race <- factor(race,exclude=NA)
  f.race <- interaction(female,race)
})

## Poll has four levels of education, so we need to combine
## the top two levels in the census data. We'll also go ahead
## and trim it down to just the variables used here.

mrp.census <- within(mrp.census,{
    age <- factor(age,exclude=NA,labels=c("18-29","30-44","45-64","65+"))
    education[education=="postgraduate"] <- "college graduate"
    education <- factor(education,exclude=NA)
    edu <- factor(education,exclude=NA,labels=c("< High School",
                                         "High School",
                                         "Some College",
                                         "Graduated College"))
    state <- factor(state,exclude=NA)
    race <- factor(race,exclude=NA)
    f.race <- interaction(sex,race)
})
mrp.census <- na.omit(mrp.census)

#CCES.complete <- subset(CCES.complete, state!="AL" )
#CCES.complete <- rbind(CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete,CCES.complete, CCES.complete)

## Ready to run simple mrp with poll and population:
mrp.simple <- mrp(ban.gaymarr ~ state+age+education+race,
                  data=CCES.complete,
                  population=mrp.census,
                  pop.weights="weighted2004")
print(100*poststratify(mrp.simple, ~ education+age), digits=2)


##################

#FITTED FOR SIM (will be included in arm soon) --  from Vince Dorie


setMethod("fitted", signature(object = "sim.mer"),
          function(object, regression)
          {
            if (missing(regression) || is.null(regression)) stop("fitted for sim.mer requires original mer object as well.");
            if (!inherits(regression, "mer")) stop("regression argument for fitted on sim.mer does not inherit from class 'mer'");
            sims <- object;
            numSimulations <- dim(sims@fixef)[1];
            numRanef <- regression@dims[["q"]];
            numLevels <- regression@dims[["nt"]];

            ranefNames <- lapply(regression@ST, colnames);
            # ranefStructure <=> list w/one element per level, each is list with first element
            # indicies into joint covar, second with name of what is varying
            ranefStructure <- lme4:::plist(lme4:::reinds(regression@Gp), ranefNames);
            wt <- lme4:::whichterms(regression);
            ranefLabel <- lme4:::plist(wt, regression@flist);

            ranef <- matrix(0, numRanef, numSimulations);

            levelIndex <- 1;
            numUniqueLevels <- length(ranefLabel);

            for (uniqueLevelNum in 1:numUniqueLevels) {
              numCombinedLevels <- length(ranefLabel[[uniqueLevelNum]][[1]]);

              simIndex <- 0;
              for (subLevelNum in 1:numCombinedLevels) {
                ranefIndices <- ranefStructure[[levelIndex]][[1]];
                simIndices <- 1:length(ranefStructure[[levelIndex]][[2]]) + simIndex;

                numCoefs <- dim(sims@ranef[[uniqueLevelNum]])[2];
                for (i in 1:length(simIndices)) {
                  ranefSubset <- 1:numCoefs + (i - 1) * numCoefs;

                  ranef[ranefIndices[ranefSubset],] <- t(sims@ranef[[uniqueLevelNum]][,,simIndices[i]]);
                }

                simIndex <- simIndex + length(simIndices);
                levelIndex <- levelIndex + 1;
              }
            }

            linkType <- regression@dims[["lTyp"]];
            invLinks <- c(function(x) { e.x <- exp(x); e.x / (1 + e.x); }, # logit
                          function(x) { pnorm(x); },                       # probit
                          function(x) { pcauchy(x); },                    # cauchit
                          function(x) { 1 - exp(-exp(x)); },               # comp log-log
                          function(x) { x; },                              # identity
                          function(x) { exp(x); },                         # log
                          function(x) { x^2; },                            # sqrt
                          function(x) { x^-0.5; },                         # 1 / mu^2
                          function(x) { 1 / x; });                         # inverse
            invLink <- invLinks[[linkType]];

            result <- tcrossprod(regression@X, sims@fixef) + crossprod(regression@Zt, ranef);
            return(invLink(result));
          });




##################
# KEY NEW BIT -- MRP WITH UNCERTAINTY (Jeffrey Lax and Andrew Guess)


poststratify.sim <- function(object, formula, n.sim=2) {

  n.sim <- max(2,n.sim)
  pop <- getPopulation(object)
  resultDimnames <- dimnames(pop)
  resultDimnames$simnum <- c(1:n.sim)
  theta.hat.sim <- array (fitted(sim(getModel(object), n.sim),getModel(object) ),
                          dim=c(dim(pop),n.sim), dimnames=) # I TOOK OUT "DIMNAMESTEMP" - AG 11/2/13
  pop.sim <- pop
  pop.sim@levels$simnum <- c(1:n.sim)
  pop.sim@.Data <- array (rep(pop.sim@.Data,n.sim),
                          dim = c(dim(pop),n.sim),
                          dimnames = resultDimnames)

  #set up the array shape to hold the full post array `n.sims` times
  c(dim(pop.sim@.Data), "simnum" = n.sim)

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



#to use it:
# poststratification without uncertainty
poststratify(mrp.simple, ~education)
#poststratification with uncertainty
poststratify.sim(mrp.simple, ~  education)  #default 2 draws
poststratify.sim(mrp.simple, ~  education,3)  #now 3 draws
     poststratify.sim(mrp.simple, ~ state,5)
    poststratify.sim(mrp.simple, ~ age + education,3)
#to produce 10 estimates of overall mean
poststratify.sim(mrp.simple,~  -simnum,10)

print(100*poststratify.sim(mrp.simple, ~ education+age,3), digits=2)


poststratify(mrp.simple, ~state)


temptemp <-     poststratify.sim(mrp.simple, ~ state,500)

#can calculate stats easily but be careful
#below are only accurate if "state" is the only dimension in the poststratification --- that is, if you are simming on one variable
apply(temptemp,c("state"),mean)
apply(temptemp,c("state"),sd)
apply(temptemp,c("state"),quantile,probs = c(.025,.5,.975))




temptemp2 <-     poststratify.sim(mrp.simple, ~ education,500)

#can calculate stats easily but be careful
#below are only accurate if "education" is the only dimension in the poststratification
apply(temptemp2,c("education"),mean)
apply(temptemp2,c("education"),sd)
apply(temptemp2,c("education"),quantile,probs = c(.025,.5,.975))

#to do it by two variables, make sure two variables are included both in the poststratify.sim call AND the apply --- otherwise it weights across the missing categories equally which is unlikely to be desired
temptemp3 <-     poststratify.sim(mrp.simple, ~ education + age,500)
apply(temptemp3,c("age", "education"),mean)
apply(temptemp3,c("age", "education"),sd)
apply(temptemp3,c("age", "education"),quantile,probs = c(.025,.5,.975))




