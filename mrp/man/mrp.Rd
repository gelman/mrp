\name{mrp}
\alias{mrp}

\title{Multilevel regression and poststratification}
\description{
  Set up survey and population data, and a multilevel regression model
  used for poststratifying by an arbitrary number of strata or \dQuote{ways}.
}

\usage{
mrp (formula,
     data, poll.weights=1,
     population=NULL, use=NULL,
     population.formula=formula,
     add=NULL, mr.formula=NULL,
     \dots)
}
\arguments{
    \item{formula}{A formula representation of the binary outcome variable
      and the desired eventual
      poststratification; i.e., the \dQuote{ways} by which to
      break down the poll and population data, given as factor
      variables with matching names and levels in both \code{poll}
      and \code{population}. By default, this formula will also be
      used to construct the multilevel model, with an intercept
      estimated for each named stratum. See \code{mr.formula}
      below to easily fit more complex models.}
    \item{data}{A \code{data.frame} representing a survey, containing (at
      least) the named variables in \code{formula}. The LHS response
      is expected to be dichotomous, and will be coerced to
      binary-logical (if factor, 1 for \sQuote{yes}, 0 for
      \sQuote{no}).}
    \item{poll.weights}{Name of variable of the survey weights for 
	  respondents in the poll. This is used to compute the effective 
	  \eqn{N}, weighted \eqn{\bar{Y}}, and Design Effect. Default is 
	  to make all weights equal.}
    \item{population}{A \code{data.frame} containing population (e.g. census)
      data with variable names and factor levels matching those in
      \code{poll} and specified in \code{formula}.}
    \item{use}{The column of the \code{population} data to use for
      poststratification.}
    \item{population.formula}{Any modifications to be made to the
      formula above. In the example below, we poststratify by poll,
      which has to be removed from the ways that we stratify the
      population. If used, should be of the form
      \code{.~.-var}. 
	  % DL: Should remove? MM: yes
	  % It is possible, but strongly discouraged, to use a matching
      % factor with a different name (e.g., \sQuote{edu} in \code{poll}
      % with \code{population.formula=.~education})
	  }
    \item{add}{A \code{list} of additional variables to add to
      the prepared \code{data.frame} used in the multilevel
      regression step. The list may contain items of two types of
      objects: a \bold{data.frame} to be left-joined onto the data via
      \code{merge}; or an \bold{expression} to be evaluated \sQuote{within} 
	  the prepared data, such as a transformation of another variable. 
	  An example of the latter would be a data.frame with
      a column containing the \sQuote{region} for every \sQuote{state}
      in the data.}
    \item{mr.formula}{A formula specification for the multilevel model
      to run in the prepared data. The left-hand side should always be
    \sQuote{response}. For convenience, the formula is handled by
      \code{update.formula} so that \code{.} indicates the 
    current formula contents on either side of the \code{~}, e.g.,
    \code{.~.+newVar}. The initial default formula is constructed as
      just an intercept term for each of the variables in the main
      formula specification (\code{(1|way1)+(1|way2)} etc.)}
    \item{\dots}{Additional arguments to be passed to the multilevel
      regression \code{\link[blme]{bglmer}} step.}
  }
\seealso{\code{\link{mrp-class}} for other
    methods on the objects produced by \code{mrp()};
      \code{\link{plotmrp}} for how to plot poststratified results  onto maps.}
    \examples{
\donttest{
library(mrpdata)
library(mrp)

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
  female <- factor(sex=="Female",exclude=NA)
  race <- factor(race,exclude=NA)
  f.race <- interaction(female,race)
})

## Poll has four levels of education, so we need to combine
## the top two levels in the census data. We'll also go ahead
## and trim it down to just the variables used here.
mrp.census <- na.omit(mrp.census[mrp.census$state \%in\% CCES.complete$state,])
mrp.census <- na.omit(mrp.census[mrp.census$education \%in\% CCES.complete$education,])

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

## Ready to run simple mrp with poll and population:
mrp.simple <- mrp(ban.gaymarr ~ state+age+education, 
                  data=CCES.complete,
                  population=mrp.census,
                  use="weighted2004",cov.prior="none")
print(100*poststratify(mrp.simple, ~ education+age), digits=2)
\dontrun{
## Fit a fuller model, adding state-level predictors:
## This model is also used in the not-run example
## for plotting maps.
mrp.statelevel <- mrp(ban.gaymarr~
                      state+f.race+age+education,
                      data=CCES.complete,
                      population=mrp.census,use="weighted2008",
                      population.formula=.~.-poll,
                      add=list(Statelevel,
                        mrp.regions,
                        expression(age.edu <- interaction(age,education)),
                        expression(z.age <- rescale(age)))
                      )
## Note: the formula is expanded from the condensed version in "formula" to
##  an expanded version.
getFormula(mrp.statelevel)

## Update the mr.formula on already-prepared mrp object and re-fit:
mrp.statelevel <- mr(mrp.statelevel, .~.+(1|region)+ (1|age.edu)+
                     z.age+p.relig.full+p.kerry.full)

## Fine plot control is shown with this example in plotmrp documentation!
}
}
}
