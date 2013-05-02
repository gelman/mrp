\name{mrp}
\alias{mrp}
\title{Multilevel regression and poststratification}
\usage{
  mrp(formula.cell, data, population = NULL,
    pop.weights = NULL, grouplevel.data.frames = NULL,
    grouplevel.expressions = NULL,
    formula.model.update = NULL, poll.weights = 1,
    formula.pop.update = formula.cell, pop.margin = NULL,
    ...)
}
\arguments{
  \item{formula.cell}{A formula representation of the
  binary outcome variable and the desired eventual
  poststratification; i.e., the \dQuote{ways} by which to
  break down the poll and population data, which should
  always be given as factor variables with matching names
  and levels in both \code{data} and \code{population}. By
  default, this formula will also be used to construct the
  multilevel model, with an intercept estimated for each
  named stratum. See \code{formula} below to easily fit
  more complex models.}

  \item{data}{A \code{data.frame} representing a survey,
  containing (at least) the named variables in
  \code{formula.cell}.  Those variables should be
  \code{\link[base]{factor}}s.  The LHS response is
  expected to be dichotomous, and will be coerced to
  binary-logical (if factor, 1 for \sQuote{yes}, 0 for
  \sQuote{no}).}

  \item{population}{A \code{data.frame} containing
  population (e.g. census) data with variable names and
  factor levels matching those in \code{data} and specified
  in \code{cell.formula}.  \emph{Note: } As in \code{data},
  the cell formula variables should be of type
  \code{\link[base]{factor}}.}

  \item{pop.weights}{character. The column of the
  \code{population} data that contains frequencies or
  proportions (of the entire population) for the cells
  defined in \code{formula}.}

  \item{grouplevel.data.frames}{A \code{list} of
  \code{data.frame}s to be left-joined onto the data via
  \code{\link[plyr]{join}}.  An example is mrp.regions,
  which contains two columns, \sQuote{state} (the matched
  key) and \sQuote{region}, a column being added. Multiple
  keys (e.g., age and education) are supported.}

  \item{grouplevel.expressions}{A \code{list} of
  \code{expression}s to be evaluated in the data (with the
  grouplevel.data.frames already joined). In the example
  below we construct two types: interaction terms and a
  linear prior mean for a factor:
  \code{expression(interaction1and2 <-
  interaction(term1,term2))} \code{expression(z.age <-
  rescale(age))}}

  \item{formula.model.update}{A formula specification for
  the multilevel model to run in the prepared data. The
  left-hand side should always be \sQuote{response}. For
  convenience, the formula is handled by
  \code{update.formula} so that \code{.} indicates the
  current formula contents on either side of the \code{~},
  e.g., \code{.~.+newVar}. The initial default formula is
  constructed as just an intercept term for each of the
  variables in the main formula specification
  (\code{(1|way1)+(1|way2)} etc.)}

  \item{formula.pop.update}{Any modifications to be made to
  the \code{formula.cell} above. In the example below, we
  control for poll, but the population is the same for all
  values of \emph{poll}.  If used, should be of the
  \code{\link[stats]{update.formula}} template form,
  \code{.~.-var}.

  This replicates the population data for all levels of the
  variable \sQuote{excluded} in this fashion.  \emph{Note:}
  This argument \bold{will change}, to
  \code{constant.population.dimensions="chr"}.}

  \item{poll.weights}{Name of variable of the survey
  weights for respondents in the poll. This is used to
  compute the effective \eqn{N}, weighted \eqn{\bar{Y}},
  and Design Effect. Default is to make all weights equal.

  Ideally, the dimensions specified by \code{formula.cell}
  would account for all of the variation of survey weights
  in respondents in all cells. Sometimes, survey
  researchers design samples that leave some variation in
  the design weights of poststratification cells. Using
  \code{poll.weights} inflates or deflates the effective
  \eqn{N} of respondents in poststratification cells based
  on the average variance of design weights in all cells
  and each cell's deviation from that overall design
  effect.

  If multiple polls are included and contain poll.weights,
  they must be normalized within each poll before
  \emph{mrp} attempts to normalize the weights across all
  polls for all cells.}

  \item{pop.margin}{Margin of population data on which to
  sum other cells. Used in modeling PartyID, but
  implementation is not stable.}

  \item{\dots}{Additional arguments to be passed to the
  multilevel regression \code{\link[blme]{bglmer}} step.}
}
\description{
  Set up survey and population data, and a multilevel
  regression model used for poststratifying by an arbitrary
  number of strata or \dQuote{ways}.
}
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

## Ready to run simple mrp with poll and population:
mrp.simple <- mrp(ban.gaymarr ~ state+age+education+race,
                  data=CCES.complete,
                  population=mrp.census,
                  pop.weights="weighted2004")
print(100*poststratify(mrp.simple, ~ education+age), digits=2)
\dontrun{
## Fit a fuller model, adding state-level predictors:
## This model is also used in the not-run example
## for plotting maps.
mrp.statelevel <- mrp(ban.gaymarr~
                      state+f.race+age+education,
                      data=CCES.complete,
                      population=mrp.census, pop.weights="weighted2008",
                      population.formula.update= .~.-age,
                      grouplevel.data.frames=list(Statelevel,
                        mrp.regions),
                      grouplevel.expressions=list(
                        expression(age.edu <- interaction(age,education)),
                        ## an ordered factor, we use a normalized
                        ## continuous z.age as the prior mean for
                        ## varying intercepts of the 'age' groups.
                        ## That is, the prior mean for
                        ## age cat 1 of 4 (18-29) becomes (-.58)
                        expression(z.age <- rescale(age)))
                      )
## Note: the formula is expanded from the condensed version in "formula" to
##  an expanded version.
getFormula(mrp.statelevel)

## Update the model.formula on already-prepared mrp object and re-fit:
mrp.statelevel <- mr(mrp.statelevel, .~.+(1|region)+ (1|age.edu)+
                     z.age+p.relig.full+p.kerry.full)

## Fine plot control is shown with this example in plotmrp documentation!
}
}
}
\seealso{
  \code{\link{mrp-class}} for other methods on the objects
  produced by \code{mrp()}; \code{\link{plotmrp}} for how
  to plot poststratified results onto maps.
}

