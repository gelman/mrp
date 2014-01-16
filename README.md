<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{MRP Primer Reloaded}
-->

# Package, Data, Examples, and this Vignette

The running example uses open access data from the Comprehensive Congressional Election Study from 2008. You can run the 'simple' models just by doing `example(mrp)` or `example(plotmrp)`. The purpose of this document is to explain what things do in a somewhat more friendly way than the usual R function documentation. For the vignette, we also show the more complete, closer-to-real-life 'not run' models, but these are also present in the R-doc 'examples'.



```r
library(mrp)
```

```
## Loading required package: Matrix
## Loading required package: blme
## Loading required package: lme4
## Loading required package: lattice
## Loading required package: arm
## Loading required package: MASS
## 
## arm (Version 1.6-10, built: 2013-11-15)
## 
## Working directory is /Users/malecki/mrp/mrp/vignettes
## 
## Loading required package: sp
## Loading required package: grid
## Loading required package: plyr
## Loading required package: reshape2
## Loading required package: mrpdata
## 
## mrp (Version 1.0-1, built: 2013-11-09)
```

```r
data(CCES.complete)
data(mrp.census)
data(mrp.regions)
```


## Data

The survey included in the package as `CCES.complete` is the 2008 Common Content from the CCES, a large national sample of adults in the US. We also include some state-level data as `Statelevel`, attributed to Lax, Kastellec, and Phillips; and Lax's "Demographically Purged State Predictor (DPSP)" as `dpsp`, and aggregated US Census/ACS data as `mrp.census`. A convenience data frame of (factor) US state abbreviations and the "regions" that researchers often use is included as `mrp.regions`. 

### mrp.census

A table of weighted cell counts for the Cartesian product of the geographic-demographic indicators. It is of course also possible to use "sampling frame" type data and MRP will sum across all of the unused dimensions.

### CCES.complete

This is actually a single realization from a multiple-imputation run over the CCES 2008 data. 

Stephen Ansolabehere, 2011, "CCES, Common Content, 2008", <a href="http://hdl.handle.net/1902.1/14003">hdl:1902.1/14003</a> UNF:5:7eeaUMPVCcKDNxK6/kd37w== V6 [Version]

### spmap.states

A projected sp-package map object also containing a data frame with state names, two-letter state abbreviations, and FIPS codes, which are common in survey or census data. In the example, we use two-letter abbreviations.

# Fitting a basic model

## Preparing the data

Almost all data will involve recoding and carefully checking categorical variables. We rely heavily on R's "factor" data type. Factors have associated "levels" (names for categories) and may be ordered. For MRP, categorical variables need to be factors, and the survey data's levels need to be a subset of those in the population data. If survey data exist for all population cells, so much the better, but MRP can only make predictions for cells that exist in the population.

To ensure matching of strata between poll and population, both should be factors with identical names and levels.

```r
CCES.complete <- within (CCES.complete, {
  education <- factor(education,exclude=NA)
  female <- factor(sex=="Female", labels=c("Male","Female"), exclude=NA)
  race <- factor(race,exclude=NA)
  f.race <- interaction(female,race)
})
```

The survey has four levels of education, so we need to combine
the top two levels in the census data. We'll also go ahead
and trim it down to just the variables used here.

```r

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
```


## Calling `mrp()` for the first time

In the msot basic setup, MRP will fit an intercept to each of the groups provided in the `formula.cell`. It is called this because it defines the "cells" of the hypercube of geographic-demographic classifications that are partially pooled by the multilevel model. In the first example, these categories are **state, age,** and **education**; all of these exist in the census data as well. The outcome variable **ban.gaymarr** is the response to the following question in the CCES common content:

> Congress considered many important bills over the past two years. For each of the following tell us whether you support or oppose the legislation in principle: Constitutional Amendment banning Gay Marriage ("Support" coded 1; "Oppose" coded as 0).

