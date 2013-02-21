

##' Get Number of Ways for MRP analysis
##' 
##' A cross-classified dataset for multilevel regression and poststratification
##' is an \eqn{N}-dimensional array. Each dimension in the array is one of the
##' \dQuote{ways} by which it can later be poststratified. For various reasons
##' it may be useful to query a \code{\link{mrp-class}} object or an
##' \code{\link{NWayData-class}} object for these dimensions and their names.
##' 
##' 
##' @name getNumberWays-methods
##' @aliases getNumberWays-methods getNumberWays,mrp-method
##' getNumberWays,NWayData-method getNumberWays
##' @docType methods
##' @section Methods: \describe{ \item{list("signature(object = \"mrp\")")}{
##' Returns a vector of length 2 with names \sQuote{poll} and \sQuote{pop}. In
##' the special case where no population data has been assigned (and an array
##' of 1s is used to return aggregated fitted values from the multilevel model)
##' the value returned for \sQuote{pop} is 0.  }
##' 
##' \item{list("signature(object = \"NWayData\")")}{ Returns the number of ways
##' and an attribute \sQuote{ways} containing the character vector of names of
##' the \dQuote{ways} by which it is constructed.  } }
##' @keywords methods
NULL





##' Run Multilevel Regression step of MRP Analysis
##' 
##' Run a (binomial) multilevel regression in survey data for later
##' poststratification.
##' 
##' 
##' @aliases mr mr,mrp-method
##' @param object A \code{mrp} object.
##' @param mr.formula A formula specification for the multilevel model to run
##' in the prepared data. The left-hand side should always be
##' \sQuote{\code{response}}. For convenience, the formula is handled by
##' \code{update.formula} so that \code{.} indicates the current formula
##' contents on either side of the \code{~}, e.g., \code{.~.+newVar}. The
##' initial default formula is constructed as just an intercept term for each
##' of the variables in the main formula specification
##' (\code{(1|way1)+(1|way2)} etc.)
##' @param \dots Additional arguments to be passed to the multilevel regression
##' step, which uses \code{\link[lme4]{glmer}} by default.
##' @seealso \code{\link{mrp-class}} for an example.  \code{\link{mrp-class}}
##' for other methods on the objects produced by \code{mrp()};
##' \code{\link{plotmrp}} for how to plot poststratified results onto maps.
NULL





##' Multilevel regression and poststratification
##' 
##' The \code{mrp} class holds N-dimensional cross-classified arrays for a
##' survey and a population, a 2-dimensional summary of survey data used to fit
##' a multilevel model, and a fitted model. Methods described here operate on
##' \code{mrp} objects which are typically created by calling the
##' \code{\link{mrp}} function.
##' 
##' 
##' @name mrp-class
##' @aliases mrp-class getData,mrp-method getResponse,mrp-method
##' getThetaHat,mrp-method getPopulation,mrp-method getPopulation
##' setPopulation,mrp-method setPopulationOnes,mrp-method getFormula,mrp-method
##' setFormula,mrp-method getFormula getResponse getThetaHat getEstimates
##' getEstimates,mrp-method setFormula setPopulationOnes setPopulation
##' @docType class
##' @section Slots: \describe{ \item{list("data")}{ the \code{data.frame} used
##' in the multilevel regression step. This is created by summarizing the
##' NWayData array in poll, and optionally joining additional columns of other
##' predictors. \bold{Note:} the row order of the data object, which will be
##' the row order of the fitted values, \bold{must be preserved} because it
##' corresponds to the NWayData arrays used in poststratification. To add new
##' columns onto the data.frame in the \code{mrp@data} slot, take care to
##' preserve this ordering. Base \code{merge()} will almost invariably permute
##' it in unpredictable ways: \bold{use \code{\link[plyr]{join}}} from plyr
##' instead.  This is done automatically when data.frames are joined by the
##' \code{add} argument to \code{\link{mrp}}.}\item{:}{ the \code{data.frame}
##' used in the multilevel regression step. This is created by summarizing the
##' NWayData array in poll, and optionally joining additional columns of other
##' predictors. \bold{Note:} the row order of the data object, which will be
##' the row order of the fitted values, \bold{must be preserved} because it
##' corresponds to the NWayData arrays used in poststratification. To add new
##' columns onto the data.frame in the \code{mrp@data} slot, take care to
##' preserve this ordering. Base \code{merge()} will almost invariably permute
##' it in unpredictable ways: \bold{use \code{\link[plyr]{join}}} from plyr
##' instead.  This is done automatically when data.frames are joined by the
##' \code{add} argument to \code{\link{mrp}}.} \item{list("formula")}{ The
##' formula used in the multilevel regression. The left-hand side is always
##' \sQuote{response}.}\item{:}{ The formula used in the multilevel regression.
##' The left-hand side is always \sQuote{response}.}
##' \item{list("multilevelModel")}{ The multilevel regression model (class
##' \code{mer} created by \code{\link[lme4]{glmer}})}\item{:}{ The multilevel
##' regression model (class \code{mer} created by \code{\link[lme4]{glmer}})}
##' \item{list("outcome")}{character. The name of the outcome
##' variable.}\item{:}{character. The name of the outcome variable.}
##' \item{list("poll")}{ An N-way array (\code{\link{NWayData-class}})
##' constructed from a survey.}\item{:}{ An N-way array
##' (\code{\link{NWayData-class}}) constructed from a survey.}
##' \item{list("population")}{ The population distribution, also of
##' \code{\link{NWayData-class}}, with dimensions matching those of
##' \code{poll}. This is intended to be a probability mass table (all entries
##' sum to 1), but the package will also work with unnormalized population
##' distributions.  The population is not used for the multilevel regression
##' step.  It is used only in poststratification.}\item{:}{ The population
##' distribution, also of \code{\link{NWayData-class}}, with dimensions
##' matching those of \code{poll}. This is intended to be a probability mass
##' table (all entries sum to 1), but the package will also work with
##' unnormalized population distributions.  The population is not used for the
##' multilevel regression step.  It is used only in poststratification.} }
##' @seealso \code{\link{mrp}}, which produces \code{mrp} objects.
##' @keywords classes
NULL





##' N-Dimensional Arrays for Multilevel Regression and Poststratification
##' 
##' Arrays used in multilevel regression and poststratification are
##' low-dimensional (usually 2 or 3) slices of survey respondents or population
##' frequencies. For the multilevel regression step, dimensions of such arrays
##' are then summarized using methods here taking into account survey weights
##' and associated design effects.
##' 
##' 
##' @name NWayData-class
##' @aliases NWayData-class getDesignEffect,NWayData-method
##' getN,NWayData-method getNEffective,NWayData-method
##' getYbarWeighted,NWayData-method is.NWayData getDesignEffect getN
##' getNEffective getYbarWeighted
##' @docType class
##' @section Slots: \describe{ \item{.Data}{An array.} \item{type}{Character,
##' either \sQuote{poll} or \sQuote{pop} indicating the origin of the data
##' contained in the array.} \item{levels}{When the array is collapsed back to
##' a two-dimensional form, array dimension labels become character and any
##' ordering of the original factor levels is lost. To restore level names and
##' order, those attributes (in a \code{list}) are stored here.} }
##' @author Andrew Gelman <gelman@@stat.columbia.edu>, Daniel Lee
##' <bearlee@@alum.mit.edu>, Yu-Sung Su <ys463@@columbia.edu>, Michael Malecki
##' <malecki@@wustl.edu>
##' @seealso \code{\link{mrp-class}}
##' @keywords classes
NULL





##' Poststratification method
##' 
##' Poststratify multilevel regression model by an arbitrary number of strata
##' or \dQuote{ways}. By default this method returns a single poststratified
##' predicted value.
##' 
##' 
##' @aliases poststratify poststratify,mrp-method poststratify,NWayData-method
##' poststratify,jagsNWayData-method
##' @param object An \code{mrp}, \code{NWayData}, or \code{jagsNWayData}
##' object.
##' @param formula A formula representation of the desired poststratification.
##' The formula is \code{NULL} on the left-hand side and right-hand side
##' variable names corresponding to the \dQuote{ways} in the population data by
##' which to poststratify.  The right-hand side can also be a character vector
##' of such names or a logical vector of length \dQuote{ways}.
##' 
##' See example in \code{\link{mrp}}.
##' @param fun The function (default=\emph{mean}) to summarize the collapsed
##' dimensions.
##' @param population An array or \code{NWayData} with dimensions matching
##' \code{object}, used to produce population-weighted estimates from
##' \code{jagsNWayData.}
##' @seealso \code{\link{mrp-class}} for an example.  \code{\link{mrp-class}}
##' for other methods on the objects produced by \code{mrp()};
##' \code{\link{plotmrp}} for how to plot poststratified results onto maps.
NULL





##' Plot poststratified estimates onto sp map objects
##' 
##' Poststratify and plot a fitted \code{\link{mrp-class}} object onto a
##' \link[sp]{SpatialPolygonsDataFrame}. Note that package \emph{mrp} masks the
##' \emph{sp} function \code{\link[sp]{panel.polygonsplot}}: in particular it
##' plots \code{NA} values in a color, and enables a custom stroke on selected
##' polygons in the map.
##' 
##' 
##' @name plotmrp
##' @aliases spplot,mrp-method panel.polygonsplot plotmrp
##' @docType methods
##' @param obj A \code{\link{mrp-class}} object.
##' @param formula The desired poststratification, where the left-hand side
##' indicates the name of the variable in the \code{mrp} data that identifies
##' the polygons to plot. In the special case where the desired
##' poststratification is \emph{only} by geographic unit, place it on both
##' sides of the equation (e.g., \code{state ~ state}).
##' @param spmap A map of class \code{\link[sp]{SpatialPolygonsDataFrame}},
##' such as is produced by \code{\link[maptools]{readShapeSpatial}} for reading
##' ESRI \dQuote{shapefiles.} This package includes a US map called
##' \link{spmap.states}, which is the default here.
##' @param FID The name of the column in the \code{spmap} to use as unique
##' \dQuote{feature IDs}, whose values to the left-hand-side of \code{formula}.
##' In the included \link{spmap.states} map, the column with two-letter state
##' abbreviations is called \sQuote{STATE}, which is the default here.
##' @param exclude A vector of feature IDs to exclude from the plot.
##' \emph{Note:} any rows not present in the \code{mrp} data itself will
##' automatically be excluded.
##' @param stroke A list of feature IDs (\bold{character}, \bold{logical}, or
##' \bold{integer} with respect to the order of \code{spmap}) to plot with a
##' distinct stroke; or an \bold{expression} to evaluate in the
##' \code{mrp@data}. Expressions allow data pertaining to geographic units,
##' joined via the \code{add} argument to \code{\link{mrp}} to be used in
##' plotting.
##' 
##' The type, width, and color of the stroke(s) are taken from the
##' \emph{superpose.line} list, which can be customized as described below. In
##' the event that a feature ID appears multiple times in the list, the last
##' value is used; e.g., if it is in both the first and second vectors of the
##' \sQuote{stroke} list, its \sQuote{type} will be
##' \code{superpose.line$lty[2]}.
##' @param add.settings A list of lattice graphical parameters, suitable for
##' \code{trellis.par.set()}. These will override the MRP defaults, which bear
##' mention for two reasons: 1) top and left strip titles are drawn only in
##' panels on the top and left of the full lattice. 2) any \code{NA}-valued
##' parts of any panel will be filled in the color given by
##' \sQuote{reference.line$col}, which in our theme defaults to a 68\% opaque
##' black, or \code{"#00000044"}.
##' @param subset Expression that evalues to a logical or integer indexing
##' vector, evaluated in the \code{\link[reshape2]{melt}}ed poststratified
##' results. N.B. because of partial argument matching, to place a subheading
##' blow the plot using the argument \sQuote{sub} you have to include
##' \code{subset=TRUE}.
##' @param at See \link[lattice]{levelplot}. If a vector of values is provided
##' to \code{at}, it overrides \code{cuts}, \code{pretty}, and \code{center}.
##' @param cuts See \link[lattice]{levelplot}. If using a
##' \code{\link[RColorBrewer]{RColorBrewer}} palette, should be one less than
##' the number of colors in the call to \code{\link[RColorBrewer]{brewer.pal}}.
##' @param pretty See \link[lattice]{levelplot}
##' @param center Value at which to center the calculation of limits.
##' Default=0.5
##' @param between Documented for \code{\link{xyplot}}. A named list giving
##' space between panels (in character heights). Default is 0.25 for both
##' \eqn{x} and \eqn{y}.
##' @param \dots further arguments to be passed to \code{\link[sp]{spplot}}
##' which in turn passes them to \code{\link[lattice]{levelplot}}.
##' @return An object of class \code{"trellis"}. See
##' \code{\link[lattice]{xyplot}} for details.
##' @seealso \code{\link[sp]{spplot}} for the parent function of the
##' \code{spplot,mrp-method}; \code{spplot} in turn extends
##' \code{\link[lattice]{levelplot}}.
##' 
##' Several packages produce gradients of colors which can be supplied to
##' \code{spplot} as \code{regions=list(col=\dots{})} in the
##' \sQuote{add.settings} list. Among them are
##' \link[RColorBrewer]{RColorBrewer}, \code{\link[colorRamps]{colorRamps}},
##' and \code{colorspace}, along with the basic grDevices palettes such as
##' \code{\link[grDevices]{heat.colors}} and the function
##' \code{\link[grDevices]{colorRampPalette}}.
##' 
##' A wrapper for all these nice palettes (including spline interpolation for
##' the Brewer palettes) is in package \emph{fBasics}; it is used in the
##' example below
##' @keywords methods
##' @examples
##' 
##' \donttest{
##' ## plot the example from mrp()
##' example(mrp)
##' ## Dimension of plot are by order of terms -- y,x so
##' ## this will plot age along x and poll along y.
##' spplot(mrp.simple, formula=state~ age,
##'        spmap=spmap.states, FID="STATE", 
##' 	   exclude=c("AK","DC","HI"),
##'        stroke=list(c("IA","NH","VT","CT","MA","ME"), "CA"),
##'        colorkey=list(
##'          space="bottom",height=.5,width=.5)
##'        )
##' \dontrun{
##' spplot(mrp.statelevel, state ~ age+edu, cuts=50,
##'        subset=TRUE,
##'        spmap.states, "STATE", exclude=c("AK","DC","HI"),
##'        stroke=list(expression(hasmarriage2010==TRUE),
##'          "CA"),
##'        sub=paste("National average:",
##'          format(poststratify(mrp.statelevel),digits=2)),
##'        add.settings=list(
##'          regions=list(col=fBasics:::divPalette(51,"BrBG"))
##'          ),
##'        colorkey=list(
##'          space="bottom",height=.5,width=.5,
##'          labels=list(at=c(.1,.5,.9),labels=c("10%","50%","90%"), cex=.7)
##'          )
##'        )
##' 
##' 
##' }
##' }
##' 
NULL





##' Run k-fold cross validations on mrp model
##' 
##' Run a (binomial) multilevel regression in survey data for later
##' poststratification.
##' 
##' 
##' @aliases xval xval,mrp-method
##' @param object A \code{mrp} object.
##' @param formula A formula specification for the multilevel model on which
##' the k-folds cross validation is run. The default is the formula of the
##' \code{mrp} object. The left-hand side should always be
##' \sQuote{\code{response}}. For convenience, the formula is handled by
##' \code{update.formula} so that \code{.} indicates the current formula
##' contents on either side of the \code{~}, e.g., \code{.~.+newVar}. The
##' initial default formula is constructed as just an intercept term for each
##' of the variables in the main formula specification
##' (\code{(1|way1)+(1|way2)} etc.)
##' @param folds Number of folds for cross valiation. Default is 4.
##' @param loss.type Type of loss measure used in cross validation. Currently
##' "logloss" is supported.
##' @param \dots Additional arguments to be passed to the multilevel regression
##' step, which uses \code{\link[blme]{bglmer}} by default.
NULL



