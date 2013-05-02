\docType{methods}
\name{spplot}
\alias{plotmrp}
\alias{spplot}
\alias{spplot,NWayData-method}
\alias{spplot,mrp-method}
\title{Plot poststratified estimates onto sp map objects}
\arguments{
  \item{obj}{A \code{\link{mrp-class}} object.}

  \item{formula}{The desired poststratification, where the
  left-hand side indicates the name of the variable in the
  \code{mrp} data that identifies the polygons to plot. In
  the special case where the desired poststratification is
  \emph{only} by geographic unit, place it on both sides of
  the equation (e.g., \code{state ~ state}).}

  \item{spmap}{A map of class
  \code{\link[sp]{SpatialPolygonsDataFrame}}, such as is
  produced by \code{\link[maptools]{readShapeSpatial}} for
  reading ESRI \dQuote{shapefiles.} This package includes a
  US map called \link{spmap.states}, which is the default
  here.}

  \item{FID}{The name of the column in the \code{spmap} to
  use as unique \dQuote{feature IDs}, whose values to the
  left-hand-side of \code{formula}. In the included
  \link{spmap.states} map, the column with two-letter state
  abbreviations is called \sQuote{STATE}, which is the
  default here.}

  \item{exclude}{A vector of feature IDs to exclude from
  the plot. \emph{Note:} any rows not present in the
  \code{mrp} data itself will automatically be excluded.}

  \item{stroke}{A list of feature IDs (\bold{character},
  \bold{logical}, or \bold{integer} with respect to the
  order of \code{spmap}) to plot with a distinct stroke; or
  an \bold{expression} to evaluate in the \code{data}.
  Expressions allow data pertaining to geographic units,
  joined via the \code{add} argument to \code{\link{mrp}}
  to be used in plotting.

  The type, width, and color of the stroke(s) are taken
  from the \emph{superpose.line} list, which can be
  customized as described below. In the event that a
  feature ID appears multiple times in the list, the last
  value is used; e.g., if it is in both the first and
  second vectors of the \sQuote{stroke} list, its
  \sQuote{type} will be \code{superpose.line$lty[2]}.}

  \item{add.settings}{A list of lattice graphical
  parameters, suitable for \code{trellis.par.set()}. These
  will override the MRP defaults, which bear mention for
  two reasons: 1) top and left strip titles are drawn only
  in panels on the top and left of the full lattice. 2) any
  \code{NA}-valued parts of any panel will be filled in the
  color given by \sQuote{reference.line$col}, which in our
  theme defaults to a 68\% opaque black, or
  \code{"#00000044"}.}

  \item{subset}{Expression that evalues to a logical or
  integer indexing vector, evaluated in the
  \code{\link[reshape2]{melt}}ed poststratified results.
  N.B. because of partial argument matching, to place a
  subheading blow the plot using the argument \sQuote{sub}
  you have to include \code{subset=TRUE}.}

  \item{at}{See \link[lattice]{levelplot}. If a vector of
  values is provided to \code{at}, it overrides
  \code{cuts}, \code{pretty}, and \code{center}.}

  \item{cuts}{See \link[lattice]{levelplot}. If using a
  \code{\link[RColorBrewer]{RColorBrewer}} palette, should
  be one less than the number of colors in the call to
  \code{\link[RColorBrewer]{brewer.pal}}.}

  \item{pretty}{See \link[lattice]{levelplot}}

  \item{center}{Value at which to center the calculation of
  limits. Default=0.5}

  \item{between}{Documented for \code{\link{xyplot}}. A
  named list giving space between panels (in character
  heights). Default is 0.25 for both \eqn{x} and \eqn{y}.}

  \item{\dots}{further arguments to be passed to
  \code{\link[sp]{spplot}} which in turn passes them to
  \code{\link[lattice]{levelplot}}.}
}
\value{
  An object of class \code{"trellis"}. See
  \code{\link[lattice]{xyplot}} for details.
}
\description{
  Poststratify and plot a fitted \code{\link{mrp-class}}
  object onto a \link[sp]{SpatialPolygonsDataFrame}. Note
  that package \emph{mrp} masks the \emph{sp} function
  \code{\link[sp]{panel.polygonsplot}}: in particular it
  plots \code{NA} values in a color, and enables a custom
  stroke on selected polygons in the map.
}
\examples{
\donttest{
## plot the example from mrp()
example(mrp)
## Dimension of plot are by order of terms -- y,x so
## this will plot age along x and poll along y.
spplot(mrp.simple, formula=state~ age,
       spmap=spmap.states, FID="STATE",
	   exclude=c("AK","DC","HI"),
       stroke=list(c("IA","NH","VT","CT","MA","ME"), "CA"),
       colorkey=list(
         space="bottom",height=.5,width=.5)
       )
\dontrun{
spplot(mrp.statelevel, state ~ age+edu, cuts=50,
       subset=TRUE,
       spmap.states, "STATE", exclude=c("AK","DC","HI"),
       stroke=list(expression(hasmarriage2010==TRUE),
         "CA"),
       sub=paste("National average:",
         format(poststratify(mrp.statelevel),digits=2)),
       add.settings=list(
         regions=list(col=fBasics:::divPalette(51,"BrBG"))
         ),
       colorkey=list(
         space="bottom",height=.5,width=.5,
         labels=list(at=c(.1,.5,.9),labels=c("10\%","50\%","90\%"), cex=.7)
         )
       )


}
}
}
\seealso{
  \code{\link[sp]{spplot}} for the parent function of the
  \code{spplot,mrp-method}; \code{spplot} in turn extends
  \code{\link[lattice]{levelplot}}.

  Several packages produce gradients of colors which can be
  supplied to \code{spplot} as
  \code{regions=list(col=\dots{})} in the
  \sQuote{add.settings} list. Among them are
  \link[RColorBrewer]{RColorBrewer},
  \code{\link[colorRamps]{colorRamps}}, and
  \code{colorspace}, along with the basic grDevices
  palettes such as \code{\link[grDevices]{heat.colors}} and
  the function \code{\link[grDevices]{colorRampPalette}}.

  A wrapper for all these nice palettes (including spline
  interpolation for the Brewer palettes) is in package
  \emph{fBasics}; it is used in the example below
}
\keyword{methods}

