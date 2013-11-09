

"mrpplot" <- function(obj, formula, spmap=spmap.states,
                      FID="STATE", exclude=NULL, stroke=NULL, subset=TRUE,
                      at, cuts=63, pretty=FALSE, center=0.5,
                      between=list(x=.25,y=.25),
                      add.settings=list(), ...) {
    plot.terms <- terms(formula, keep.order=TRUE)
    obj.p <- melt(poststratify(obj, all.vars(formula)))

    if(length(all.vars(plot.terms))==1){
        obj.p <- data.frame(fid=rownames(obj.p), value=obj.p[,1])
        names(obj.p)[1] <- all.vars(plot.terms)
    }

    if(is(obj,"mrp")) {
        obj.p <- restoreNWayLevels(obj.p, getPopulation(obj))
        post.levels <- getPopulation(obj)@levels
    } else if (is(obj,"NWayData")) {
        obj.p <- restoreNWayLevels(obj.p, obj)
        post.levels <- obj@levels
    }
    subset <- eval(substitute(subset), obj.p)
    obj.p <- obj.p[subset,]
    if(length(all.vars(plot.terms))>1) {
        plotdf <- dcast(obj.p, formula)
    } else { plotdf <- obj.p }

    names(plotdf) <- make.names(names(plotdf))

    if(is.factor(plotdf[,1])) {
        plotdf[,1] <- as.character(levels(plotdf[,1])[plotdf[,1]])
    }
    if(is.factor(spmap@data[,FID])) {
        spmap@data[,FID] <- as.character(levels(spmap@data[,FID])[spmap@data[,FID]])
    }

    ## remove excludes list from fitted model, and then pare the sp obj to match.
    if(!is.null(exclude)) {
        plotdf <- subset(plotdf,
                         !(plotdf[,1] %in% exclude))
    }
    spmap <- spmap[spmap@data[,FID] %in% plotdf[,1],]
    ## set feature ids and rownames to match
    spmap <- spChFIDs(spmap,
                      spmap@data[,FID])
    rownames(spmap@data) <- spmap@data[,FID]
    rownames(plotdf) <- plotdf[,1]

    startcol <- ncol(spmap@data)+2
    spmap@data <- cbind(spmap@data,
                        plotdf[rownames(spmap@data),])
    endcol <- ncol(spmap@data)

    ## Process the stroke list
    if(!is.null(stroke)) {
        stroke <- doStrokeList(stroke, spmap@data, obj@data,
                               ## lc fid is column name in data
                               fid=names(plotdf)[1])
    } else {
        stroke <- rep(NA, nrow(spmap))
    }
    dimlist <- post.levels[attr(plot.terms,"term.labels")]
    dimlabels <- lapply(names(dimlist),
                        function(x) { labels <- dimlist[[x]]
                                      return(labels[labels %in% levels(obj.p[,x])])
                                  })
    if(length(dimlabels)==1) {
        dimlabels[[2]] <- ""
        dimlabels <- rev(dimlabels)
    }
    ## Quietly swap the order.
    ## will need to be better if we allow
    ## multipage plots.
    if(length(all.vars(plot.terms))==1){
        dimlabels=list(c(""),c(""))
    }

    ## set up 'at', centered at center
    centered.range <- lattice:::extend.limits(center + c(-1,1)*{max(abs(range(as.matrix(spmap@data[,startcol:endcol])-center,
                                                                              finite=TRUE)))})
    if(diff(centered.range)>1){ centered.range <- c(0,1) }
    if (missing(at))
        at <-
            if (pretty) pretty(centered.range, cuts)
            else seq(centered.range[1], centered.range[2], length.out = cuts + 2)

    theplot <- spplot(spmap,
                      startcol:endcol,
                      layout=c(length(dimlabels[[2]]), length(dimlabels[[1]])),
                      panel=panel.polygonsplot,
                      stroke=stroke,
                      at=at,
                      strip=strip.custom(
                      factor.levels=rep(dimlabels[[2]],
                      length(dimlabels[[1]]))),
                      strip.left=strip.custom(horizontal=FALSE,
                      factor.levels=rep(dimlabels[[1]],
                      each=length(dimlabels[[2]]))),
                      between=between,
                      par.settings=lattice:::updateList(mrp.theme(
                      length(dimlabels[[1]]), # row length
                      length(dimlabels[[2]])),
                      add.settings)# col length
                      ,subset=TRUE,...
                      )

    return(theplot)
}

"doStrokeList" <- function(stroke, plotdf, datadf, fid) {
    if (!is.list(stroke)) { stroke <- list(stroke) }
    stroke <- lapply(stroke, function(scol, type=class(scol)[1]) {
        switch(type,
               "logical" = {
                   if(length(scol) != nrow(plotdf)){
                       warning("Length of logical stroke column not equal to length of polygons.")}
                   return(scol)
               },
               "integer" = {
                   ans <- rep(FALSE, nrow(plotdf))
                   ans[scol] <- TRUE
                   return(ans)
               },
               "character" = {
                   ans <- rownames(plotdf) %in% scol
                   return(ans)
               }
               ,"expression" = {
                   stroke <- eval(scol, envir=datadf)
                   ans <- data.frame(fid=rownames(plotdf))
                   scol <- tapply(stroke, datadf[,fid], any)
                   scol <- data.frame(fid=names(scol), scol)
                   ans <- join(ans, scol, by="fid")
                   return(as.logical(ans$scol))

               }
               )
    })

    stroke <- matrix(unlist(stroke),ncol=length(stroke))
    stroke <- apply(stroke,1, function(x) {
        ifelse(any(x),
               max(which(x)),
               NA )
    })
    return(stroke)
}

##' Generate a lattice theme with top and left strip titles
##'
##'
##' @param rowlength numeric, number of rows in a matrix of maps
##' @param collength numeric, number of columns in a matrix of maps
##' @return a list suitable for \code{par.settings} or \code{trellis.par.set}
##' @export
"mrp.theme" <- function(rowlength,collength){
    list(
    strip.background = list(col = "transparent"),
    reference.line = list(col="#00000044"),
    superpose.line=list(col=c("#00000066", trellis.par.get()$superpose.line$col)),
    ##par.main.text=list(fontfamily="gotham"),
    add.line=list(col="#00000022",lwd=0), # state borders
    add.text=list(cex=.7,fontface="italic"),
    ##  fontfamily="gotham"),
    axis.line=list(lwd=0),
    ## Here we are going to do some
    ## strip and strip.left magic.
    layout.heights=list(strip =
    c(rep(0, rowlength-1),
      1)), # should by dynamic for linebreaks
    layout.widths=list(strip.left=
    c(1, rep(0,collength-1))),
    strip.border=list(col="transparent"),
    regions=list(col=timPalette(),alpha=1)
    )}


"panel.polygonsplot" <- function (x, y, z, subscripts, at = pretty(z), shrink, labels = NULL,
                                  label.style = c("mixed", "flat", "align"), contour = FALSE,
                                  region = TRUE, col = add.line$col, lty = add.line$lty, lwd = add.line$lwd,
                                  cex = add.text$cex, font = add.text$font, fontfamily = add.text$fontfamily,
                                  fontface = add.text$fontface, col.text = add.text$col, ...,
                                  col.regions = regions$col, alpha.regions = regions$alpha,
                                  grid.polygons, sp.layout,stroke=stroke)
{
    regions <- trellis.par.get("regions")
    add.line <- trellis.par.get("add.line")
    add.text <- trellis.par.get("add.text")
    numcol <- length(at) - 1
    numcol.r <- length(col.regions)
    col.regions <- if (numcol.r <= numcol)
        rep(col.regions, length = numcol)
    else col.regions[floor(1 + (1:numcol - 1) * (numcol.r - 1)/(numcol - 1))]
    zcol <- rep(NA, length(z))
    for (i in seq(along = col.regions)) zcol[!is.na(x) & !is.na(y) &
                  !is.na(z) & z >= at[i] & z < at[i + 1]] <- i
    label.style <- match.arg(label.style)
    x <- as.numeric(x[subscripts])
    y <- as.numeric(y[subscripts])
    z <- as.numeric(z[subscripts])
    zcol <- as.numeric(zcol[subscripts])

                                        #EJP,2010-10-8:
                                        #if (is(grid.polygons, "SpatialLines"))
                                        #	sp.panel.layout(sp.layout, panel.number())
    sp:::sp.panel.layout(sp.layout, panel.number(), first = TRUE)
    if (any(subscripts)) {
        if (is(grid.polygons, "SpatialLines")) {
            sp.lines3 = function(x, col, ...) panel.lines(coordinates(x), col = col, ...)
            sp.lines2 = function(x, col, ...) lapply(x@Lines, sp.lines3, col, ...)
            for (i in 1:length(grid.polygons@lines))
                sp.lines2(grid.polygons@lines[[i]], col = col.regions[zcol[i]], lwd = lwd, lty = lty, ...)
        } else {
            pls = slot(grid.polygons, "polygons")
            pO = slot(grid.polygons, "plotOrder")
            for (i in pO) {

                Srs <- slot(pls[[i]], "Polygons")
                pOi <- slot(pls[[i]], "plotOrder")
                for (j in pOi) {
                    coords = slot(Srs[[j]], "coords")
                    if (slot(Srs[[j]], "hole")) {
                        bg = trellis.par.get()$background
                        if (bg$col == "transparent")
                            fill = "white"
                        else
                            fill = bg$col
                        alpha = bg$alpha
                    } else {
                        if(is.na(zcol[i])) {
                            fill <- trellis.par.get()$reference.line$col
                        } else {
                            fill = col.regions[zcol[i]]
                            alpha = alpha.regions
                        }
                        gp = grid:::gpar(fill = fill, alpha = alpha, col = col, lwd = lwd, lty = lty)
                        if(!is.na(stroke[i])) {
                            mystroke <- trellis.par.get()$superpose.line
                            gp <- grid:::gpar(fill = fill, alpha = alpha,
                                              col = mystroke$col[stroke[i]],
                                              lwd = mystroke$lwd[stroke[i]],
                                              lty = mystroke$lty[stroke[i]])
                        }
                    }
                    grid.polygon(coords[,1], coords[,2], default.units = "native",
                                 gp = gp)
                }
            }
        }
    }
                                        # EJP, 2010-10-8
                                        #if (!is(grid.polygons, "SpatialLines"))
                                        #	sp.panel.layout(sp.layout, panel.number())
    sp:::sp.panel.layout(sp.layout, panel.number())
}


timPalette <-
    function(n = 64)
{
### FIX THIS TO USE colorspace PACKAGE INSTEAD FOR HCL


    ## A function implemented by Diethelm Wuertz
    ## taken from package fBasics. This is the default palette for MRP plots

    ## Description:
    ##   Creates a cyan, yellow, to orange palette

    ## Notes:
    ##   'Tim.colors' in 'fields' package goes from blue to red, and passes
    ##   through the colors cyan, yellow, and orange. Also known as Jet
    ##   color-map in Matlab. You can also easily design your own color map
    ##   using 'rgb' function from 'gdDevices'.
    ##   From:  <Jaroslaw.W.Tuszynski@saic.com>
    orig = c(
    "#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
    "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
    "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
    "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
    "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
    "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
    "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
    "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
    "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
    "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
    "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
    "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
    "#AF0000", "#9F0000", "#8F0000", "#800000")
    if (n == 64) return(orig)
    rgb.tim = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , 64)
    xg = seq(0, 1, , n)
    for (k in 1:3) {
        hold = spline(x, rgb.tim[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    ans = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    ans
}


##' Plot poststratified estimates onto sp map objects
##'
##' Poststratify and plot a fitted \code{\link{mrp-class}} object onto a
##' \link[sp]{SpatialPolygonsDataFrame}. Note that package \emph{mrp} masks the
##' \emph{sp} function \code{\link[sp]{panel.polygonsplot}}: in particular it
##' plots \code{NA} values in a color, and enables a custom stroke on selected
##' polygons in the map.
##'
##'
##' @rdname spplot-methods
##' @aliases plotmrp
##' @export
##' @name spplot
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
##' \code{data}. Expressions allow data pertaining to geographic units,
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
##' spplot(mrp.statelevel, state ~ age + education, cuts=50,
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
##' @aliases spplot,mrp-method
##' @export
setMethod("spplot", signature("mrp"), mrpplot)
##' @rdname spplot-methods
##' @aliases spplot,NWayData-method
##' @export
setMethod("spplot", signature("NWayData"), mrpplot)
