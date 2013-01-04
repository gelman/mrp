pkgname <- "mrp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mrp')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("mrp")
### * mrp

flush(stderr()); flush(stdout())

### Name: mrp
### Title: Multilevel regression and poststratification
### Aliases: mrp

### ** Examples




cleanEx()
nameEx("spplot")
### * spplot

flush(stderr()); flush(stdout())

### Name: plotmrp
### Title: Plot poststratified estimates onto sp map objects
### Aliases: spplot,mrp-method panel.polygonsplot plotmrp
### Keywords: methods

### ** Examples




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
