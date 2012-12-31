.showmrp <- function(object) {
    message("An object of class ", sQuote(class(object)))
    dims <- dimnames(getPopulation(object))
    l <- dim(getPopulation(object))
    l <- paste0("(",l,")")
    message("Stratified on: ", serialPaste(paste(names(dims),l)))
    #message("Group-level data added: ", serialPaste(getAdded(object)))
    message("Multilevel regression formula: ", getOutcome(object), " ~")
    print(getFormula(object)[[3]], showEnv=FALSE)
}
setMethod("show", signature("mrp"), .showmrp)
