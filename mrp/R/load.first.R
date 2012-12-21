.onAttach <- function(...) {
  mylib <- dirname(system.file(package = "mrp"))
  pkgdesc <- packageDescription("mrp") 
  packageStartupMessage(paste("\nmrp (Version ", pkgdesc$Version, ", built: ", pkgdesc$Date, ")\n", sep = ""))
}
