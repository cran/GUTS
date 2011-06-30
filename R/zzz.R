#modguts <- Module( "modguts" )
.onLoad <- function(libname, pkgname){
#  modguts <- Module( "modguts", PACKAGE="GUTS", mustStart=TRUE)
#  GUTS <- modguts$GUTS
  require("methods", character = TRUE, quietly = TRUE)
  loadRcppModules()
}
