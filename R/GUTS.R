# GUTS R definitions.
# carlo.albert@eawag.ch, soeren.vogel@uzh.ch
# 30 Jun 2011
# GPL-2

# Method show, print, showObject
.guts.print <- function(x) {
  x$showObject()
}
print.guts <- function(x, ...) {
  out <- .guts.print(x)
  return(invisible(out))
}
setMethod("show", "Rcpp_GUTS", function(object) .guts.print(object) )
