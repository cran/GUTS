# GUTS R definitions.
# carlo.albert@eawag.ch, soeren.vogel@uzh.ch
# Version 0.1.45
# GPL-2

# Method show, print, showObject
.guts.print <- function( x ) {
  cat(
    "\n",
    "GUTS object with the following attributes:\n\n",

    "Vector of concentrations (C, ",
    length(x$C),
    " elements):\n\t",
    paste(x$C, sep="", collapse=", "), "\n",
    
    "Vector of time points of concentrations (Ct, ",
    length(x$Ct),
    " elements):\n\t",
    paste(x$Ct, sep="", collapse=", "), "\n",

    "Vector of survivors (y, ",
    length(x$y),
    " elements):\n\t",
    paste(x$y, sep="", collapse=", "), "\n",

    "Vector of time points of survivors (yt, ",
    length(x$yt),
    " elements):\n\t",
    paste(x$yt, sep="", collapse=", "), "\n",

    "Parameters (par, ",
    length(x$par),
    " elements):\n\t",
    paste(x$par, sep="", collapse=", "), "\n",

    "Time grid points (M) : ", x$M, "\n",
    "Distribution (dist)  : ", dQuote(x$dist), "\n",
    "Sample length (N)    : ", x$N, "\n",
    
    "Sample (z, ",
    length(x$z),
    " elements):\n\t",
    "Min=", min(x$z), ", Max=", max(x$z), ", Mean=", mean(x$z), ", Sigma=", sd(x$z), "\n",
    
    "\n\n",
    sep=""
  )
}

print.GUTS <- function( x, ... ) {
  out <- .guts.print( x )
  return( invisible( out ) )
}
setMethod("show", "Rcpp_GUTS", function( object ) .guts.print( object ) )
