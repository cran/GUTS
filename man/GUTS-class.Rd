\encoding{UTF-8}
\name{GUTS-class}
\Rdversion{1.1}
\docType{class}

\alias{GUTS-class}
\alias{GUTS-method}
\alias{GUTS-package}
\alias{GUTS}
\alias{Rcpp_GUTS-class}
\alias{Rcpp_GUTS}
\alias{modguts}
\alias{print.guts}

\title{Class \code{GUTS} (\code{"Rcpp_GUTS"})}

\description{\code{GUTS} represents survivor model experiments for which the logarithm of the likelihood can be calculated. \code{GUTS} is the constructor object of the C++ class “GUTS” exposed to R through the \link{Rcpp} module \code{modguts}. The reference class for \code{GUTS} is \code{"Rcpp_GUTS"}.}

\section{Objects from the Class}{
New objects can be created by calls to the \bold{reference method} \code{new()} (see \link{ReferenceClasses}), or the \bold{function} \code{new("Rcpp_GUTS")}. Currently, no arguments must be provided, and a “factory fresh” \code{GUTS} object is filled with default non-sense values. A \code{GUTS} object has no public fields. Setter methods must be used to set values of an object.
}

\section{Slots}{The class has no slots. Private fields can be set using methods (see below).}

\section{Methods}{
  \describe{
    \item{new}{\code{signature()}: Create a new \code{GUTS} object with non-sense values.}

    \item{setConcentrations}{\code{signature(C = "double", Ct = "double")}:
      Set the vectors of concentrations (\code{C}) and concentration time points (\code{Ct}).
      \code{C} and \code{Ct} must have the same length.
      \code{Ct} and \code{yt} must use the same time units.}

    \item{setSurvivors}{\code{signature(y = "integer", yt = "double")}:
      Set the vectors of survivors (\code{y}) and survivor time points (\code{yt}).
      \code{y} must be a vector of integers (counts).
      \code{y} and \code{yt} must have the same length.
      \code{Ct} and \code{yt} must use the same time units.
    }

    \item{setParameters}{\code{signature(par = "double")}:
      Set parameters (\code{length(par) > 2}) for the experiment.
      The first three parameters are used for the calculation of the logarithm of the likelihood.
      More parameters may be appended and are used to construct samples of length \code{N} (see \code{setSampleLength()})
      from distribution \code{dist} (see \code{setDistribution()}).
    }

    \item{setTimeGridPoints}{\code{signature(M = "integer")}:
      Set the number of grid points on the time axis (numerical exactness).
    }
    
    \item{setDistribution}{\code{signature(dist = "character")}:
      Set the \bold{name} of the distribution to create.
      See section \bold{Distributions} for details.
      Do not append distribution parameters here!
      Ignored if \code{setSample} is used.
    }

    \item{setSampleLength}{\code{signature(N = "integer")}:
      Set the length of the sample that is drawn from the distribution (numerical exactness).
      Ignored if \code{setSample} is used.
    }

    \item{setSample}{\code{signature(z = "double")}:
      Set an empirical sample \bold{instead} of a theoretical distribution to sample from.
      \code{setSample} will override \code{setDistribution} and \code{setSampleLength}, and
        sets \code{dist} to \"empirical\", and \code{N} to the length of \code{z}.
    }

    \item{calcLoglikelihood}{\code{signature()}:
      Calculate the logarithm of the likelihood of a properly set-up \code{"GUTS"} object.
    }

    \item{showObject}{\code{signature()}:
      Display current values of a GUTS object.
    }
    \item{show}{\code{signature()}:
      S4 method wrapper for method \code{showObject}.
    }
  }
}

\section{Distributions}{%
\code{GUTS} creates samples internally according to \code{dist}, \code{N}, and \code{pars[4:length(pars)]}. If no parameters for the creation of the distribution are given, default values are applied. \code{dist} must be one of the following alternatives:
  \describe{
    \item{\code{"normal"}}{Create a sample of size \code{N} from a normal distribution
      with \code{mean = pars[4]}, and \code{sigma = pars[5]}.}
    \item{\code{"lognormal"}}{Create a sample of size \code{N} from a lognormal distribution
      with \code{mean = pars[4]}, and \code{sigma = pars[5]}. Note that mean and sigma are parameters for the
      \bold{lognormal} distribution, not its logarithm (differs from \code{rlnorm} in R)!}
  }
}

\section{Functions}{
  \describe{
    \item{print.guts(object)}{Simple wrapper for \code{showObject}.}
  }
}

\section{Notes and warnings}{
  \describe{
    \item{Compiling from source}{Compilation from source would require an installation of the C++ boost library (\href{http://www.boost.org/}{http://www.boost.org/}). The library is expected in \code{/opt/local/include} under UNIX-likes, and in \code{/Boost} under Windows. Change the makefiles in \code{GUTS/src} to fit your installation.}
    \item{Serialisation}{Serialisation is currently not implemented. Objects of class \code{GUTS} cannot be saved across sessions. Clean out every \code{GUTS} objet from your workspace before leaving R. On restart, rebuild \code{GUTS} objects if necessary.}
    \item{Development status}{The package is under current development. \bold{Save your work before working with it.}}
  }
}

\references{Jager T., Albert C., Preuss T., Ashauer R., General unified theory of survival – a toxicokinetic toxicodynamic framework for ecotoxicology, published online in Env. Sci. and Tech. (2011) (\href{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/\#}{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/#}).}

\author{%
Carlo Albert & Sören Vogel \email{guts.package@gmail.com}

Maintainer: Carlo Albert & Sören Vogel \email{guts.package@gmail.com}%
}

\seealso{\code{\link{ReferenceClasses}}, \code{\link{Rcpp}}, \code{\link{new}}}

\examples{
library("GUTS")

getRefClass("Rcpp_GUTS")
GUTS # is a C++ class
getRefClass(GUTS) # is a generator object

obj0 <- new("Rcpp_GUTS")
obj0$calcLoglikelihood() # 0

obj1 <- GUTS$new()
obj1$calcLoglikelihood() # 0

obj1$setConcentrations(
   C = c(102.65, 97.59, 0.00, 0.00, 103.88, 98.19, 0.00, 0.00,
         0.00, 0.00),
  Ct = c(0.00, 1.02, 1.03, 2.99, 3.01, 4.01, 4.02, 11.01,
         18.01, 22.01)
)
obj1$setSurvivors(
   y = c(70, 66, 61, 55, 31, 31, 29, 26, 24, 22, 21, 19, 17,
         14, 14, 13, 11, 11, 10, 9, 8, 8, 8),
  yt = 0:22
)
obj1$setParameters(
 par = c(0.05084761, 0.12641525, 1.61840054, 19.09911,
         6.495246)
)
obj1$setTimeGridPoints( M = 10000 )
obj1$setDistribution( dist = "lognormal" )
obj1$setSampleLength( N = 10000 )
obj1$calcLoglikelihood() # abt -183

obj2 <- obj1
obj2$calcLoglikelihood() # abt -183

z <- sort( rlnorm(10000,
  2.8943280, 0.33420543) ) # *meanlog* and *meansd*!!!
obj2$setSample( z = z )
obj2$calcLoglikelihood() # abt -183
#\dontrun{
#}
}

\keyword{classes}
\keyword{GUTS-class}
