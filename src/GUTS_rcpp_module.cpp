/** GUTS Rcpp module.
 * @file    GUTS_rcpp_module.cpp
 * @author  carlo.albert@eawag.ch, soeren.vogel@eawag.ch
 * @date    30 June 2011
 * @license: GPL-2
 */

#ifndef guts_h
#include "GUTS.h"
#endif
#include <Rcpp.h>

using namespace Rcpp;

RCPP_MODULE(modguts)
{
    class_<GUTS>( "GUTS" )
        .constructor()
        
        .method( "setConcentrations",   &GUTS::setConcentrations, "Set time series vector of concentrations." )
        .method( "setSurvivors",        &GUTS::setSurvivors, "Set time series vector of survivors." )
        .method( "setParameters",       &GUTS::setParameters, "Set parameter vector of the object." )
        .method( "setTimeGridPoints",   &GUTS::setTimeGridPoints, "Set number of grid points on the time-axis." )
        .method( "setDistribution",     &GUTS::setDistribution, "Set distribution to sample from." )
        .method( "setSampleLength",     &GUTS::setSampleLength, "Set length of sample (numerical exactness)." )
        
        .method( "setSample", &GUTS::setSample, "Set ordered sample vector." )

        .method( "calcLoglikelihood", &GUTS::calcLoglikelihood, "Returns calculated logarithm of the likelihood from complete and valid object." )

        .property( "C",     &GUTS::getC, "Vector of concentrations." )
        .property( "Ct",    &GUTS::getCt, "Time vector of concentrations." )
        .property( "y",     &GUTS::gety, "Vector of survivors." )
        .property( "yt",    &GUTS::getyt, "Time vector of survivors." )
        .property( "par",   &GUTS::getpar, "Parameter vector." )
        .property( "M",     &GUTS::getM, "Grid points on time axis." )
        .property( "dist",  &GUTS::getdist, "Distribution." )
        .property( "N",     &GUTS::getN, "Sample length (numerical exactness)." )
        .property( "z",     &GUTS::getz, "Sample." )

//        .method( "showObject", &GUTS::showObject, "Prints object information to stdout." )
    ;
}
