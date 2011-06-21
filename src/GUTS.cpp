/** GUTS class implementation.
 * @file    GUTS.cpp
 * @author  carlo.albert@eawag.ch, soeren.vogel@eawag.ch
 * @date    30 Jun 2011
 *
 * C++ class GUTS:
 * Fast Calculation of the Likelihood of a Stochastic Survival Model
 * License: GPL-2
 */

#ifndef guts_h
#include "GUTS.h"
#endif

using namespace std;

/** GUTS class constructor.
 *
 * Defaults do not represent a working example, they
 * are only present to avoid calculation errors.
 */
GUTS::GUTS()
{
    /*
     * Attributes and public fields.
     */
    mC.assign( 10, 1.0 );
    mCt.assign( 10, 1.0 );
    my.assign( 11, 0 );         // appended field with 0
    myt.assign( 10, 1.0 );
    mpar.assign( 5, 1.0 );
    mM = 10000;
    mdist = "none";
    mN = 10000;
    mz.assign( 10000, 1.0 );

    /*
     * Current state.
     */
    D.assign( mM, 0.0 );        // length of mM, all elements 0.0
    S.assign( my.size(), 0.0 ); // length of my, any elements
    ee.assign( mN, 0.0 );       // length of mN, all elements 0.0
    ff.assign( mN, 0 );         // length of mN, all elements 0.0

    /*
     * dtau is updated in setSurvivors(), AND setTimeGripPoints()
     */
    dtau = 0.1;
    
    /*
     * Helpers.
     */
    mustSample = 0;             // default: false, do not sample
    zPars.assign( 10, 1.0 );    // long enough
    gMsg.assign( 7, 0.0 );
} // end GUTS::GUTS()

/*
 * Destructor.
 */
GUTS::~GUTS()
{
    mC.clear();
    mCt.clear();
    my.clear();
    myt.clear();
    mpar.clear();
    mz.clear();
    D.clear();
    S.clear();
    ee.clear();
    ff.clear();
    zPars.clear();
    gMsg.clear();
} // end GUTS::~GUTS()

/*
 * mC and mCt
 */
void GUTS::setConcentrations(const vector<double>& C, const vector<double>& Ct)
{
    // Minimum length, same size, first Ct = 0
    if (C.size() > 2 && C.size() == Ct.size() && Ct.at(0) == 0.0)
    {
        // No Ct.at(k+1) - Ct.at(k) <= 0! Prevent division by 0.
        for ( unsigned int i = 0; i < (Ct.size() - 1); ++i)
        {
            if ( (Ct.at(i+1) - Ct.at(i)) <= 0 )
            {
                gMsg.at(0) = 0.0;
                return;
            }
        }

        // Okay, lets assign
        mC = C;
        mCt = Ct;
        gMsg.at(0) = 1.0;
    }
    else
    {
        gMsg.at(0) = 0.0;
    }
} // end GUTS::setConcentrations()

/*
 * Setter for my and myt. my needs a 0 appended. S.back() must be 0.0.
 * S and dtau are also set.
 */
void GUTS::setSurvivors(const vector<int>& y, const vector<double>& yt)
{
    // Minimum length, same size, first yt = 0, last > first
    if (y.size() > 2 && y.size() == yt.size() && yt.at(0) == 0.0)
    {
        // No yt.at(k+1) - yt.at(k) <= 0!
        for ( unsigned int i = 0; i < (yt.size() - 1); ++i)
        {
            if ( (yt.at(i+1) - yt.at(i)) <= 0 )
            {
                gMsg.at(1) = 0.0;
                return;
            }
        }

        my = y;
        my.push_back( 0 ); // my[last + 1] is always 0
        myt = yt;

        S.resize( my.size() );
        S.back() = 0.0;

        dtau = ( myt.back() - myt.front() ) / mM;

        gMsg.at(1) = 1.0;
    }
    else
    {
        gMsg.at(1) = 0.0;
    }
} // end GUTS::setSurvivors()

/*
 * Setter for mpar, zPars.
 * zPars are checkt in doSample, according to sample to be drawn.
 */
void GUTS::setParameters(const vector<double>& par)
{
    // at least 3 elements required.
    if ( par.size() > 2 )
    {
        /*
         * If any of the first three pars is lower than 0.000...
         * was: std::numeric_limits<double>::epsilon(), but no need for
         * this precision.
         * Second par must be > 0
         * Check zPars in doSample, not here.
         */
        if ( *min_element( par.begin(), par.begin()+3 ) < 0.0 )
        {
            gMsg.at(2) = -1000.0;
            return;
        }
        else if ( par.at( 1 ) <= 0.0 )
        {
            gMsg.at(2) = 0.0;
            return;
        }

        mpar = par;
        gMsg.at(2) = 1.0;
        
        // Overwrite zPar defaults, if mpar.size() > 3
        for ( unsigned int i = 3; i < mpar.size(); ++i )
        {
            mustSample = 1; // no extra if-check, just say do sample
            zPars.at( i - 3 ) = mpar.at( i ); // 0 <- 3 ...
        }
    }
    else
    {
        gMsg.at(2) = 0.0;
    }
} // end GUTS::setParameters()

/*
 * Setter for mM, D, and dtau.
 */
void GUTS::setTimeGridPoints(const int& M)
{
    // M > 4 to be safe.
    if (M > 4)
    {
        mM = M;
        D.assign( mM, 0.0 ); // default
        dtau = ( myt.back() - myt.front() ) / mM;

        gMsg.at(3) = 1.0;
    }
    else
    {
        gMsg.at(3) = 0.0;
    }
} // end GUTS::setTimeGridPoints()

/*
 * Set distribution.
 */
void GUTS::setDistribution( const string dist )
{
    if ( dist == "lognormal" )
    {
        mdist = dist;
        mustSample = 1; // true, do sample
        gMsg.at(4) = 1.0; // no error
    }
    else if ( dist == "empirical" )
    {
        mdist = dist;
        mustSample = 0; // false, do not sample
        gMsg.at(4) = 1.0; // no error
    }
    else
    {
        mdist = "none"; // the default
        gMsg.at(4) = 0.0; // error
        mustSample = 1; // true, fill all with 1.1
    }
} // end GUTS::setDistribution()

/*
 * Setter for mN.
 */
void GUTS::setSampleLength(const int& N)
{
    /*
     * N > 4 to be safe. Sizes of ee and ff are N-dependent, but
     * these two vectors are reset in calcLoglikelihood anyway.
     */
    if (N > 4)
    {
        mN = N;
        gMsg.at(5) = 1.0;
    }
    else
    {
        gMsg.at(5) = 0.0;
    }
} // end GUTS::setSampleLength()

/*
 * The user method for passing a ready to run vector.
 */
void GUTS::setSample( const vector<double>& z )
{
    if ( z.size() > 4 )
    {
        mustSample = 0; // false, do not sample
        mz = z;
        sort( mz.begin(), mz.end() ); // Always sort!
        setDistribution( "empirical" );
        setSampleLength( mz.size() );
        doSample(); // setting error there to 0
        gMsg.at(6) = 1.0;
    }
    else
    {
        gMsg.at(6) = 0.0;
    }
} // end GUTS::setSample()

/*
 * Printout. Currently not very nice :-).
 */

void to_print( double i )
{
    cout << i << ", ";
}
void GUTS::showObject()
{
    std::cout << left << "Concentrations (C, " << mC.size()
        << " elements):" << endl;
    for_each( mC.begin(), mC.end()-1, to_print);
    cout << mC.back() << endl;
    
    cout << left << "Concentration times (Ct, " << mCt.size()
        << " elements):" << endl;
    for_each( mCt.begin(), mCt.end()-1, to_print);
    cout << mCt.back() << endl;

    cout << left << "Suvivors (y, " << myt.size()
        << " elements):" << endl;
    for_each( my.begin(), my.end()-2, to_print);
    cout << my.at(my.size()-2) << endl;

    cout << left << "Suvivor times (yt, " << myt.size()
        << " elements):" << endl;
    for_each( myt.begin(), myt.end()-1, to_print);
    cout << myt.back() << endl;
    
    cout << left << "Parameters (par, " << mpar.size()
        << " elements):" << endl;
    for_each( mpar.begin(), mpar.end()-1, to_print);
    cout << mpar.back() << endl;
    

    cout << left << "Time grid points (M): " << mM << endl;
    cout << left << "Distribution (dist) : " << mdist << endl;
    cout << left << "Sample length (N)   : " << mN << endl;
/*    
    cout << "mz: first: " << mz.front() << ", last: " << mz.back()
        << ", size: " << mz.size() << endl;
*/
/*
    cout << left << "gMsg (" << gMsg.size()
        << " elements):" << endl;
    for_each( gMsg.begin(), gMsg.end()-1, to_print);
    cout << gMsg.back() << endl;
*/
    cout << "\n" << endl;
}
