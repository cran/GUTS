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
    my.assign( 11, 0 ); // appended field with 0
    myt.assign( 10, 1.0 );
    mpar.assign( 5, 1.0 );
    mM = 10000;
    mdist = "lognormal";
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
    dtau = ( myt.back() - myt.front() ) / mM;
    
    /*
     * Helpers.
     */
    sampleByUser = 0; // default: false, do sample
    zPars.assign( 10, 1.0 ); // long enough

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
} // end GUTS::~GUTS()

/*
 * Setter of mC and mCt
 */
bool GUTS::setConcentrations(const vector<double>& C, const vector<double>& Ct)
{
    // Minimum length, same size, first Ct = 0
    if (C.size() > 2 && C.size() == Ct.size() && Ct.at(0) == 0.0)
    {
        // No Ct.at(k+1) - Ct.at(k) <= 0! Prevent division by 0.
        for ( unsigned int i = 0; i < (Ct.size() - 1); ++i)
        {
            if ( (Ct.at(i+1) - Ct.at(i)) <= 0 )
                return 0;
        }

        // Okay, lets assign
        mC = C;
        mCt = Ct;
        return 1;
    }
    else
    {
        return 0;
    }
} // end GUTS::setConcentrations()

/*
 * Setter for my and myt. my needs a 0 appended. S.back() must be 0.0.
 * S and dtau are also set.
 */
bool GUTS::setSurvivors(const vector<int>& y, const vector<double>& yt)
{
    // Minimum length, same size, first yt = 0, last > first
    if (y.size() > 2 && y.size() == yt.size()
        && yt.at(0) == 0.0 && yt.back() > yt.front() )
    {
        my = y;
        myt = yt;
        my.push_back( 0 ); // my[last + 1] is always 0
        S.resize( my.size() );
        S.back() = 0.0;
        dtau = ( myt.back() - myt.front() ) / mM;
        return 1;
    }
    else
    {
        return 0;
    }
} // end GUTS::setSurvivors()

/*
 * Setter for mpar, zPars.
 * zPars are checkt in doSample, according to sample to be drawn.
 */
bool GUTS::setParameters(const vector<double>& par)
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
        if ( *min_element( mpar.begin(), mpar.begin()+3 ) < 0.0 )
            return 0;
        else if ( par.at( 1 ) <= 0.0 )
            return 0;

        mpar = par;
        
        // Overwrite zPar defaults, if mpar.size() > 3
        for ( unsigned int i = 3; i < mpar.size(); ++i )
        {
            sampleByUser = 0; // no extra if-check, just do it 2 times or so
            zPars.at( i - 3 ) = mpar.at( i ); // 0 <- 3 ...
        }
        return 1;
    }
    else
    {
        return 0;
    }
} // end GUTS::setParameters()

/*
 * Setter for mM, D, and dtau.
 */
bool GUTS::setTimeGridPoints(const int& M)
{
    // M > 4 to be safe.
    if (M > 4)
    {
        mM = M;
        D.assign( mM, 0.0 ); // default
        dtau = ( myt.back() - myt.front() ) / mM;
        return 1;
    }
    else
    {
        return 0;
    }
} // end GUTS::setTimeGridPoints()

/*
 * Set distribution.
 */
bool GUTS::setDistribution( const string dist )
{
// FIXME: kann dict abfragen, ob dist enthalten ist!!!
/*
    if ( mdist == "" || mdist == "none" || mdist == "empirical" )
    {
        sampleByUser = 1;
    }
*/
    // if not in dict, return -1

    mdist = dist;
    sampleByUser = 0;
    return 1;
} // end GUTS::setDistribution()

/*
 * Setter for mN.
 */
bool GUTS::setSampleLength(const int& N)
{
    /*
     * N > 4 to be safe.
     * Sizes of ee and ff are N-dependent, but these two vectors are
     * reseted in calcLoglikelihood anyway.
     */
    if (N > 4)
    {
        mN = N;
        sampleByUser = 0;
        return 1;
    }
    else
    {
        return 0;
    }
} // end GUTS::setSampleLength()

/*
 * The user method for passing a ready to run vector.
 */
bool GUTS::setSample(const vector<double>& z )
{
    if ( z.size() > 4 )
    {
        sampleByUser = 1; // remember that the user sampled already
        mz = z;
        sort( mz.begin(), mz.end() ); // Always sort!
        mN = mz.size();
        mdist = "empirical";
        return 1;
    }
    else
    {
        return 0;
    }
} // end GUTS::setSample()

/*
 * for debugging only
 */
void GUTS::showObject()
{
    cout << "mC: first: " << mC.front() << ", last: " << mC.back()
        << ", size: " << mC.size() << endl;
    cout << "mCt: first: " << mCt.front() << ", last: " << mCt.back()
        << ", size: " << mCt.size() << endl;
    cout << "my: first: " << my.front() << ", last: " << my.back()
        << ", size: " << my.size() << endl;
    cout << "myt: first: " << myt.front() << ", last: " << myt.back()
        << ", size: " << myt.size() << endl;
    cout << "mM: " << mM << "\nmN: " << mN << endl;
    cout << "dist: " << mdist << endl;
    cout << "mpar: first: " << mpar.front() << ", last: " << mpar.back()
        << ", size: " << mpar.size() << endl;
    cout << "mz: first: " << mz.front() << ", last: " << mz.back()
        << ", size: " << mz.size() << endl;
    cout << "D: first: " << D.front() << ", last: " << D.back()
        << ", size: " << D.size() << endl;
    cout << "S: first: " << S.front() << ", last: " << S.back()
        << ", size: " << S.size() << endl;
    cout << "ee: first: " << ee.front() << ", last: " << ee.back()
        << ", size: " << ee.size() << endl;
    cout << "ff: first: " << ff.front() << ", last: " << ff.back()
        << ", size: " << ff.size() << endl;
    cout << "dtau: " << dtau << endl;
}
