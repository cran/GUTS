/** GUTS method doSample
 * @file    GUTS_doSample.cpp
 * @author carlo.albert@eawag.ch, soeren.vogel@uzh.ch
 * @date 30 Jun 2011
 * @license GPL-2
 *
 * Private method for generating the sample, if done from GUTS.
 * ee and ff need to have the length of the sample, however, they are
 * reset to the length mN in calcLoglikelihood anyway, so no need here.
 */

#include "GUTS.h"

using namespace std;

/*
 * Declaration of sample functions to invoke in doSample.
 */
std::vector<double> sampleNormal     ( int N, double mean, double sigma );
std::vector<double> sampleLognormal  ( int N, double mean, double sigma );

/*
 * The method.
 */
bool GUTS::doSample()
{
    /*
     * Choose distribution and sample.
     * Currently only normal + lognorm, more future work.
     */
    if ( mdist == "normal" )
    {
        if ( zPars.at(1) <= 0 )
            return 0;
        else
            mz = sampleNormal( mN, zPars.at(0), zPars.at(1) );
    }
    else if ( mdist == "lognormal" )
    {
        if ( zPars.at(0) <= 0 || zPars.at(1) <= 0 )
            return 0;
        else
            mz = sampleLognormal( mN, zPars.at(0), zPars.at(1) );
    }
    else
    {
        mz.assign( mN, 0.0);
        return 0;
    }

    /*
     * Always sort!
     */
    sort( mz.begin(), mz.end() );
    return 1;
} // end GUTS::doSample()

/*
 * Sample functions implementations.
 */

using namespace boost;

// Create a Mersenne twister random number generator
// that is seeded once with #seconds since 1970
static mt19937 rng(static_cast<unsigned> (std::time(0)));

std::vector<double> sampleNormal (int N, double mean, double sigma)
{
//    using namespace boost;

    /*
     * temporary vector
     */
    vector<double> ret;

    // Create a Mersenne twister random number generator
    // that is seeded once with #seconds since 1970
//    static mt19937 rng(static_cast<unsigned> (std::time(0)));

    // select Gaussian probability distribution
    normal_distribution<double> norm_dist(mean, sigma);

    // bind random number generator to distribution, forming a function
    variate_generator<mt19937&, normal_distribution<double> >
        normal_sampler(rng, norm_dist);
    
    // sample from the distribution
    for ( unsigned int i = 0; i < N; ++i )
    {
        ret.push_back( normal_sampler() );
    }
    return ret;
}

std::vector<double> sampleLognormal (int N, double mean, double sigma)
{
    vector<double> ret;
    lognormal_distribution<double> lognorm_dist(mean, sigma);
    variate_generator<mt19937&, lognormal_distribution<double> >
        lognormal_sampler(rng, lognorm_dist);
    for ( unsigned int i = 0; i < N; ++i )
    {
        ret.push_back( lognormal_sampler() );
    }
    return ret;
}
