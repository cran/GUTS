/** GUTS method calcLoglikelihood
 * @file    GUTS_calcLoglikelihood.cpp
 * @author  carlo.albert@eawag.ch, soeren.vogel@uzh.ch
 * @date    30 Jun 2011
 * @license GPL-2
 *
 * Calculate the logarith of the likelihood of data present in GUTS object.
 */

#include "GUTS.h"

using namespace std;

double GUTS::calcLoglikelihood()
{
    double outDefault = -1000.0; // the default return

    /*
     * Check some errors.
     */
    if ( myt.back() > mCt.back() )
        return 0.0; //exit( 1 ); // error
    if ( dtau <= 0.0 )
        return outDefault;
    if ( doSample() == 0 )
        return 0.0;
    for ( unsigned int i = 0; i < gMsg.size(); ++i )
        if ( gMsg.at(i) < 1.0 ) return gMsg.at(i); // any error
//    if ( *min_element( gMsg.begin(), gMsg.end() ) < 0 )
//        return outDefault;
    /*
     * End error check.
     */

    /*
     * Reset vectors.
     */
    D.at( 0 ) = 0.0;
    ee.assign( mN, 0.0 );
    ff.assign( mN, 0 );

    /*
     * Iterators:
     *
     * j: iterator for D, ends at mM
     * ii: iterator for D, is the last j
     * zpos: current position in z
     * tau: 0 -- max(yt), 10000 elements, incremented by value of dtau
     * k: index for C
     * i: iterator for yt
     */
    int j = 0;
    int ii = 0;
    int zpos = 0;
    double tau = 0.0;
    double dk = mpar.at(2) * dtau; // killing rate
    int k = 0;

    /*
     * Loop over yt_i's.
     */
    for ( unsigned int i = 0; i < myt.size(); ++i )
    {
        /*
         * must control for tau
         * must control for j to not exceed dim
         * of D due to calculation inaccuracy
         */
        while ( tau < myt.at(i) && j < mM )
        {
            // Calculate D
            double tmp = exp( -mpar.at(1) * ( tau - mCt.at(k) ) );
            D.at(j) = D.at(ii) * tmp
                      + mC.at(k) * ( 1-tmp )
                      + ( mC.at(k+1) - mC.at(k) )
                      * ( tau - mCt.at(k) - ( 1 / mpar.at(1) ) * ( 1-tmp ) )
                      / ( mCt.at(k+1) - mCt.at(k) );

            // Increment or decrement zpos
            if ( D.at(j) < mz.at(zpos) )
            {
                while ( D.at(j) < mz.at(zpos) && zpos > 0 )
                {
                    --zpos;
                }

                /*
                 * zpos may have changed, so mz.at(zpos) may also
                 * D.at(j) can now be larger than mz.at(zpos)
                 */
                if ( D.at(j) > mz.at(zpos) )
                {
                    ++zpos;
                }
            }
            else
            {
                while ( D.at(j) > mz.at(zpos) && zpos < (mN-1) )
                {
                    ++zpos;
                }
            }

            // Update ee and ff
            if ( D.at(j) > mz.at(mN-1) ) //
            {
                ee.at(mN-1) += D.at(j);
                ff.at(mN-1)++;
            }
            else if ( D.at(j) > mz.at(0) )
            {
                ee.at(zpos-1) += D.at(j);
                ff.at(zpos-1)++;
            }

            // Increment or decrement j, tau, k
            j++;
            tau += dtau;
            if ( tau > mCt.at(k+1) )
            {
                k++;
                ii = j-1;
            }
        } // end while ( tau < myt.at(i) )

        // E = sum of all ee
        // F = sum of all ff
        double E = accumulate( ee.begin(), ee.end(), 0.0);
        double F = accumulate( ff.begin(), ff.end(), 0.0);
        
        // Calculate S.at(i):
        S.at(i) = exp( dk * ( mz.at(0) * F - E ) );
        for ( int u=1; u < mN; ++u )
        {
            E -= ee.at(u-1);
            F -= ff.at(u-1);
            S.at(i) += exp( dk * ( mz.at(u) * F - E ) );
        }
        S.at(i) *= exp( -mpar.front() * myt.at(i) ) / mN;

    } // end for ( unsigned int i = 0; i < myt.size (); i++ )

    // Calculate loglikelihood
    double out = 0.0;
    for ( unsigned int i=0; i < myt.size(); ++i )
    {
        double diffS = S.at(i) - S.at(i+1); // always positive
        int diffy = my.at(i) - my.at(i+1);
        if ( diffS < 0.001 )
        {
            if ( diffy != 0 ) return outDefault;
        }
        else
        {
            out += diffy * log( diffS );
        }
    } // end caluclate loglikelihood

    return out;
} // end GUTS::calcLoglikelihood()
