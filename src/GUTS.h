/** GUTS class.
 * @file    GUTS.h
 * @author  carlo.albert@eawag.ch, soeren.vogel@eawag.ch
 * @date    30 Jun 2011
 *
 * C++ class GUTS:
 * Fast Calculation of the Likelihood of a Stochastic Survival Model
 * License: GPL-2
 */

#ifndef guts_h
#define guts_h guts_h

#include <vector>
#include <valarray>
#include <numeric>
#include <algorithm>
#include <string>
#include <map>
#include <ctime> // for rng
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>

using namespace std;

/** Class GUTS
 * GUTS objects represent survival models for which the logarithm
 * of the likelihood can be calculated.
 */
class GUTS
{
public:

    GUTS();
    ~GUTS();
    
    /** Time series vector of concentrations.
     * @param   C   Vector of concentrations.
     * @param   Ct  Vector of time points of concentrations (must start at 0).
     *
     * Must be of equal length. Set mC and mCt.
     */
    void setConcentrations(const vector<double>& C, const vector<double>& Ct);

    /** Time series vector of survivors.
     * @param   y   Vector of survivors.
     * @param   yt  Vector of time points of survivors (must start at 0).
     *
     * Must be of equal length. Set y/my and yt/myt, appends value 0 to my.
     */
    void setSurvivors(const vector<int>& y, const vector<double>& yt);

    /** Parameter vector of the object.
     * @param   par     Vector of numerics (> 2).
     * @return  True/false where false indicates an error.
     *
     * Vector must be of size >2. First three are required for the calculation
     * of the loglikelihood, the latter are used for the construction of sample
     * vector z. Set par/mpar.
     */
    void setParameters(const vector<double>& par);

    /** Number of grid points on the time-axis.
     * @method  setTimeGridPoints
     * @param   M   Positive non-zero integer.
     * @return  True/false where false indicates an error.
     *
     * Indicates the grid points of the time axis used for the integration.
     * Set M/mM.
     */
    void setTimeGridPoints(const int& M);

    /** The (name of the) distribution to sample from.
     * @param:  dist    character (name)
     * @return: True/false where false indicates an error.
     */
    void setDistribution( const string dist );

    /** Size of the sample.
     * @param   N   Positive non-zero integer.
     * @return  True/false where false indicates an error.
     */
    void setSampleLength(const int& N);

    /** Sample vector.
     * @param   z    Vector of non-negative ordered numerics.
     * @param   byUser  Sample provided by user or from class parameters.
     * @return  True/false where false indicates an error.
     *
     * Public method takes a vector of non-negative numerics. The vector must
     * be in ascending order. Sets mz and hbyUser.
     */
    void setSample( const vector<double> &z);

    /** Calculate loglikelihood.
     * @return  Numeric, the loglikelihood.
     *
     * Public method takes all private fields of the object and calculates
     * the loglikelihood. Various consistency checks are done, see
     * implementation for details. Returns -1000 per default, else a
     * negative numeric.
     */
    double calcLoglikelihood();

    /*
     * for debugging only
     */
    void showObject();
    
    /*
     * Getters.
     */
    vector<double> getC() { return mC; }
    vector<double> getCt() { return mCt; }
    vector<int> gety() { return my; }
    vector<double> getyt() { return myt; }
    vector<double> getpar() { return mpar; }
    int getM() { return mM; }
    string getdist() { return mdist; }
    int getN() { return mN; }
    vector<double> getz() { return mz; }

private:

    /** Construct the sample, sort, and save result in z.
     * @param   none
     * @return  ordered sample of length mN from mdist saved in mz
     *
     * The doSample method uses algorithms from the boost library.
     */
    bool doSample();

    /** Attributes of a GUTS object.
     *
     * Each variable from the methods above is prefixed with character m.
     * See method description for details on the variables.
     */
    vector<double> mC;      // Concentrations
    vector<double> mCt;     // Time points of concentratinos
    vector<int> my;         // Survivors (= y, appended 0)
    vector<double> myt;     // Survivor time points
    vector<double> mpar;    // Parameters
    int mM;                 // Time grid points
    string mdist;           // Name of distribution
    int mN;                 // Length of z (sample)
    vector<double> mz;      // z (sample)

    /** Current state of a GUTS object.
     *
     * Fields represent the current state of the object.
     * D: time series of damages in the organism
     * S: time series of survivor probabilities
     * ee: sum of all Ds between z_j and z_j+1
     * ff: number of Ds between z_j and z_j+1
     */
    vector<double> D;       // length of mM, first element 0.0
    vector<double> S;       // length of my, any elements
    vector<double> ee;      // length of mN, all elements 0.0
    vector<int> ff;         // length of mN, all elements 0
    double dtau;            // time difference of myt divided by mM
    
    /*
     * Helpers.
     */
    bool mustSample;        // Do we need to sample?
    vector<double> zPars;   // Parameters for sampling
    vector<double> gMsg;
    //std::map<std::string, int> gError;
};

#endif
