/*
 *  PulseShape.h
 *  C_ppu
 *
 *  Created by vchaplin on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PULSESHAPE_H
#define PULSESHAPE_H

#include <math.h>
#include "ChannelEnergy.h"

// !! NEW... required for GSL splines.
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

	
const double dpu_sampling_MHz = 9.6;
const double dpu_period_mus = 1.0 / dpu_sampling_MHz; //period in micro-seconds
const double dpu_peak_wait = 4.0*dpu_period_mus;


namespace paramSet0 {
	//Parameters from N.B.; Jan 2012
	const double a = 3;
	const double b = 2;
	const double c = 1.6;
	const double al = 1.25;
	const double be = 2;
	const double pulse_max = 0.342317;
	const double pulse_maxinv = 1.0/pulse_max;
	const double t0p = 0.484861; //Peak time of a single pulse; computed seperately
	const double positivePulseWidth = 1.71707; //Total width of v(t) > 0
    
    const double t0minimum = 2.43765;
    
    const double t1Offset = 0.01;
    const double t1VoltConst = 2.0 / 9.0;
    const double t1GrowthConst = c*c/7.0;
	
    
	const double TailStart = 2.6;
    
	const double TailEnd = 5.0;
	
	const double TauPEAK = t0p + dpu_peak_wait;
    //       const double TailStart = 0.9;
	const double TauTAIL = TailEnd - TailStart;
	const double TauDEAD = TailStart - TauPEAK;
};

namespace paramSet1 {
	//Parameters from N.B.; March 2012
	const double a = 28;
	const double b = 19;
	const double c = 2.6;
	const double al = 1.2;
	const double be = 3.5;
	const double pulse_max = 3.026147451168;
	const double pulse_maxinv = 1.0/pulse_max;
	const double t0p = 0.387896642254; //Peak time of a single pulse; computed seperately
	const double positivePulseWidth = 1.1836; //Total width of v(t) > 0
    
    const double t0minimum = 1.844648;
    
    const double t1Offset = 0.01;
	const double t1VoltConst = 2.0 / 9.0;
	const double t1GrowthConst = c*c/16.0;
	

	const double TailStart = 2.6;
    
	const double TailEnd = 5.0;
	
	const double TauPEAK = t0p + dpu_peak_wait;
 //       const double TailStart = 0.9;
	const double TauTAIL = TailEnd - TailStart;
	const double TauDEAD = TailStart - TauPEAK;
};

namespace paramSet2 {
	//Parameters from N.B.; June 2, 2012
	const double a = 26;
	const double b = 31;
	const double c = 2.6;
	const double al = 1.27;
	const double be = 3.5;
	const double pulse_max = 2.446026964264095;
	const double pulse_maxinv = 1.0/pulse_max;
	const double t0p = 0.3649115034990959; //Peak time of a single pulse; computed seperately
	const double positivePulseWidth = 0.924155678; //Total width of v(t) > 0
    
    const double t0minimum = 1.6630257;
    
    const double t1Offset = 0.00;
    const double t1VoltConst = 2.09;//2.0 / 9.0;
    const double t1GrowthConst = 0.40;//c*c/16.0;
	
    
	const double TailStart = 2.6;
    
	const double TailEnd = 4.5;
	
	const double TauPEAK = t0p + dpu_peak_wait;
    //       const double TailStart = 0.9;
	const double TauTAIL = TailEnd - TailStart;
	const double TauDEAD = TailStart - TauPEAK;
};



/*****************************

	New parameter set with your values
	
*****************************/
namespace NarayanParams {

	//whatever's necessary for the spline data...
	double nodesX[] = { 0.0, 0.1, 0.2, }; //time will be passed in microseconds so the spline nodes should be in us
	double nodesY[] = { 0.0, .... }; //Y amplitude is in arbitrary units
	const long numNodes = .. ; //number of nodes
	
	//now construct spline functions using GSL...
	static gsl_interp_accel * accelator = gsl_interp_accel_alloc ();
    static gsl_spline * myspline = gsl_spline_alloc (gsl_interp_cspline, numNodes);
	gsl_spline_init (myspline, nodesX, nodesY, numNodes);
	//I think that's enough for the spline data to persist between function calls.  
	
	
	
	//put in the new values where there blank...
	
	//whatever the units of nodesY, 'pulse_max' must be the max-value of the pulse shape in these units
	//The pulse is normalized such that it equals one at the maximum, ie: 'pulse(t0p) = 1.0'
	
	const double pulse_max = ;  //<----  max value of the interpolating function.  Will be scaled for different energies
	
	const double pulse_maxinv = 1.0/pulse_max;
	const double t0p =  ; // <----  Peak time of a single pulse; computed seperately
	
	//these are the parameters for the function t1p() below.  They will likely need to change but not sure how
    const double t1Offset = 0.00;
    const double t1VoltConst = 2.09;
    const double t1GrowthConst = 0.40;
	
    //TailStart is actually equal to the deadtime, so if that's still 2.6us use that
	const double TailStart = 2.6; //micro seconds
	
	//TailEnd = end point of a single pulse (somewhat arbitrary...probably between 4.5us - 6us)
	//  In the paper I found that the quality of the model depends somewhat on this value.  In
	//  general, if it looks like there's too much tail effect in the output spectrum 
	//  (aka too many low energy counts), then TailEnd might be too large.  If there's not enough
	//  low energy counts, TailEnd might be too small.  The analytical method is formally an 
	//  approximation so the question is: what reasonable value of TailEnd best suits the assumptions.
	
	const double TailEnd = ; // <----    microseconds
	
	//THis is all unchanged
	const double TauPEAK = t0p + dpu_peak_wait;
	const double TauTAIL = TailEnd - TailStart;
	const double TauDEAD = TailStart - TauPEAK;
	
	//These aren't used anymore but probably need to be defined in order to compile easily
	const double a = 26;
	const double b = 31;
	const double c = 2.6;
	const double al = 1.27;
	const double be = 3.5;
};




// changed to the name of the new namespace to use new parameters
#define PARAM_SET NarayanParams  


using namespace PARAM_SET;

/**************
 Possible example of new pulse() function
*/
inline double pulse( double time ) 
{
	if (time <= 0.0) {
		return 0.0;
	}
	
	//now pulse() will call the spline evaluation function.
	//The pulse_maxinv is there so the pulse is a "unit" pulse (max value of 1.0)
	
	return pulse_maxinv * gsl_spline_eval (myspline, time, accelator);
};
/**************
 OLD pulse() function
*/
/*
inline double pulse( double time ) 
{
	if (time <= 0.0) {
		return 0.0;
	}
	return pulse_maxinv * exp( -c*time ) * (a*pow(time, al) - b*pow(time, be));
};
*/

inline double t1p( double v0, double v1, double s1 )
{
	
	double t0_avg = (v0*t0p + v1*(t0p+s1)) / ( v0+v1 );
    //if ((v0/v1) >= 3.0 && (t0_avg < s1)) return t0p;
	
	double t1_approx  = t0_avg + 0.0 + exp(-t1VoltConst * pow((v0-v1)/(v0+v1), 2) ) * (exp( t1GrowthConst * pow(s1,2) ) - 1.0 );
	
	
	/*
	if ( (s1 > t1_approx) || (s1 > TauPEAK)) return(t0p);
    if ( (v1/v0) >= 3.0 ) return t0p+s1;
	if ((v0/v1) >= 3.0 && (t1_approx > TauPEAK)) t1_approx = t0p;
	*/
    
    /*
     12 dec 2012 - removed these and moved to PeakE function in MeasredEnergies.cpp
    if ( t1_approx < t0p*exp(1.7*s1*s1) ) return (t0p);
    if ( s1 > TauPEAK ) return (t0p);
     */
    
	return t1_approx;
};


//These two functions, 't1p_dt' and 'pulse_dt', aren't used anymore but might need to 
//be here so everything compiles easily
inline double t1p_dt( double v0, double v1, double s1 )
{
	return v1/(v0+v1) + exp(-t1VoltConst*pow((v0-v1)/(v0+v1), 2)) * exp( t1GrowthConst * pow(s1,2) ) * 2*s1*t1GrowthConst;	
};


inline double pulse_dt( double time )
{
	if (time <= 0.0) {
		return 0.0;
	}
	
	return pulse_maxinv * exp( -c*time ) * ( 
				a*al*pow(time, al-1) - b*be*pow(time, be-1) - c*(a*pow(time, al) - b*pow(time, be)) 
									   );
};

#endif

