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


// other parameter values possible if they're in a namespace, e.g. 'namespace paramSet2 {...};'
// PARAM_SET macro defines which is used

#define PARAM_SET paramSet2
using namespace PARAM_SET;

inline double pulse( double time ) 
{
	if (time <= 0.0) {
		return 0.0;
	}

	return pulse_maxinv * exp( -c*time ) * (a*pow(time, al) - b*pow(time, be));
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

inline double t1p_dt( double v0, double v1, double s1 )
{
	return v1/(v0+v1) + exp(-t1VoltConst*pow((v0-v1)/(v0+v1), 2)) * exp( t1GrowthConst * pow(s1,2) ) * 2*s1*t1GrowthConst;	
};

#endif

