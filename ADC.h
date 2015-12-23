//
//  ADC.h
//  C_ppu
//
//  Created by vchaplin on 1/23/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef C_ppu_ADC_h
#define C_ppu_ADC_h

#include "PulseShape.h"

#include <iostream>
#include <cstdlib>
#include <math.h>

using namespace std;


inline float summedPulse( float& vmin, float& vmax, float * times, float * volts, size_t& num, float& t, float& t0 )
{
    size_t pk=0;
    float sum=0;
    float summand;
    
    while ( pk < num ) {
        
        // cout << pk << " " << pulse( t - (times[pk] - t0) ) << "  tn: " <<(times[pk] - t0)*1e6 << " vn: "  << volts[pk] << endl;
        
        summand = volts[pk]*pulse( (t - (times[pk] - t0))*1e06 );
        if (summand < vmin) {
            summand = vmin;
        } else if (summand > vmax) {
            summand = vmax;
        }
        sum += summand;
        
        pk++;
    }
    
    //if ( sum > vmax ) sum = vmax;
    //else 
    //if ( sum < vmin ) sum = vmin;
    
    return sum;
}





/*******
 
 FUNCTION:  adcMeasurement
 
 Takes a list of input times, voltages, and emulates ADC / PHA measurement done by GBM.  Returns a list of measured times and pulse heights.
 The output arrays must be allocated and have enough space to hold the counted pulses ( <= # input events ).  The number counted is not known exactly.
 
 ARGUMENTS:
 
 threshold - voltage LLT, below which no pulse-heights are measured (GBM = 0.0104 Volts)
 
 vfloor - lower clipping voltage (min.) value for a single pulse (not the sum) [GBM = -0.5 V]
 
 vmax - upper clipping voltage [GBM= 5.0 V]
 
 inputTimes - array of times in SECONDS of the input events

 inputVolts - array of peak-pulse hieghts (in volts) of the input events
 
 numPulse - number of input events
 
 outputTimes - measured times
 
 outputPulseHeights - measured pulse heights (volts)
 
 dSt - ADC sampling period (default = 50 ns, about half of GBM's value)
 
 PostPeakBuffer - size of additional deadtime after the 4-sample buffer (default = 0.104 * 21 ).
 
 
 RETURN VALUE:
 
 Number of measured counts (# valid elements in output arrays)
 
********/


size_t adcMeasurement( float& threshold, float& vfloor, float& vmax, float * inputTimes, float * inputVolts, size_t numPulses, float * outputTimes, float * outputPulseHieghts, float dSt = 5e-08, float PostBufferTime=2.1e-06 );


/*
 
 
 inputTimes - input times in seconds
 inputVolts - voltage of the input events ( overflow channel=5V, min channel = 0V)
 
 numInputs - number of input events
 
 
 window_width = 6e-06 seconds (this is the default)
 
 
 
 */

size_t ProcessPulsePileup( float *& inputTimes, float *& inputVolts, size_t numInputs, float *& outputTimes, float *& outputVolts, float window_width = 6.0e-06, float threshold=0.01, float vfloor=-0.5, float vmax = 5.0 );

size_t ProcessPPUtoQuota( size_t detectionQuota, size_t& numNeeded, float *& inputTimes, float *& inputVolts, size_t numInputs, float *& outputTimes, float *& outputVolts, float window_width = 6.0e-06, float threshold=0.01, float vfloor=-0.5, float vmax = 5.0   );

size_t FilterDeadtime( float *& inputTimes, float *& inputPH, size_t numIn, float *& outputTimes, float *& outputPulseHeights, float deadTime=2.6e-06 );






#endif
