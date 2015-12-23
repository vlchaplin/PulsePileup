//
//  ADC.cpp
//  C_ppu
//
//  Created by vchaplin on 1/23/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "ADC.h"


size_t adcMeasurement( float& threshold, float& vfloor, float& vmax, float * inputTimes, float * inputVolts, size_t numPulses, float * outputTimes, float * outputPulseHieghts, float dSt, float PostBufferTime ) 
{
    float tr = inputTimes[0];
    float t0 = tr;//nearest(tr, dSt);
    float tmax = inputTimes[numPulses-1] + (1e-06)*6 - t0;
    float ts,Va,Vb, tp, Vpeak;
    size_t outK=0;
    bool have_peak=0;
    bool rising=1;
    
    Vb=0;
    ts=0;
    
    //cout << "Beginning ADC loop @ t=" << tr <<endl;;
    
    while( ts < tmax ) {
        
        Va = summedPulse(vfloor, vmax, inputTimes, inputVolts, numPulses, ts, t0);
        
        if ( Va < threshold || Vb < threshold ) {
            ts += dSt;
            Vb=Va;
            continue;
        }
        
        if ( Va > Vb ) {
            Vb = Va;
            ts += dSt;
            rising=1;
            continue;
        }
        
        if ( Va < Vb && rising ) {
            //cout << "peak" << endl;
            tp = ts -dSt;
            Vpeak = Vb;
            tr=0;
            rising=0;
            have_peak=1;
        }
        //cout << ts << endl;
        while ( have_peak && Va < Vb && tr < (1e-06)*dpu_peak_wait )
        {
            Vb = Va;
            ts += dSt;
            tr += dSt;
            Va = summedPulse(vfloor, vmax, inputTimes, inputVolts, numPulses, ts, t0);
            
            //if ( Va - Vb  > 0 ) cout << "next count..." << endl;
        }
        
        
        
        if (have_peak && tr >= (1.0e-06)*dpu_peak_wait ) {
            if( Vpeak >= threshold ){
                outputTimes[outK] = tp + t0;//- (1e-06)*PARAM_SET::t0p;
                outputPulseHieghts[outK] = Vpeak;
                outK++;
            }
            
            if ( Vpeak >= 5.0 ) {
                ts += 10e-06;
            } else {
                ts += PostBufferTime;
            }
            
            have_peak=0;
            ts-=dSt;
            Vb=summedPulse(vfloor, vmax, inputTimes, inputVolts, numPulses, ts, t0);
            ts+=dSt;
            Va=0;
        } else {
            have_peak=0;
            ts += dSt;
            Vb = Va;
        }
        
    }
    
    return outK;
    
};



size_t ProcessPulsePileup( float *& inputTimes, float *& inputVolts, size_t numInputs, float *& outputTimes, float *& outputVolts, float window_width, float threshold, float vfloor, float vmax  )
{
    
    size_t counts, totCounts=0, PPUorder;
    size_t j,k;
    
    float * ptr_INtime, * ptr_OUTtime;
    float * ptr_INvolt, * ptr_OUTvolt;
    
    ptr_OUTtime = outputTimes;
    ptr_OUTvolt = outputVolts;
    
    k=0;
    while ( k < numInputs ) {
        
        j=k;
        while ( ((j+1) < numInputs) && ((inputTimes[j+1] - inputTimes[j]) <= window_width) ) j++;
        
        ptr_INtime = inputTimes+k;
        ptr_INvolt = inputVolts+k;
        PPUorder = j-k;
        
        counts = adcMeasurement( threshold, vfloor, vmax, ptr_INtime, ptr_INvolt, PPUorder+1, ptr_OUTtime, ptr_OUTvolt );
        
        //cout << PPUorder << " -> " << counts << endl; 
        
        ptr_OUTtime += counts;
        ptr_OUTvolt += counts;
        
        totCounts += counts;
        
        k+=(PPUorder+1);
        
        
    }
    
    return totCounts;
}

size_t ProcessPPUtoQuota( size_t detectionQuota, size_t& numNeeded, float *& inputTimes, float *& inputVolts, size_t numInputs, float *& outputTimes, float *& outputVolts, float window_width, float threshold, float vfloor, float vmax  )
{
    
    size_t counts, totCounts=0, PPUorder;
    size_t j,k;
    
    float * ptr_INtime, * ptr_OUTtime;
    float * ptr_INvolt, * ptr_OUTvolt;
    
    ptr_OUTtime = outputTimes;
    ptr_OUTvolt = outputVolts;
    
    numNeeded=0;
    k=0;
    while ( (k < numInputs) && (totCounts < detectionQuota) ) {
        
        j=k;
        while ( ((j+1) < numInputs) && ((inputTimes[j+1] - inputTimes[j]) <= window_width) ) j++;
        
        ptr_INtime = inputTimes+k;
        ptr_INvolt = inputVolts+k;
        PPUorder = j-k;
        
        counts = adcMeasurement( threshold, vfloor, vmax, ptr_INtime, ptr_INvolt, PPUorder+1, ptr_OUTtime, ptr_OUTvolt );
        
        //cout << PPUorder << " -> " << counts << endl; 
        
        if ( counts + totCounts > detectionQuota ) {
            counts = detectionQuota - totCounts;
        }
        
        
        
        ptr_OUTtime += counts;
        ptr_OUTvolt += counts;
        
        totCounts += counts;
        
        k+=(PPUorder+1);
        
        
    }
    
    numNeeded=k;
    
    
    return totCounts;
}


size_t FilterDeadtime( float *& inputTimes, float *& inputPH, size_t numInputs, float *& outputTimes, float *& outputPH, float deadTime )
{
    
    size_t counts=0, totCounts=0, PPUorder;
    size_t j,k;
    
    k=0;
    counts=0;
    while ( k < numInputs ) {
        
        j=k;
        while ( ((j+1) < numInputs) && (fabs(inputTimes[j+1] - inputTimes[k]) <= deadTime) ) j++;
        PPUorder = j-k;
        outputPH[counts] = inputPH[k];
        outputTimes[counts++] = inputTimes[k];
        
        totCounts += 1;
        
        k+=(PPUorder+1);
        
        
    }
    
    return totCounts;
}


