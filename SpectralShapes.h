//
//  SpectralShapes.h
//  C_ppu
//
//  Created by vchaplin on 1/23/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef C_ppu_SpectralShapes_h
#define C_ppu_SpectralShapes_h

#include <math.h>
#include "EdgeSet.hh"

const float PI = 4.0*atan(1.0);

inline float boxCarPdf(float energy, float& mean, float& width )
{
    float emin = mean-width/2.0;
    float emax = mean+width/2.0;
    
    if ( energy < emin || energy > emax ) return 0;
    
    return 1.0/width;
}

inline float power_law_pdf(float energy, float& norm, float& plIndex, float& E0)
{
    return powf( (energy/E0), plIndex)*norm;
}

inline float power_law_norm(float plExp, float E0, float emin, float emax)
{
    if ( plExp == -1.0 ) 
        return 1.0 / ( E0 * ( log(emax) - log(emin) )  );
    else
        return (1.0 + plExp) / (emax * powf( emax / E0, plExp ) - emin * powf(emin / E0, plExp ));
}

inline float cutoffpl_emax(float index, float e0, float ecut, float emin, float emax)
{
    if (index > 0) return index*ecut;
    else {
        float norm=1;
        return power_law_pdf(emin,norm,index,e0);
    }
};

inline float gauss_pdf(float energy, float& mean, float& sig)
{
    return (1.0/(sig*sqrt(2.0*PI)))*exp(-powf(energy - mean, 2.0) / (2.0*sig*sig) );
}

inline float hardness_spectrum( float energy, EdgeFunction<float>& inputEnergy, float r1, float r2, float r3)
{
    static float bands[] = { 150.0, 1000.0, 10000.0, 45000.0 };
    static float center1 = (bands[0] + bands[1])/2.0;
    static float center2 = (bands[1] + bands[2])/2.0;
    static float center3 = (bands[2] + bands[3])/2.0;
    
    static float width1 = bands[1] - bands[0];
    static float width2 = bands[2] - bands[1];
    static float width3 = bands[3] - bands[2];
    
    
    return width1*0.995*(r1*boxCarPdf(energy, center1, width1) + r2*boxCarPdf(energy, center2, width2) + r3*boxCarPdf(energy, center3, width3)) / (r1+r2+r3);
    
}

inline float hardness_spectrum( float energy, float * bandEdges, float * weights, int nbands, bool updateMax=false)
{
    int i=0;
    float funcValue=0;
    float totWeight=0;
    float center, width;
    
    static bool HaveMax=false;
    static float widthOfMaxBand=1.0;
    
    for(i=0; i<nbands; i++) 
    {
        width = bandEdges[i+1] - bandEdges[i];
        center = bandEdges[i] + width / 2.0;
        
        totWeight += weights[i];
        
        funcValue += weights[i]*boxCarPdf(energy, center, width);
    }
    
    if ( updateMax) {
        
        widthOfMaxBand=0;
        
        for(i=0; i<nbands; i++) 
        {
            width = bandEdges[i+1] - bandEdges[i];
            if ( width > widthOfMaxBand ) widthOfMaxBand = width;
        }
        
        widthOfMaxBand *= 0.995;
    }
    
    
    return widthOfMaxBand*funcValue;
}



#endif
