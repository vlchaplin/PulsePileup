/*
 *  GapStatistics.h
 *  C_ppu
 *
 *  Created by vchaplin on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GapStatistics_H
#define GapStatistics_H

#include <math.h>
#include <stdlib.h>
#include "PulseShape.h"


inline double dPdS_peak( double& s1, int& k ){
    if ( s1 > PARAM_SET::TauPEAK || s1 < 0 ) 
        return 0;
    else
        return k*pow( PARAM_SET::TauPEAK - s1, k-1 ) / pow( PARAM_SET::TauPEAK, k );	
};

inline double dPdS_tail( double& s1, int& k ){
	return k*pow( PARAM_SET::TauTAIL - s1, k-1 ) / pow( PARAM_SET::TauTAIL, k );	
};

inline double dPdS( double s1, int& k, double interval ){
    
    if ( s1 > interval || s1 < 0 ) return 0;
    
	return k*pow( interval - s1, k-1 ) / pow( interval, k );	
};
inline double Pr_s( double a, double b, double normalizeInterval, int& k ){
    
    double integrationInterval = b-a;
    if ( integrationInterval > normalizeInterval || integrationInterval < 0 ) return 0;
    
	return ( pow( normalizeInterval - a, k ) - pow( normalizeInterval - b, k ) ) / pow( normalizeInterval, k );	
};


inline double s1_peak_average( int order )
{
    if (order <= 0) return 0;
    else return PARAM_SET::TauPEAK / ( 1.0 + order );
};

/*
inline double dPdS_dead( double s1, int order ){
	return k*pow( PARAM_SET::TauDEAD - s1, k-1 ) / pow( PARAM_SET::TauDEAD, k );	
};*/



#endif