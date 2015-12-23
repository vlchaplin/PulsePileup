/*
 *  MeasuredEnergies.cpp
 *  C_ppu
 *
 *  Created by vchaplin on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#include "MeasuredEnergies.h"

//using namespace paramSet;

ChannelEnergy chan2en_default( 150.0, 45000.0 );

double PeakE( double e0, double e1, double s1, ChannelEnergy& chan2en ) 
{
	if ( s1 < 0 ) return 0;
    
    double f0 = chan2en.EtoV(e0);
	double f1 = chan2en.EtoV(e1);
	double peakV;
	
	if ( s1 <= TauPEAK ) {
		
		double t1 = t1p( f0, f1, s1 );
        
        if ( t1 < t0p*exp(1.7*s1*s1) ) t1 = t0p;
		
		peakV = f0*pulse(t1) + f1*pulse(t1 - s1);
		
	} else {
		peakV = f0*pulse(t0p);
	}
	
	return chan2en.VtoE(peakV);
};

double TailE( double e0, double e1, double s1, ChannelEnergy& chan2en ) 
{
	if ( s1 < 0 ) return 0;
    
    double peakV;
	
    double f0 = chan2en.EtoV(e0);
    double f1 = chan2en.EtoV(e1); 
    double t1 = t1p( f0, f1, s1 );
    double v0,v1;
    
    const double BreakPos = 1.0;
    const double BreakScale = 30.0;
    
    double logisticJoin = 1.0 / (1.0 + exp(-BreakScale*(s1 - BreakPos)));
    
    t1 = t1*( 1.0 - logisticJoin ) + (t0p + s1)*logisticJoin;
    
    v0 = f0*pulse(t1);
    v1 = f1*pulse(t1 - s1);
    //if (v0 < Vfloor) v0 = Vfloor;   
    peakV = v0 + v1;
    
    
	
	if (peakV < chan2en.vthresh ) return 0;
		
	return chan2en.VtoE(peakV);
};


double MixE( double e0, double e1, double e2, double s2, ChannelEnergy& chan2en ) 
{
	if ( s2 < 0 ) return 0;
    
    double v0 = chan2en.EtoV(e0);
	double v1 = chan2en.EtoV(e1);
    double v2 = chan2en.EtoV(e2);
	double peakV;
    
    peakV = (v0+v1)*pulse( s2 + t0p );
	
	//if (peakV < Vfloor) peakV = Vfloor;
	
	peakV = peakV+v2;
	if (peakV < chan2en.vthresh ) return 0;
    
	return chan2en.VtoE(peakV);
};


double dPeakEdS( double e0, double e1, double s1, ChannelEnergy& chan2en )
{
	if ( s1 < 0 ) return 0;
    
    double v0 = chan2en.EtoV(e0);
	double v1 = chan2en.EtoV(e1);
	double t1 = t1p( v0, v1, s1 );
	
	double dVdT = t1p_dt( v0, v1, s1 )*( v0*pulse_dt(t1) + v1*pulse_dt(t1 - s1) ) - v1*pulse_dt(t1 - s1);
	return chan2en.VtoE(dVdT);
};

double dTailEdS( double e0, double e1, double s1, ChannelEnergy& chan2en )
{
	if ( s1 < 0 ) return 0;
    
    double v0 = chan2en.EtoV(e0);
	//double v1 = chan2en.EtoV(e1);
	double dVdT = v0*pulse_dt(s1 + t0p);// + 0*pulse_dt(t0p) ;
	
	return chan2en.VtoE(dVdT);
};

double PDFPeakE( double E0, double E1, double s1, int order, double interval, ChannelEnergy& chan2en )
{

	if ( s1 <= 0.0 ) return 0;
    
    //eR = PeakE( E0, E1, s1, chan2en );
	double dEdS = fabs(dPeakEdS(E0, E1, s1, chan2en));
	double dPdE;
	
	if ( dEdS < 1e-03 ) return 0;
	if ( s1 <= interval ) 
		dPdE = dPdS(s1, order, interval) / dEdS;
	else
		dPdE = 0.0;
	
	return dPdE;
};

double PDFTailE( double E0, double E1, double s1, int order, double interval, ChannelEnergy& chan2en )
{
    
	//eR = PeakE( E0, E1, s1, chan2en );
	double dEdS = fabs(dTailEdS(E0, E1, s1, chan2en));
	double dPdE;
	
	if ( dEdS < 1e-03 ) return 0;
    s1 -= PARAM_SET::TailStart;
	dPdE = dPdS(s1, order, interval) / dEdS;
	
	return dPdE;
};










