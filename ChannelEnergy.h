/*
 *  ChannelEnergy.h
 *  C_ppu
 *
 *  Created by vchaplin on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ChannelEnergy_H
#define ChannelEnergy_H

#include <iostream>

//const double max_V = 5.0;
const double Vfloor = -0.5;

class ChannelEnergy {

private:
	double max_kev;
	double llt_kev;
public:
	
    double max_V;
	double adc_llt;
	double vthresh;
	
    ChannelEnergy(){};
    
	ChannelEnergy(double minE, double maxE, double MaxVoltage=5.0){
		
        max_V = MaxVoltage;
        
		this->max_kev = maxE;
		
		adc_llt = energy2adc(minE);
		
		vthresh = EtoV(minE);
		
	};
	
    virtual void Set(double minE, double maxE, double MaxVoltage=5.0)
    {
        max_V = MaxVoltage;
        
		this->max_kev = maxE;
		
		adc_llt = energy2adc(minE);
		
		vthresh = EtoV(minE);
    };
    
	virtual double adc2energy( double adc ) {
		if( adc < 4095 ) 
			return adc* max_kev/4095.0;
		else 
			return max_kev;
	};
	virtual double energy2adc( double keV ) {
		if ( keV < max_kev ) 
			return (keV/max_kev)*4095;
		else 
			return 4095.0;
	};
	
	double volts2adc( double volts ) { return volts*4095.0/max_V; };
	double adc2volts( double chan ) { return chan*max_V/4095.0; };
	
	
	double EtoV( double keV ) {
		return adc2volts(energy2adc(keV));
	};
	double VtoE( double volts ) {
		return adc2energy(volts2adc(volts));
	};
};

#endif
	
	
