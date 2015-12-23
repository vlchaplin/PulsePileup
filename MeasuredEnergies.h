/*
 *  MeasuredEnergies.h
 *  C_ppu
 *
 *  Created by vchaplin on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MEASUREDENERGY_H
#define MEASUREDENERGY_H

#include <math.h>
#include <iostream>
#include "PulseShape.h"
#include "GapStatistics.h"
#include "EdgeSet.hh"

using namespace std;
using namespace PARAM_SET;

extern ChannelEnergy chan2en_default;

double PeakE( double E0, double E1, double s1, ChannelEnergy& chan2en = chan2en_default );
double TailE( double E0, double E1, double s1, ChannelEnergy& chan2en = chan2en_default );
double MixE( double e0, double e1, double e2, double s2, ChannelEnergy& chan2en = chan2en_default );
double dPeakEdS( double E0, double E1, double s1, ChannelEnergy& chan2en = chan2en_default );
double dTailEdS( double E0, double E1, double s1, ChannelEnergy& chan2en = chan2en_default );

double PDFPeakE( double E0, double E1, double s1, int order, double interval=PARAM_SET::TauPEAK, ChannelEnergy& chan2en = chan2en_default  );
double PDFTailE( double E0, double E1, double s1, int order, double interval=PARAM_SET::TauTAIL, ChannelEnergy& chan2en = chan2en_default  );

//, double(*mixingFunc)( double, double, double, ChannelEnergy&), ChannelEnergy& chan2en = chan2en_default
//int CalcProbTensor3( void * PrTensor, size_t * dims, EdgeSet<double>& basis0, EdgeSet<double>& basis1, EdgeSet<double>& basis2, double t0, double tEnd, double dt, int order  );


inline int CalcProbTensor3( float * PrTensor, size_t * dims, EdgeSet<double>& basis0, EdgeSet<double>& basis1, EdgeSet<double>& basis2, 
                            double t0, double tEnd, double dt, int order,
                            double(*mixingFunc)( double, double, double, ChannelEnergy&), ChannelEnergy& chan2en = chan2en_default)
{
    double s_ij,s_ij_next, e_i, e_j, eR, integralStep;
    size_t i,j,oChan,imax,jmax,omax,long0,long1;
    
    double interval = tEnd - t0;
    
    imax = basis0.size();
    jmax = basis1.size();
    omax = basis2.size();
    
    for (i=0; i<imax; i++) {
		
		e_i = basis0.binMid(i);
        
        long0 = i * jmax * omax;
        
		for (j=0; j<jmax; j++) {
            
			e_j = basis1.binMid(j);
            
            long1 = j * omax;
            
            for(s_ij=t0; s_ij <= tEnd ; s_ij += dt) {
				
                s_ij_next = s_ij+dt;
                
				eR = mixingFunc(e_i, e_j, s_ij, chan2en);
                
                if ( eR >= basis2.binLower(0) ) {
                    
                    oChan = (size_t) basis2.lookup(eR);
                    if ( oChan >= omax ) oChan = omax-1;
                    
                    integralStep = Pr_s( s_ij - t0, s_ij_next - t0, interval, order );
                    
                    /*Next line is equivalent to:
                     
                     PrTensor[i][j][oChan] += integralStep
                     
                     ...3-d array subscripts can't be used in this context
                     */
                    
                    *(PrTensor + long0 + long1 + oChan) += integralStep;
                }
			}
        }
        
    }
    
    return 0;
}



inline int CalcDegenProbTensor3( float * PrTensor, size_t * dims, EdgeSet<double>& basis0, EdgeSet<double>& basis1, EdgeSet<double>& basis2, 
                        double t0, double tEnd, double dt, double normalizationInt, int order,
                        double(*mixingFunc)( double, double, double, ChannelEnergy&),
                        double(*mixingDfDt)( double, double, double, ChannelEnergy&),
                        ChannelEnergy& chan2en = chan2en_default)
{
    double s_ij,s_ij_next, e_i, e_j, eR,dEds,dEdsb;
    size_t i,j,oChan,imax,jmax,omax,long0,long1;
    
    //double interval = tEnd - t0;
    float functionalStep;
    
    imax = basis0.size();
    jmax = basis1.size();
    omax = basis2.size();
    
    float * energyPDF = new float[omax];
    float timePDFsampA, timePDFsampB;
    
    for(oChan=0; oChan < omax; oChan++) {
        energyPDF[oChan]=0;
    }
    
    for (i=0; i<imax; i++) {
		
		e_i = basis0.binMid(i);
        
        long0 = i * jmax * omax;
        
		for (j=0; j<jmax; j++) {
            
			e_j = basis1.binMid(j);
            
            long1 = j * omax;
            
            //In the first step we define the energy pdf. It is implicitly defined,
            //and may have more than one value so the '+=' in this loop addresses that (i.e., oChan may be visited more than once)
            for(s_ij=t0; s_ij <= tEnd ; s_ij += dt) {
				
                s_ij_next = s_ij+dt;
                
				eR = mixingFunc(e_i, e_j, s_ij, chan2en);
                
                if ( eR >= basis2.binLower(0) ) {
                    
                    oChan = (size_t) basis2.lookup(eR);
                    if ( oChan >= omax ) oChan = omax-1;
                    
                    dEds  = fabs(mixingDfDt(e_i, e_j, s_ij, chan2en));
                    dEdsb = fabs(mixingDfDt(e_i, e_j, s_ij_next, chan2en));
                    
                    timePDFsampA = dPdS(s_ij - t0, order, normalizationInt);
                    timePDFsampB = dPdS(s_ij_next - t0, order, normalizationInt);
                    
                    if ( (timePDFsampA/dEds) > 0.1) {
                        functionalStep=0.1;
                    } else 
                        functionalStep = 0.5*(  dPdS(s_ij - t0, order, normalizationInt)/dEds + dPdS(s_ij_next - t0, order, normalizationInt)/dEdsb   );
                    
                    //functionalStep = 0.5*(  dPdS(s_ij - t0, order, normalizationInt)/dEds + dPdS(s_ij_next - t0, order, normalizationInt)/dEdsb   );
                    
                    energyPDF[oChan] += functionalStep;
                }
			}
            
            //Next integrate the pdf and store
            for(oChan=0; oChan < omax; oChan++) {
                *(PrTensor + long0 + long1 + oChan) += (energyPDF[oChan] * basis2.binSize(oChan));
                energyPDF[oChan]=0;
            }
        }
        
    }
    
    delete [] energyPDF;
    
    return 0;
}





#endif