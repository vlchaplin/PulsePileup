//
//  PPUProb.cpp
//  C_ppu
//
//  Created by vchaplin on 10/15/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

/******
 
 Program to compute pileup likelihood functions Pr( E | E0, E1 ) up to arbitrary orders
 
*******/



#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <sstream>

#include "PulseShape.h"
#include "GapStatistics.h"
#include "MeasuredEnergies.h"
#include "PrRspsHDF5.h"

#include "EdgeSet.hh"

#include "spoccExeUtilities.h"

using namespace std;

void _print_banner() {
    
    cout << endl;
    cout << endl;
    cout << "$$$$$$$  MakeH5Response  $$$$$$$$$" << endl;
    cout << "$" <<endl;
    cout << "$ This program is used to compute the 'pileup response' file required for the PPU model code." << endl;
    cout << "$ Output is a single HDF5 file, based on user inputs via STDIN" << endl;
    cout << "$" <<endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    
    cout << endl;
    cout << endl;
    
}

void _print_usage() {
    
    cout << "Takes no arguments." << endl << endl;
    
}

int main (int argc, char * const argv[]) {

    H5std_string file;
    
    _print_banner();
    
    char ** argitr = getCmdOption( (char**)argv, (char**)argv+argc, "--help", 0);
	if ( argitr != NULL ) {
        _print_usage();
        return 0;
	}
    
    int minOrder = 1; // >=1. 
    int maxOrder = 7;
    
    //nIn & nOut must be the same at the moment (because of iterative approx)
    const size_t nInBins = 128;
	const size_t nOutBins= 128;
	size_t oChan, i, j;    
	clock_t begin, end;
    
    static float PrFunc_1[nInBins][nInBins][nOutBins];
	static float PrFunc_k[nOutBins][nInBins][nOutBins];
    static float TailPr_k[nOutBins][nInBins][nOutBins];
    static float PMixPr_k[nOutBins][nInBins][nOutBins];

	static float outBins[nOutBins+1];
	static float inBins[nInBins+1];
    
    size_t prDims_1[] = { nInBins, nInBins, nOutBins };
    size_t prDims_k[] = { nOutBins, nInBins, nOutBins };

    int k;
    //double ds = 0.01; noisy
    //double ds = 0.001; good
    double ds = 0.0002; //really really good (takes a while though; ~ 8 min per order, 2.2 Ghz i7)
    
    double deadStart = TauPEAK + ds;
    
    EdgeFunction<double> inEdgeFunc;
	EdgeFunction<double> outEdgeFunc;
    
    bool doVoltage;
    
    bool ok=false;
    string userInput;
    
    
    cout << "Select option for temporal step size in integration:" << endl;
    cout << setw(10) << "1" << ": 0.01 microseconds (results in noisy spectra, but fast... good for testing)" << endl;
    cout << setw(10) << "2" << ": 0.001 microseconds (good; some slight noise in spectra)" << endl;
    cout << setw(10) << "3" << ": 0.0002 microseconds (excellent; results are smooth for 128 channel data. Takes ~ 8 min per order, 2.2 GhZ i7)" << endl;
    cout << "Enter your choice or leave blank for default [Default is 3, which was used for the paper] :" << endl;
    
    int opt;
    while (! ok) {
        getline( cin, userInput );
        
        if ( userInput == "q" ) return 0;
        
        opt = atoi(userInput.c_str());
        
        ok=true;
        if ( opt==1 ) {
            ds = 0.01;
        }
        else if (opt == 2)
        {
            ds = 0.001;
        }
        else if (opt == 3)
        {
            ds = 0.0002;
        }
        else {
            cout << "Choose 1-3" << endl;
            ok = false;
        }
        
    }
    
    
    cout << endl;
    ok = false;
    while (! ok) {
        cout << "Compute in voltage units? (y/n, 'q' quits):" << endl;
        getline( cin, userInput );
        
        if ( userInput == "q" ) {
            return 0;
        }
        else if ( userInput == "y" ) {
            doVoltage=true;
            ok=true;
        }
        else if ( userInput == "n" ) {
            doVoltage=false;
            ok=true;
        }
    }
    
    ChannelEnergy chan2en;
    
    if ( doVoltage ) {
    
        //file.assign( "/gbm/TGF/PPU/C_ppu/PPUprob_Volts2.h5" );
        
        cout << "Voltage basis" << endl;
        
        chan2en.Set(0.0, 5.0, 5.0);
        
        inEdgeFunc.SetParameters(0.0045, 6.0, nInBins);
        inEdgeFunc.UseLog(1);
        outEdgeFunc.SetParameters(0.0045, 6.0, nOutBins);
        outEdgeFunc.UseLog(1);

    }
    else 
    {
        cout << endl;
        cout << "Enter the approximate LLT and overlow channel boundary in keV (e.g., for NAIs type: '4.0  1000.0'. For B0 type: '104.0 45000.0', for B1: '150.0 45000.0':" << endl;
        cout << endl;
        
        getline( cin, userInput );
        
        stringstream ss;
        
        ss << userInput;
        
        double llt_kev, thresh_kev;
        
        ss >> llt_kev;
        ss >> thresh_kev;
        
        
        
        cout << llt_kev << " , " << thresh_kev << endl;
        
        
        chan2en.Set( llt_kev, thresh_kev, 5.0);
        
        inEdgeFunc.SetParameters(0.95*llt_kev, 1.05*thresh_kev, nInBins);
        inEdgeFunc.UseLog(1);
        outEdgeFunc.SetParameters(0.95*llt_kev, 1.05*thresh_kev, nOutBins);
        outEdgeFunc.UseLog(1);
        
        
    }
    
    while ( file.length()==0 ) {

        cout << endl;
        cout << "Enter the name of the output .h5 file:" << endl;
        getline( cin, userInput );
        
        file = userInput;
    }

    
    cout << endl;
    ok=false;
    
    cout << "Maximum order to calculate? Must be >=1, <= 10. [default is:" << maxOrder << "]" << endl;
    cout << "(higher orders will still be available via an approximation, but order >10 here is beyond diminishing returns):" << endl;
    getline( cin, userInput );
    
    while ( !ok ) {
        
        if ( userInput == "q" ) return 0;
        
        if ( userInput.length() > 0 ) {
            maxOrder = atoi(userInput.c_str());
            if (maxOrder >= 1 && maxOrder <= 10) {
                ok=true;
            }
        } else {
            ok = true;
        }
        
        if ( !ok) {
            cout << "Maximum order to calculate? Must be >=1, <= 10. [default is:" << maxOrder << "]" << endl;
            getline( cin, userInput );
        }
    }
    
    cout << "Maximum order stored in the file will be: " << maxOrder << endl;
    
    //initialize Pr arrays to zero
    for (oChan=0; oChan<nOutBins; oChan++) {

        outBins[oChan] = outEdgeFunc[oChan];
		for (j=0; j<nInBins; j++) {
			for (i=0; i<nInBins; i++) {
				PrFunc_1[i][j][oChan] = 0.0;
			}
            
            for (i=0;i<nOutBins;i++) {
                PrFunc_k[i][j][oChan] = 0.0;
                TailPr_k[i][j][oChan] = 0.0;
                PMixPr_k[i][j][oChan] = 0.0;
            }
            
            if (oChan == 0) inBins[j] = inEdgeFunc[j];
		}
	}
    
    inBins[nInBins] = inEdgeFunc[nInBins];
	outBins[nOutBins] = outEdgeFunc[nOutBins];
    
    /************
     
     * FIRST ORDER
     
     ************/ 
    if ( minOrder == 1 ) {
     
        cout << endl << "Calculating first order..." << endl;
        begin = clock();
        k=1;
        
        CalcProbTensor3( (float*)PrFunc_1, prDims_1, inEdgeFunc, inEdgeFunc, outEdgeFunc, 0.0, TauPEAK, ds, 1, &PeakE, chan2en );
        CalcProbTensor3( (float*)PMixPr_k, prDims_k, outEdgeFunc, inEdgeFunc, outEdgeFunc, deadStart, TailStart+ds, ds, k, &TailE, chan2en );
        CalcProbTensor3( (float*)TailPr_k, prDims_1, inEdgeFunc, inEdgeFunc, outEdgeFunc, TailStart, TailEnd, ds, 1, &TailE, chan2en );
        
        end = clock();
        
        cout << " First order tensors took (sec): " << double( end - begin )/ CLOCKS_PER_SEC  << endl << endl;
        cout << "Writing " << file << endl;
        new_prsp_hdf5( file, (hsize_t)nInBins, (hsize_t)nInBins, (hsize_t)nOutBins, (float*)PrFunc_1, inBins, outBins );
        add_prsp_hdf5( file, (hsize_t)nOutBins, (hsize_t)nInBins, (hsize_t)nOutBins, (float*)TailPr_k, 1, TAIL_PPU_TENSOR );
        add_prsp_hdf5( file, (hsize_t)nOutBins, (hsize_t)nInBins, (hsize_t)nOutBins, (float*)PMixPr_k, 1, MIX_PPU_TENSOR );
        
        minOrder=2;
    }
    
    
    /************
     
     * HIGHER ORDERs
     
     ************/ 
    for (k=minOrder; k <= maxOrder; k++) 
    {
        cout << endl << "Calculating order " << k << endl;
        begin = clock();
        
        for (oChan=0; oChan<nOutBins; oChan++) {
            for (j=0; j<nInBins; j++) {
                for (i=0;i<nOutBins;i++) {
                    PrFunc_k[i][j][oChan] = 0.0;
                    TailPr_k[i][j][oChan] = 0.0;
                    PMixPr_k[i][j][oChan] = 0.0;
                }
            }
        }
        
        CalcProbTensor3( (float*)PrFunc_k, prDims_k, outEdgeFunc, inEdgeFunc, outEdgeFunc, 0.0, TauPEAK, ds, k, &PeakE, chan2en );
        CalcProbTensor3( (float*)PMixPr_k, prDims_k, outEdgeFunc, inEdgeFunc, outEdgeFunc, deadStart, TailStart+ds, ds, k, &TailE, chan2en );
        CalcProbTensor3( (float*)TailPr_k, prDims_k, outEdgeFunc, inEdgeFunc, outEdgeFunc, TailStart, TailEnd, ds, k, &TailE, chan2en );
        
        end = clock();
        cout << "Order " << k << " tensors took (sec): " << double( end - begin )/ CLOCKS_PER_SEC  << endl << endl;
        cout << "Writing " << file << endl;
        add_prsp_hdf5( file, (hsize_t)nOutBins, (hsize_t)nInBins, (hsize_t)nOutBins, (float*)PrFunc_k, k );
        add_prsp_hdf5( file, (hsize_t)nOutBins, (hsize_t)nInBins, (hsize_t)nOutBins, (float*)TailPr_k, k, TAIL_PPU_TENSOR );
        add_prsp_hdf5( file, (hsize_t)nOutBins, (hsize_t)nInBins, (hsize_t)nOutBins, (float*)PMixPr_k, k, MIX_PPU_TENSOR );
    }
    
    //Finally, write the max order calculated into the file
    const IntType type( PredType::ALPHA_I32 );
	DataSpace scalars;
    
    H5File h5obj(file, H5F_ACC_RDWR );
    DataSet miscDset;
    Attribute attr;
    
    // This is really stupid. In the HDF5 library the only way to check for a data set's existence is 
    // to catch the IO exception error after trying to open it.
    // Trying to create a data set that already exists causes a run-time error.
    //
    //  Lots of user have commented on the idiocy of this method, so it will change at some point...
    try {  
        cout << "Opening existing '/misc' data set" << endl;
        miscDset = h5obj.openDataSet("misc");
    }
    catch( FileIException not_found_error )
    {
        cout << "Creating '/misc' data set" << endl;
        miscDset = h5obj.createDataSet("misc", type, scalars);
    }
    
    
    //Fortunately there is a method for checking attribute existence before creating
    if ( H5Aexists(miscDset.getId(), "MaxOrder") ){
        cout << "Updating 'MaxOrder' attribute = " << maxOrder << endl;
        attr = miscDset.openAttribute("MaxOrder");
    } else {
        cout << "Creating 'MaxOrder' attribute = " << maxOrder << endl;
        attr = miscDset.createAttribute("MaxOrder", type, scalars);
    }

    
    attr.write(PredType::NATIVE_INT, &maxOrder);
    
    
    h5obj.close();
    
    
    
    cout << "Done" << endl;
    
    return 0;

}