//
//  PPUExample.cpp
//  C_ppu
//
//  Created by vchaplin on 7/8/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "ppu_static_interface.h"


int main (int argc, char * const argv[]) {

    // The PPU model requires the following input for each detector:  count rate, channel edges, counts spectrum
    //Begin by assuming these have been defined for two detectors, (n9 and n11).
    
    double n9_counts_per_microsecond = 0.15; //150,000 cps
    double n11_counts_per_microsecond= 0.05; // 50,000 cps
    
    double * n9_channel_edges_keV;
    double * n11_channel_edges_keV;
    
    double * n9counts;
    double * n11counts;
    
    long n9_nChan = 128;
    long n11_nChan = 128;
    int i;
    
    //   Code to define the edge arrays and counts spectra...(i.e., folding photon model through DRMS)
    //   ...
    //   ...
    
    
    
    // --------------------  PPU set-up ------------------------------- ||
    
    //Specify the PPU file.  THis is a specially formatted HDF5-type 
    //file containing information on the "pileup response" of GBM detectors.
    
    // This file is for NAIs, calculated assuming 4.0 keV -> 0 volts and 1000.0 keV -> 5.0 volts
    char NAI_ppu_file[] = "/gbm/TGF/PPU/C_ppu/PPUprob_NAI_better4.h5";
    
    
    
    
    //Prior to first read, use this routine to indicate the maximum order one might wish to calculate.  This is necessary for
    //memory to be allocated in the right amount when the file is loaded. Orders up to 8 have been tested.
    //The default is 8, but this will be overkill for some sources. This max order can't be increased after ppu_load_h5_file() is called.
    
    //See the comments in ppu_static_interface.h for a guide to approximation orders. 
    int maxOrder = 6;
    ppu_set_maximum_order(maxOrder);
    
    
    //This routines opens the HDF5 file and loads it's specially formatted contents into memory (and then closes the file).
    //If the file has already been loaded, this is a no-op.  Enough memory is allocated to accomadate the max approx. order, 
    //which was set by ppu_set_approximation_order().
    ppu_load_h5_file(NAI_ppu_file);
    
    
    
    
    
    // Now "bind" a model instance to this data file, via a string tag (any alpha-numeric string will work).
    //  
    // Binding allows multiple detectors to share certain pieces of data, improving the
    // code's performance
    
    //In this case the tag is "n9"
    char tag1[] = "n9";
    ppu_create_instance_from_file(tag1, NAI_ppu_file);
    
    //In this case the tag is "n11"
    char tag2[] = "n11";
    ppu_create_instance_from_file(tag2, NAI_ppu_file);
    

    
    //Now specify what approximation order to calculate for each detector.  They must be <= the maximum approx order.
    //Using higher orders for brighter detectors results in more accuracy, but takes longer to compute.
    
    ppu_set_approximation_order(tag1, 3);
    ppu_set_approximation_order(tag2, 1);
    
    //The each data set must be re-initialized after the order is set.  These orders can be changed, if for example the true rate might be higher than
    //estimated at the start of the program.  
    
    // --------------------  End of set-up ------------------------------- ||
    
    
    // COMPUTING THE MODEL
    
    //If the model has not been initialized, if the input spectral shape is changed, or if the approximation order is changed, ppu_initialize_* must be called.
    
    //The following routine initializes internal data and copies the input counts spectrum.  The input is counts per channel not counts per kev.
    //The input count spectrum is re-normalized to 1.0
    
    //do for n9 (tag1) (slow)
    ppu_initialize_zeroth_order_counts(tag1, n9_channel_edges_keV, n9counts, n9_nChan);
    
    //do for n11 (tag2...faster since approx order is smaller)
    ppu_initialize_zeroth_order_counts(tag2, n11_channel_edges_keV, n11counts, n11_nChan);
    
    
    
    //Final step
    
    //arrays to hold the computed spectra
    double n9_output_counts[n9_nChan];
    double n11_output_counts[n11_nChan];
    
    //Once initialized the model can be computed (fast)
    ppu_compute_distorted_spectrum(tag1, n9_counts_per_microsecond, n9_channel_edges_keV, n9_nChan, n9_output_counts);
    
    //the output array is actually a probability spectrum.  Multiply by the number of input counts to get the count spectrum. Normalization is <1.0 due to deadtime & PPU losses
    
    //Rate * bin size (for example a 256 ms bin)
    double n9_expected_counts = (n9_counts_per_microsecond * 1e06) * 0.256;
    
    for (i=0; i<n9_nChan; i++) {
        n9_output_counts[i] *= n9_expected_counts;
    }
    
    
    //do the same for n11 and other detectors...
    ppu_compute_distorted_spectrum(tag2, n11_counts_per_microsecond, n11_channel_edges_keV, n11_nChan, n11_output_counts);
    
    double n11_expected_counts = (n11_counts_per_microsecond * 1e06) * 0.256;

    for (i=0; i<n11_nChan; i++) {
        n11_output_counts[i] *= n11_expected_counts;
    }
    
    
    
    
    //Now n9_output_counts & n11_output_counts contain the count model with PPU in each detectro
    
    
    
    return 0;

}