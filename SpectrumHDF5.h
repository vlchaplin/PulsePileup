//
//  SpectrumHDF5.h
//  C_ppu
//
//  Created by vchaplin on 7/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef C_ppu_SpectrumHDF5_h
#define C_ppu_SpectrumHDF5_h



#include "H5Cpp.h"
#include <iostream>

//#ifndef H5_NO_NAMESPACE
using namespace H5;
//#endif

int write_spectrum_hdf5( H5File& h5file, hsize_t nbins, float * data, float * edges );

int read_spectrum_hdf5( H5File& h5file, hsize_t& nbins, float *& data, float *& edges, float& fracRec, long& thrown );

int read_spectrum_hdf5( std::string h5file, hsize_t& nbins, float *& data, float *& edges, float& fracRec, long& thrown);

#endif
