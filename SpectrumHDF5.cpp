//
//  SpectrumHDF5.cpp
//  C_ppu
//
//  Created by vchaplin on 7/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "SpectrumHDF5.h"
#include <iostream>


int write_spectrum_hdf5( H5File& h5file, hsize_t nbins, float * data, float * edges )
{
    try
    //if(1)
	{
		Exception::dontPrint();
		
		DataSpace spectrum( 1, &nbins );
		DataSpace ebounds( 1, &(++nbins) );
		
		const FloatType type( PredType::IEEE_F32LE );
		
		DataSet spectrumDset = h5file.createDataSet( (char*)"SPECTRUM", type, spectrum  );
		DataSet einDataSet = h5file.createDataSet( (char*)"EBOUNDS", type, ebounds  );
		
        
		spectrumDset.write( data, PredType::NATIVE_FLOAT );
		einDataSet.write( edges, PredType::NATIVE_FLOAT );
		
		//attach_pulse_attr(prbDataSet);
		
	//	file.close();
	}
    
	
	catch( FileIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSet operations
	catch( DataSetIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error )
	{
		error.printError();
		return -1;
	}
	
    
	return 0;
}

int read_spectrum_hdf5( std::string filename, hsize_t& nbins, float *& data, float *& edges, float& fracRec, long& thrown )
{
    H5File file( filename, H5F_ACC_RDONLY );
    
    int rc = read_spectrum_hdf5( file, nbins, data, edges, fracRec, thrown );
    
    file.close();
    
    return rc;
}


int read_spectrum_hdf5( H5File& h5file, hsize_t& nbins, float *& data, float *& edges, float& fracRec, long& thrown )
{
    try
    //if(1)
	{
		
        Exception::dontPrint();
        
        DataSet spectrumDset = h5file.openDataSet( (char*)"SPECTRUM" );
		DataSet einDataSet = h5file.openDataSet( (char*)"EBOUNDS" );
        DataSet inputs = h5file.openDataSet( (char*)"inputs" );
        
        DataSpace specSpace = spectrumDset.getSpace();
        DataSpace edgeSpace = einDataSet.getSpace();
        
        specSpace.getSimpleExtentDims(&nbins);
        
        hsize_t nedges = nbins+1;
        
        DataSpace memspace1(1, &nbins);
        DataSpace memspace2(1, &nedges);
        
        data = new float[nbins];
        edges = new float[nedges];
        
        spectrumDset.read( data, PredType::NATIVE_FLOAT, memspace1, specSpace );
        einDataSet.read(  edges, PredType::NATIVE_FLOAT, memspace2, edgeSpace );
        
        Attribute attA = inputs.openAttribute((const char*)"NumThrown");
        Attribute attN = inputs.openAttribute((const char*)"NumRecorded");
        
        long numCounted;
        
        attN.read( PredType::NATIVE_LONG, &(numCounted) );
        attA.read( PredType::NATIVE_LONG, &thrown);
        
        
        fracRec = (float)numCounted/(float)thrown;
        
    }
    
	
	catch( FileIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSet operations
	catch( DataSetIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error )
	{
		error.printError();
		return -1;
	}
	
    
	return 0;
}




