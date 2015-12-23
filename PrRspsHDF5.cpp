/*
 *  PrRspsHDF5.cpp
 *  C_ppu
 *
 *  Created by vchaplin on 3/29/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "PrRspsHDF5.h"

int attach_pulse_attr( DataSet& dset )
{
	const FloatType type( PredType::IEEE_F32LE );
	DataSpace scalars;
	
	Attribute att = dset.createAttribute( (const char*)"PeakEnd", type, scalars );
	att.write( PredType::NATIVE_DOUBLE, &(PARAM_SET::TauPEAK) );
	
	Attribute att2 = dset.createAttribute( (const char*)"TailStart", type, scalars );
	att2.write( PredType::NATIVE_DOUBLE, &(PARAM_SET::TailStart) );
	
	Attribute att3 = dset.createAttribute( (const char*)"TailEnd", type, scalars );
	att3.write( PredType::NATIVE_DOUBLE, &(PARAM_SET::TailEnd) );
	
	return 0;
};

/*****************************
 
 // n1 x n2 x n3  dimensions.   n3 = # bins in output spectra
 // 1st dim -> 'primary' count
 // 2nd dim -> '1st order' count 
 // 3rd dim -> 'recorded count' 
 
 ******************************/

int new_prsp_hdf5( H5std_string filename, hsize_t n1, hsize_t n2, hsize_t n3, float * data, float * edge1, float *edge3, int type  )
{
	
    H5std_string probTableName;
    H5std_string energyInputsName="eIN";
    H5std_string energyOutputsName="eOUT";
    
    if (type==PEAK_PPU_TENSOR) probTableName="Pr1";
    else if (type==TAIL_PPU_TENSOR) probTableName="TailPr1";
    else if (type==MIX_PPU_TENSOR) probTableName="MixPr1";
    else {
        return -1;
    }
    
	try
	{
		Exception::dontPrint();
		
		H5File file( filename, H5F_ACC_TRUNC );
		
		hsize_t dims[3];
		dims[0] = n1;
		dims[1] = n2;
		dims[2] = n3;
		
		DataSpace probTableSpace( 3, dims );
		DataSpace energyInputs( 1, &(++n1) );
		DataSpace energyOutputs( 1, &(++n3) );
		
		const FloatType type( PredType::IEEE_F32LE );
		
		DataSet prbDataSet = file.createDataSet( probTableName, type, probTableSpace  );
		DataSet einDataSet = file.createDataSet( energyInputsName, type, energyInputs  );
		DataSet eotDataSet = file.createDataSet( energyOutputsName, type, energyOutputs  );
		
		einDataSet.write( edge1, PredType::NATIVE_FLOAT );
		eotDataSet.write( edge3, PredType::NATIVE_FLOAT );
		
		prbDataSet.write( data, PredType::NATIVE_FLOAT );
		
		attach_pulse_attr(prbDataSet);
		
		file.close();
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

/*****************************
 
 // n1 x n2 x n3  dimensions.   nR -> # bins in output spectra
 // 1st dim -> 'primary' count
 // 2nd dim -> '1st order' count 
 // 3rd dim -> 'recorded count' 
 
 ******************************/

int add_prsp_hdf5( H5std_string filename, hsize_t n1, hsize_t n2, hsize_t n3, float * data, int order, int type )
{
    
    char probTableName[10];
    
    if (type==PEAK_PPU_TENSOR) sprintf(probTableName, "Pr%d", order);
    else if (type==TAIL_PPU_TENSOR) sprintf(probTableName, "TailPr%d", order);
    else if (type==MIX_PPU_TENSOR) sprintf(probTableName, "MixPr%d", order);
    else {
        return -1;
    }
    
    
	try
	{
		Exception::dontPrint();
		
		H5File file;
		file.openFile(filename, H5F_ACC_RDWR);
		
		hsize_t dims[3];
		dims[0] = n1;
		dims[1] = n2;
		dims[2] = n3;
				
		DataSpace probTableSpace( 3, dims );
		
		const FloatType type( PredType::IEEE_F32LE );
		
		DataSet prbDataSet = file.createDataSet( probTableName, type, probTableSpace  );
		
		prbDataSet.write( data, PredType::NATIVE_FLOAT );
		
		attach_pulse_attr(prbDataSet);
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
//get_prsp_dims_hdf5( H5std_string filename, int * ndims, hsize_t * dims, int order, int type )
int get_prsp_dims_hdf5( H5std_string filename, int * ndims, hsize_t *& dims, int order, int type )
{
    char probTableName[100];
    
    if (type==PEAK_PPU_TENSOR) sprintf(probTableName, "Pr%d", order);
    else if (type==TAIL_PPU_TENSOR) sprintf(probTableName, "TailPr%d", order);
    else if (type==MIX_PPU_TENSOR) sprintf(probTableName, "MixPr%d", order);
    else {
        return -1;
    }
    
    try
	{
		Exception::dontPrint();
        
        H5File file( filename, H5F_ACC_RDONLY );
        DataSet dataset = file.openDataSet( probTableName );
        DataSpace dataspace = dataset.getSpace();

        *ndims = dataspace.getSimpleExtentNdims();
        
        dims = new hsize_t[*ndims];
        
        dataspace.getSimpleExtentDims( dims, NULL);
        
        file.close();
        
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

/******************************
 
dims - already defined (array of dimensions)
data - already allocated for size dims
 
 // n0 x n1 x nR  dimensions.   nR -> # bins in output spectra
 // 1st dim -> 'primary' count
 // 2nd dim -> '1st order' count 
 // 3rd dim -> 'recorded count' 
 
******************************/
int read_prsp_hdf5( H5std_string filename, hsize_t * dims, float * data, int order, int type, float * eIN, float * eOUT )
{
    char probTableName[100];
    
    if (type==PEAK_PPU_TENSOR) sprintf(probTableName, "Pr%d", order);
    else if (type==TAIL_PPU_TENSOR) sprintf(probTableName, "TailPr%d", order);
    else if (type==MIX_PPU_TENSOR) sprintf(probTableName, "MixPr%d", order);
    else {
        return -1;
    }
    
    try
	{
		Exception::dontPrint();
    
        H5File file( filename, H5F_ACC_RDONLY );
        DataSet dataset = file.openDataSet( probTableName );
        DataSpace dataspace = dataset.getSpace();
        DataSpace mspace1(3, dims);
        
        //int ndims = dataspace.getSimpleExtentDims( dims, NULL);
        
        //int rank = dataspace.getSimpleExtentNdims();
        
        
        dataset.read( data, PredType::NATIVE_FLOAT, mspace1, dataspace );
        
        
        read_edges(file, eIN, eOUT);
        
        
        
        file.close();
        
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

int read_edges(H5File& fileObj, float * eIN, float * eOUT, bool alloc)
{
    try
	{
		Exception::dontPrint();
    
        if ( eIN != NULL ) {
            
            H5std_string energyInputsName="eIN";
            DataSet eset = fileObj.openDataSet( energyInputsName );
            DataSpace espace = eset.getSpace();
            hsize_t nedges;
            
            espace.getSimpleExtentDims(&nedges);
            
            DataSpace mspace(1, &nedges);
            
            if ( alloc ) eIN = new float[nedges];
            
            eset.read( eIN, PredType::NATIVE_FLOAT, mspace, espace );
            
        }
        
        if ( eOUT != NULL ) {
            
            H5std_string energyInputsName="eOUT";
            DataSet eset = fileObj.openDataSet( energyInputsName );
            DataSpace espace = eset.getSpace();
            hsize_t nedges;
            
            espace.getSimpleExtentDims(&nedges);
            
            DataSpace mspace(1, &nedges);
            
            if ( alloc ) eOUT = new float[nedges];
            
            eset.read( eOUT, PredType::NATIVE_FLOAT, mspace, espace );
            
        }
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


int get_pulse_attr( H5std_string filename, double * TauA, double * TauB, double * TauC, int order, int type )
{
    char probTableName[100];
    
    if (type==PEAK_PPU_TENSOR) sprintf(probTableName, "Pr%d", order);
    else if (type==TAIL_PPU_TENSOR) sprintf(probTableName, "TailPr%d", order);
    else if (type==MIX_PPU_TENSOR) sprintf(probTableName, "MixPr%d", order);
    else {
        return -1;
    }
    
    try
	{
		Exception::dontPrint();
        
        H5File file( filename, H5F_ACC_RDONLY );
        DataSet dataset = file.openDataSet( probTableName );
        DataSpace dataspace = dataset.getSpace();
    

        Attribute attA = dataset.openAttribute((const char*)"PeakEnd");
        Attribute attB = dataset.openAttribute((const char*)"TailStart");
        Attribute attC = dataset.openAttribute((const char*)"TailEnd");
        
        double ta,tb,tc;
        
        attA.read(PredType::NATIVE_DOUBLE, &ta);
        attB.read(PredType::NATIVE_DOUBLE, &tb);
        attC.read(PredType::NATIVE_DOUBLE, &tc);
        
        
        *TauA = ta;
        *TauB = tb - ta;
        *TauC = tc - tb;
        
        file.close();
        
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



