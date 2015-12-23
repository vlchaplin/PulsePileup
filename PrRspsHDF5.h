/*
 *  PrRspsHDF5.h
 *  C_ppu
 *
 *  Created by vchaplin on 3/29/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef PRRSP_HDF5_H
#define PRRSP_HDF5_H


#include "H5Cpp.h"
#include "PulseShape.h"
#include <iostream>
#include <iomanip>

using namespace std;

//#ifndef H5_NO_NAMESPACE
using namespace H5;
//#endif

#define PEAK_PPU_TENSOR 1
#define TAIL_PPU_TENSOR 2
#define MIX_PPU_TENSOR 3


/****************************
 
 Create2D - Allocate a 2D array with dimensions N1xN2 in contigious heap memory
 
 In the client program, the array is allocated by calling
 
 float *** array = Create2D(N1, N2);
 
 and then array can be accessed with familiar three-bracket syntax:  array[i][j]  = element (i,j)
 
 
 
 Since the array is in contiguous space, it can also be used as a 'flat' array:
 
 float * flatArr = array[0];  //return pointer to start of data-block
 
 flatArr[i*N2 + j]  = element (i,j)
 
 *****************************/

//0-based indexing!!
template <class T> T **Create2D(int n, int m)
{
    T **array = new T *[n + 1];
    array[0] = new T[n * m + 1]; 
    for(int i = 1; i < n; i++) array[i] = array[i - 1] + m;
    
    return array;
};
template <class T> void Delete2D(T **array) {delete[] array[0]; delete[] array;};




/****************************

Create3D - Allocate a 3D array with dimensions N1xN2xN3 in contigious heap memory

In the client program, the array is allocated by calling
 
 float *** array = Create3D(N1, N2, N3);
 
and then array can be accessed with familiar three-bracket syntax:  array[i][j][k] is equivalent to element (i,j,k)
 

 
Since the array is in contiguous space, it can also be used as a 'flat' array:
 
 float * flatArr = array[0][0];  //return pointer to start of data-block
 
 flatArr[i*N3*N2 + j*N3 + k]  = element (i,j,k)

 
 
 //  Note,  Create3D<T>(N1,N2,N3)  returns a triple pointer to an array of type T with size N1 x N2 x N3.
 //  It also creates two light-wieght arrays which provide triplet access operators ([i][j][k] sytax). 
 //  These are allocated dynamically on the heap.
 //
 //  Delete3D must be used on the returned T*** object to fully free all of this memory.
 //
 //  The returned pointer is structured so that it can be accessed using convenient 3D subscripts, but also so that 
 //  the allocated memory is contiguous (C++ doesn't do this automatically).  
 //  e.g.,
 //      T *** a3D = Create3D<T>(N1,N2,N3);
 //      a3D[0][0][0] = ...; a3D[0][0][1] = ...; a3D[0][1][0] = ...;  etc. for all a3D[i][j][k], i = 0 .. N1-1, etc.
 //  
 //      Delete3D(a3D);  !!! Must be freed with this routine.  
 //
 //  An equivalent 'flat' indexed array can be obtained from the returned object:
 //
 //      T *** a3D = Create3D<T>(N1,N2,N3);
 //      T * aFlat = a3D[0][0];
 //  which can be accessed like:
 //      aFlat[i*N2*N3 + j*N3 + k] <---  equivalent to a3D[i][j][k]
 //
 
*****************************/


//0-based indexing!!
template <class T> T ***Create3D(int N1, int N2, int N3)
{
    T *** array = new T ** [N1];
    
    array[0] = new T * [N1*N2];
    
    array[0][0] = new T [N1*N2*N3];
    
    int i,j;
    int k;
    
    for( i = 0; i < N1; i++) {
        
        if (i < N1 -1 ) {
            
            array[0][(i+1)*N2] = &(array[0][0][(i+1)*N3*N2]);
            
            array[i+1] = &(array[0][(i+1)*N2]);
            
        }
        
        for( j = 0; j < N2; j++) {
            
            if (j > 0) array[i][j] = array[i][j-1] + N3;
            
            for( k = 0; k < N3; k++) {
                //*(array[0][0] + i*N3*N2 + j*N3 + k) = i*N3*N2 + j*N3 + k;
                array[i][j][k] = 0;
            }
            
        }
        
    }

    return array;
};
template <class T> void Delete3D(T ***array) {
    delete[] array[0][0]; 
    delete[] array[0];
    delete[] array;
};




int new_prsp_hdf5( H5std_string filename, hsize_t n1, hsize_t n2, hsize_t n3, float * data, float * edge1, float *edge2, int type=PEAK_PPU_TENSOR );
int add_prsp_hdf5( H5std_string filename, hsize_t n1, hsize_t n2, hsize_t n3, float * data, int order, int type=PEAK_PPU_TENSOR );
int attach_pulse_attr( DataSet& dset );

int get_prsp_dims_hdf5( H5std_string filename, int * ndims, hsize_t *& dims, int order, int type );
int read_prsp_hdf5( H5std_string filename, hsize_t * dims, float * data, int order, int type=PEAK_PPU_TENSOR, float * eIN=NULL, float * eOUT=NULL );

int read_edges(H5File& fileObj, float * eIN = NULL, float * eOUT = NULL, bool alloc=false);

int get_pulse_attr( H5std_string filename, double * TauA, double * TauB, double * TauC, int order=1, int type=PEAK_PPU_TENSOR );




#endif