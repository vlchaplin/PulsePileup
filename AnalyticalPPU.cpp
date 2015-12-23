//
//  AnalyticalPPU.cpp
//  C_ppu
//
//  Created by vchaplin on 3/2/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "AnalyticalPPU.h"


/***
Requires that the zeroth order spectrum, P_components[0], be already defined
****/
int PPUModelData::CalcOrderComponents( float **& P_components, float ****& Q_components, int maxOrder )
{
    if (maxOrder < 0) maxOrder = maxOrdToCalc;
    
    if (readLevel < 3) return -1;
    
    int ord,a,b,c,i,j,k;
    int N1 = dims[0];
    int N2 = dims[1];
    int N3 = dims[2];
    
    float *** Pr_AA;
    float *** Pr_AC;
    float *** Pr_BC;
    
    float P_norm, Q_001_norm;
    
#ifdef PRINT_STATES_COMP
    static char msg[100];
#endif
    
    for (ord = 1; ord <= maxOrder; ord++) {
        
        //cout << "read status " << read_prsp_hdf5( filename, dims, Pr_AA[0][0], ord, PEAK_PPU_TENSOR, NULL, NULL ) << endl;
        
#ifdef PRINT_STATES_COMP
        sprintf(msg, "Doing P[%d], Q[0,0,%d]", ord, ord);
        cout << msg << endl;
#endif
    
        if ( ord> orderLoaded ) {
            Pr_AA = Pr_AA_set[orderLoaded];
        }
        else {
            Pr_AA = Pr_AA_set[ord];
        }
        
        
        
        // P(a>0)
        for (i=0; i<N1; i++) {
            for (j=0; j<N2; j++) {
                for (k=0; k<N3; k++) {
                    //Eqn (36)
                    P_components[ord][k] += (Pr_AA[i][j][k]*P_components[ord-1][i]*P_components[0][j]);
                }
            }
        }
        
        //The peak components must all normalize to 1.0. This accounts for small loss of probability due to 
        //the discrete integration used to calculate likelihood functions (say ~0.991 instead of 1.0).  
        //Such small discrepancies exponentiate since this is an iterative calculation
        P_norm=0.0;
        for (k=0; k<N3; k++)
            P_norm += P_components[ord][k];
        for (k=0; k<N3; k++)
            P_components[ord][k] /= P_norm;
        
        // Q(0,0,1)
        if ( ord == 1 ) {
            
            Q_001_norm=0;
            //read_prsp_hdf5( filename, dims, Pr_AC[0][0], ord, TAIL_PPU_TENSOR, NULL, NULL );
            Pr_AC = Pr_AC_set[ord];
            
            for (i=0; i<N1; i++) {
                for (j=0; j<N2; j++) {
                    for (k=0; k<N3; k++) {
                        //Eqn (41)
                        Q_components[0][0][ord][k] += (Pr_AC[i][j][k]*P_components[0][i]*P_components[0][j]);
                        Q_001_norm += (Pr_AC[i][j][k]*P_components[0][i]*P_components[0][j]);
                    }
                }
            }
            
            //Q<001> must be re-normalized so it's 1.0
            for (k=0; k<N3; k++)
                Q_components[0][0][ord][k] /= Q_001_norm;
            
        }
        
        if ( ord < maxOrder ) {                
            // Q(0,0,c>1)
            for (i=0; i<N1; i++) {
                for (j=0; j<N2; j++) {
                    for (k=0; k<N3; k++) {
                        //Eqn (43)
                        Q_components[0][0][ord+1][k] += (Pr_AA[i][j][k]*Q_components[0][0][ord][i]*Q_components[0][0][1][j]);
                    }
                }
            }
        }
    }
    
    float * secondary_spec;
    
    
    for (c=1; c <= maxOrder; c++) {
        //read_prsp_hdf5( filename, dims, Pr_AC[0][0], c, TAIL_PPU_TENSOR );
        //read_prsp_hdf5( filename, dims, Pr_BC[0][0], c, MIX_PPU_TENSOR );
        
        if ( c > orderLoaded )
        {
            Pr_AC = Pr_AC_set[orderLoaded];
            Pr_BC = Pr_BC_set[orderLoaded];
        } else {
            Pr_AC = Pr_AC_set[c];
            Pr_BC = Pr_BC_set[c];
        }
        
        if (c==1) 
            secondary_spec = P_components[0]; //for Eqn (42)
        else
            secondary_spec = Q_components[0][0][c-1]; //for Eqn (44)
        
        
        // Q(a>0,0,c>=1)
        b=0;
        for (a=1; a <= maxOrder-c ; a++) {
#ifdef PRINT_STATES_COMP
            sprintf(msg, "Doing Q[%d,%d,%d]", a, b, c);
            cout << msg << endl;
#endif
            for (i=0; i<N1; i++) {
                for (j=0; j<N2; j++) {
                    for (k=0; k<N3; k++) {
                        
                        //Eqn (42 or 44)
                        Q_components[a][0][c][k] += (Pr_AC[i][j][k] * P_components[a][i] * secondary_spec[j]);
                        
                    }
                }
            }
        }
        
        // Q(a>0,b>0,c>=1)
        for (ord=2; ord <= maxOrder; ord++) {
            
            for (b=1; b <= ord-c ; b++) {
                
                a = ord - (b+c);
#ifdef PRINT_STATES_COMP               
                sprintf(msg, "Doing Q[%d,%d,%d]", a, b, c);
                cout << msg << endl;
#endif                
                for (i=0; i<N1; i++) {
                    for (j=0; j<N2; j++) {
                        for (k=0; k<N3; k++) {
                            //Eqn (47)
                            Q_components[a][b][c][k] += ( Pr_BC[i][j][k] * P_components[b-1][i] * Q_components[a][0][c][j] );
                        }
                    }
                }  
            }
        }
        
    }
    
    return 0;
};

int PPUModelData::CalcSpectrumExpansion( double cpus, float ** P_components, float **** Q_components, double *& ans )
{
        
    if (readLevel < 3) return -1;
    
    int o=0,m,a,b,c;
    double DT = TA+TB+TC;
    double loopProb;
    
    int maxOrd = this->maxOrdToCalc;
    int nel = (int)this->numChannels();
    
    
    //Equation (61) & (63) from the PPU paper are combined below
    
    double TailOverProb=2.0;
    for (o=1; o <= 20; o++)
        TailOverProb -= this->overlap_prob(o, cpus);
    
    //cout << endl << "here 1" << endl;
    for (o=0; o <= maxOrd; o++) {
        
        //cout << endl << "order " << o <<  endl;
        if (o==0) {
            for (m=0; m<nel; m++) {
                ans[m] = exp(-cpus*DT)*(double)P_components[o][m];
                
                //cout << ans[m] << " ";
                
            }
        } else {
            
            for (a=0; a<=o; a++) {
                
                //Peak part
                for (b=0; b<=o-a; b++) {
                    c = o - (b+a);
                    loopProb = (1.0/(o+1.0))*this->state_prob(a,b,c,cpus);
                    for (m=0; m<nel; m++) {
                        ans[m] += loopProb*(double)P_components[a][m];
                    }
                }
                
                //Tail part
                for (c=1; c<=o-a; c++) {
                    b = o - (c+a);
                    
                    if ( Q_components[a][b][c] == NULL ) continue;
                    
                    loopProb = TailOverProb * (1.0/(o+1.0))*this->state_prob(a,b,c,cpus);
                    
                    for (m=0; m<nel; m++) {
                        ans[m] += loopProb*(double)Q_components[a][b][c][m];
                    }
                }
            }   
        }
        
        for (m=0; m<nel; m++) {
            ans[m] -= (1.0/(o+1.0))*this->overlap_prob(o+1, cpus)*(double)P_components[o][m];
        }
        
    }
 
    
    
    
    
    
    return 0;
}



//------------------------------------
//
// Member function:   PPUModelData::load_prsps()
//
//
//Read the data sets.  the function read_prsp_hdf5() is defined in PrRspsHDF5.cpp.
// Create3D is defined in PrRspsHDF5.h
//
//  Create3D<T>(N1,N2,N3)  returns a triple pointer to an array of type T with size N1 x N2 x N3.
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
//  However, `delete [] aFlat` is insufficient to free all the acquired memory. Delete3D(a3D) must be used instead.
//  
//  The PPUModelData class encapsulates all of these operations.  load_prsps() allocates and read.  delete_prsps() frees.
//
//
//  Because the HDf5 library expects a flat array when it reads a file, the flat address is passed to read_prsp_hdf5().
//------------------------------------

int PPUModelData::load_prsps()
{
    if (readLevel < 1) return -1;
    
    int N1 = dims[0];
    int N2 = dims[1];
    int N3 = dims[2];
    
    int rc =0;
    int maxOrdToRead = maxOrdInFile;
    
    if (maxOrdToCalc < maxOrdInFile) maxOrdToRead = maxOrdToCalc; 
    
    
    orderLoaded=0;
    // no 'zeroth order' Pr functions (so i=1 ...)
    for (int i=1; i <= maxOrdToRead; i++)
    {
        Pr_AA_set[i]=Create3D<float>(N1, N2, N3);
        Pr_AC_set[i]=Create3D<float>(N1, N2, N3);
        Pr_BC_set[i]=Create3D<float>(N1, N2, N3);
        
        rc |= read_prsp_hdf5(this->h5file, this->dims, Pr_AA_set[i][0][0], i, PEAK_PPU_TENSOR);
        rc |= read_prsp_hdf5(this->h5file, this->dims, Pr_AC_set[i][0][0], i, TAIL_PPU_TENSOR);
        rc |= read_prsp_hdf5(this->h5file, this->dims, Pr_BC_set[i][0][0], i, MIX_PPU_TENSOR);
        
        orderLoaded++;
    }
    
    if (rc == 0 && readLevel<3) {
        a_spectrum = new float[N3];
        readLevel=3;
    }
    
    return rc;
}

void PPUModelData::delete_prsps()
{
    // no 'zeroth order' Pr functions (so i=1 ...)
    for (int i=1; i <= orderLoaded; i++)
    {
        
        if ( Pr_AA_set[i] != NULL ) Delete3D(Pr_AA_set[i]);
        if ( Pr_AC_set[i] != NULL ) Delete3D(Pr_AC_set[i]);
        if ( Pr_BC_set[i] != NULL ) Delete3D(Pr_BC_set[i]);

    }
    
    if ( edges != NULL ) delete [] edges;
    if ( ebounds != NULL ) delete ebounds;
    
    if ( a_spectrum != NULL ) delete [] a_spectrum;
    
    orderLoaded=0;
}

int PPUModelData::read_and_set_maxorder( )
{
    H5std_string filename(this->get_filename());
    
    int maxOrd=-1;
    
    
    try
	{
		Exception::dontPrint();
        
        H5File h5FileObj( filename, H5F_ACC_RDONLY );
        
        DataSet miscDset = h5FileObj.openDataSet("misc");
        
        Attribute atr = miscDset.openAttribute("MaxOrder");
        
        atr.read(PredType::NATIVE_INT, &maxOrd);
        
        h5FileObj.close();
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
    
    
    this->SetMaxOrdInFile(maxOrd);
    return maxOrd;
}

int PPUModelData::read_file_dims( int order, int type )
{
    
    int ndims;
    
    hsize_t* dptr;
    
    H5std_string filename(this->get_filename());
    
    int status = get_prsp_dims_hdf5( filename, &ndims, dptr, order, type);
    
    if (status != 0) return status;
    
    readLevel=1;
    
    this->dims[0] = dptr[0];
    this->dims[1] = dptr[1];
    this->dims[2] = dptr[2];
    
    
    delete [] dptr;
    
    return status;
    
}

int PPUModelData::read_file_ebounds( )
{
    H5std_string filename(this->get_filename());
    int status;
    if ( readLevel < 1 ) return -1;
    
    try
	{
		Exception::dontPrint();
        
        H5File file( filename, H5F_ACC_RDONLY );
        
        edges = new float[this->dims[2]+1];
        
        status = read_edges(file, edges, NULL, 0 );
        
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
    
    if (status == 0 && readLevel<2) readLevel=2;
    return status;
    
}
int PPUModelData::read_file_pulse_width( int order, int type)
{
 
    H5std_string filename(this->get_filename());
    
    return get_pulse_attr( filename, &TA, &TB, &TC, order, type );
    
}





int PPUModelData::interp_spec2modelspace( float * inSpectrum, float * inEdges, int numIn, float * outSpectrum )
{
    
    if (readLevel < 2) return -1;
    
    size_t N3 = (size_t)this->dims[2];
    float minX = inEdges[0];
    float maxX = inEdges[numIn];
    float mX;
    if (minX > this->edges[N3] ) return -1;
    if (maxX < this->edges[0] ) return -1;
    
    int i;
    size_t k;
    double * X = new double[numIn];
    double * Y = new double[numIn];
    
    for (i=0; i< numIn; i++) {
        X[i] = (inEdges[i] + inEdges[i+1])/2.0;
        Y[i] = inSpectrum[i];
    }
    
    gsl_interp * gsl_interp_obj = gsl_interp_alloc(gsl_interp_linear, (size_t)numIn );
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp_init(gsl_interp_obj, X, Y, (size_t)numIn);
    
    for (k=0; k<N3; k++)
    {
        mX = this->edges[k];
        if ( mX >= minX && mX <= maxX )
            outSpectrum[k] = gsl_interp_eval(gsl_interp_obj, X, Y, mX, acc);
        else 
            outSpectrum[k] = 0;
    }
        
    delete [] X;
    delete [] Y;
    gsl_interp_free (gsl_interp_obj);
    gsl_interp_accel_free (acc);
    return 0;
}



// Convert an input differential spectrum (represented in discrete form) to an array of channel values.
// If the input spectrum is a PDF (normalized to 1.0), the output is equivalent to a discrete probability distribution.
//
// The input is assumed to be the average value of the diff. spectrum or PDF in each channel
//  
int PPUModelData::interp_spec2modelspace_integ( double * inSpectrum, double * inEdges, int numIn, double * outSpectrum )
{
    if (readLevel < 2) return -1;
    
    size_t N3 = (size_t)this->dims[2];
    double minX = inEdges[0];
    double maxX = inEdges[numIn];
    double mXa,mXb;
    if (minX > this->edges[N3] ) return -1;
    if (maxX < this->edges[0] ) return -1;
    

    size_t k;
    double * X = new double[numIn+1];
    double * Y = new double[numIn+1];
    
    for (k=0; k< (size_t)numIn; k++) {
        
        //Y[i] is the average value of the PDF in the interval X[i] to X[i+1]
        
        X[k] = inEdges[k];
        Y[k] = inSpectrum[k];
    }
    X[k] = inEdges[k];
    Y[k] = inSpectrum[k-1]; //repeat the last one for the upper boundary
    
    gsl_interp * gsl_interp_obj = gsl_interp_alloc(gsl_interp_linear, (size_t)numIn+1 );
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp_init(gsl_interp_obj, X, Y, (size_t)numIn+1);
    
    for (k=0; k<N3; k++)
    {
        mXa = this->edges[k];
        mXb = this->edges[k+1];
        if ( mXa >= minX && mXb <= maxX )
            outSpectrum[k] = gsl_interp_eval_integ(gsl_interp_obj, X, Y, mXa, mXb, acc);
        else 
            outSpectrum[k] = 0;
    }
    
    delete [] X;
    delete [] Y;
    gsl_interp_free (gsl_interp_obj);
    gsl_interp_accel_free (acc);
    return 0;
}


