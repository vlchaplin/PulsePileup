//
//  AnalyticalPPU.h
//  C_ppu
//
//  Created by vchaplin on 3/2/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef C_ppu_AnalyticalPPU_h
#define C_ppu_AnalyticalPPU_h

#include <iostream>
#include <cstdio>
#include "SpectrumHDF5.h"
#include "PrRspsHDF5.h"

#include "ADC.h"

#include "ChannelEnergy.h"
#include "EdgeSet.hh"

//#include "SpectralShapes.h"

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_interp.h>


using namespace std;


//way more than will actually be used... little memory overhead (just extra pointers)
//Actual number is equal to the max order calculated by PPUProb.cpp and stored in the Hdf file, (around 5-7).
//Just a convenience, so the handles don't have to be dynamically allocated 
#define PR_SET_MAX 20

//#define PRINT_STATES_COMP

class PPUModelData {
    
protected:
    //These represent 3 sets of 'Pileup likelihood functions', which are the
    //  Pr(E | E0, E1) objects derived in the paper. Since we're in the discrete
    //  case they are 3D arrays, calculated by PPUProb.cpp, and stored in the HDF5 file 
    float *** Pr_AA_set[PR_SET_MAX]; //peak
    float *** Pr_AC_set[PR_SET_MAX]; //tail
    float *** Pr_BC_set[PR_SET_MAX]; //'mix'
    
    int readLevel;
    int orderLoaded;
    
public:
    
    int maxOrdInFile;
    int maxOrdToCalc;
    
    string h5file;
    
    EdgeSet<float> * ebounds;
    float * edges;
    float * a_spectrum;
    
    hsize_t dims[3];
    
    double TA;
    double TB;
    double TC;
    
public:
    PPUModelData(){
        
        for (int i=0; i < PR_SET_MAX; i++)
        {
            Pr_AA_set[i]=NULL;
            Pr_AC_set[i]=NULL;
            Pr_BC_set[i]=NULL;
        }
        ebounds=NULL;
        edges=NULL;
        a_spectrum=NULL;
        readLevel=0;
    };
    
    PPUModelData(string file){
        for (int i=0; i < PR_SET_MAX; i++)
        {
            Pr_AA_set[i]=NULL;
            Pr_AC_set[i]=NULL;
            Pr_BC_set[i]=NULL;
        }
        ebounds=NULL;
        edges=NULL;
        a_spectrum=NULL;
        readLevel=0;
        
        h5file = file;
    };
    ~PPUModelData() {
        delete_prsps();
        
    }
    
    inline double state_prob(int& a, int& b, int& c, double& cpus) {
        return pow(cpus, a+b+c)*pow(TA, a)*pow(TB, b)*pow(TC, c)*exp( -cpus*(TA+TB+TC) ) / ( gsl_sf_fact(a)*gsl_sf_fact(b)*gsl_sf_fact(c) );
    };
    inline double overlap_prob(int c, double& cpus) {
        if ( c <=0 ) return 0;
        else
            return pow(cpus, c+c-1)*pow(TC, c)*pow(TA, c-1)*exp( -cpus*(TA+TC) ) / ( gsl_sf_fact(c)*gsl_sf_fact(c-1) );
    };
    
    //MaxOrdInFile is the maximum order of likelihood functions-- Pr(E| E0,E1)-- that were
    //calculated and stored by the program PPUProb
    //
    //...MaxOrdToCalc is the approximation order that will be used in calculating the pileup model. For an input rate
    //    of 1 MHz, the average pileup order is (rate x pulse-width) = (1 MHz x 4.5us) = 4.5, 
    //    and is poisson distributed about that average. Since sources at / above 1MHz probably cause additional instrumental
    //    effects and are not well-modelled, a Max. order of 9-10 should be sufficient for most cases (~2.1 sigma for 1Mhz).
    //
    //  If MaxOrdToCalc is greater than MaxOrdInFile, the highest order likelihood is used at each order above MaxOrdInFile.
    void SetMaxOrdInFile(int m) {
        maxOrdInFile = m;
    };
    void SetMaxOrdToCalc(int n) {
        maxOrdToCalc = n;
    };
    void set_filename(string f) {
        this->h5file = f;
    }
    string get_filename() {
        return h5file;
    }
    
    hsize_t numChannels() {
        return dims[2];
    }
    
    int read_and_set_maxorder( );
    int read_file_dims( int order, int type );
    int read_file_ebounds( );
    int read_file_pulse_width( int order=1, int type=PEAK_PPU_TENSOR );
    
    float *& edges_ptr() { return this->edges; };
    
    //Interface to the code which allocates and structures the 3D arrays and reads them from the file
    //! Loads the contents of the "pileup response" file, which is a specially formatted HDF5 file with numerous data extensions and keywords.
    /*! 
    The file contains probability functions for peak-peak, peak-tail, and dead-tail pulse-pileup interactions, describing how pulse-heights are restributed.
    For each type of interaction there is an extension corresponding to the pileup order assumed.  Typically the max order stored in the file, so in this case there are 3 x 7 extensions.
     
    They are stored in the Pr_AA_set, Pr_AC_set, and Pr_BC_set, respectively.  (For example, Pr_AA_set[1] is the address of the 3-D array representing the 1st-order peak-peak pileup restribution function). Thus, <tt> Pr_AA_set[1][i][j][k] </tt>
    is the probability that, given two counts in region A, the first having pulse height in channel i, the next in channel j, a single count is measured in channel k.
     
     The functions CalcSpectrumExpansion() handles such computations once the input spectrum has been defined.
     
     */
    int load_prsps();
    void delete_prsps();
    
    //! Calculate the pulse-state spectral components referenced by P_, Q_ components. Requires that the 0th order (aka, "input") spectrum, whose address is P_components[0], be initialized.
    /*!
     The two pointer arguments are multi-dimensional arrays returned by allocate_sized_components(). Prior to calling this function the 0th order input spectrum must be normalized to 1.0 and written into P_components[0], where P_components[0][i] is the ith channel. The spectrum channel edges must be the same as those defined in the instance of PPUModelData, which can be obtained with PPUModelData::edges_ptr() and PPUModelData::numChannels().
     On return from this function, P_components[k] is the address for the kth order spectral state (a probability function), <tt> and P_components[k][i] </tt> is the ith channel of that array.
     Likewise, <tt> Q_components[a][b][c] </tt> is the address of the tail spectrum for state <abc>, and <tt> Q_components[a][b][c][i] </tt> is the ith channel of the <abc> spectrum.  
     Since there are no tail measurements for c=0, all Q_components[a][b][0] are set to NULL.
     \param P_components Allocated array to hold peak spectral components. Use allocate_sized_components() to allocate the P and Q arrays 
     \param Q_components Allocated array to hold tail spectral components
     \param maxOrder Order to calculate (if absent, <tt> correspondingData.maxOrdToCalc </tt> is used)
     \sa allocate_sized_components
     */
    int CalcOrderComponents( float **& P_components, float ****& Q_components, int maxOrder=-1 );
    
    int CalcSpectrumExpansion( double cpus, float ** P_components, float **** Q_components, double *& ans );
    
    int interp_spec2modelspace( float * inSpectrum, float * inEdges, int numIn, float * outSpectrum );
    int interp_spec2modelspace_integ( double * inSpectrum, double * inEdges, int numIn, double * outSpectrum );
    
    
    
    
    //! Allocates memory to hold two sets of spectra (in the form of multi-dimensional arrays) for the model computation.
    /*!
     The number of channels in each spectrum is obtained from correspondingData.numChannels().  The arrays are specially formatted and 
     will be defined via CalcOrderComponents().
     
     They can be set zero using zero_sized_components(), and freed using free_sized_components().
     
     \param correspondingData An object of type PPUModelData with which the P, Q components are associated (they will be handled according the number of channels in the object and its maxOrdToCalc if the order argument isn't used).
     \param P_components Allocated array to hold peak spectral components. 
     \param Q_components Allocated array to hold tail spectral components
     \param ordToAlloc Order to allocate (if absent, <tt> correspondingData.maxOrdToCalc </tt> is used)
     \sa zero_sized_components, free_sized_components
     */
    template<typename T>
    friend int allocate_sized_components( PPUModelData& correspondingData,  T **& P_components, T ****& Q_components, int ordToAlloc=-1 ) {
        
        int ord,a,b,c,maxOrder;
        
        if ( correspondingData.readLevel < 2 ) return -1;
        
        size_t N3 = correspondingData.numChannels();
        
        if ( ordToAlloc < 0 ) {
            maxOrder = correspondingData.maxOrdToCalc;
        } else {
            maxOrder = ordToAlloc;
        }
        
        if ( maxOrder <= 0 ) return -1;
        
        Q_components = Create3D<float *>(maxOrder+1, maxOrder+1, maxOrder+1);
        P_components = Create2D<float>(maxOrder+1, N3);

        
        for (ord = 1; ord <= maxOrder; ord++) {
            for (a=0; a<=ord; a++) {
                
                b = ord - a;
                Q_components[a][b][0] = NULL;
                
                for (c=1; c<=ord-a; c++) {
                    b = ord - (c+a);
                    Q_components[a][b][c] = new float[N3];
                }
            }  
        }
     
        return 0;
    };
    
    //! Frees the arrays allocated by allocate_sized_components()
    /*!
     \param correspondingData An object of type PPUModelData with which the P, Q components are associated (they will be handled according the number of channels in the object and its maxOrdToCalc if the order argument isn't used).
     \param P_components Peak components to be freed. 
     \param Q_components Tail components to be freed.
     \param ordAllocated Order to allocated in the arrays when allocate_sized_components was called  (if absent, <tt> correspondingData.maxOrdToCalc </tt> is used)
     \sa zero_sized_components, free_sized_components
     */
    template<typename T>
    friend int zero_sized_components( PPUModelData& correspondingData,  T **& P_components, T ****& Q_components, int ordAllocated=-1 ) {
        
        size_t ord,maxOrder;
        
        if ( correspondingData.readLevel < 2 ) return -1;
        
        size_t a,b,c,k;
        size_t N3 = correspondingData.numChannels();
        
        maxOrder = correspondingData.maxOrdToCalc;
        
        if ( ordAllocated < 0 ) {
            maxOrder = correspondingData.maxOrdToCalc;
        } else {
            maxOrder = ordAllocated;
        }
        
        for (k=0;k<N3;k++)
            P_components[0][k] = 0;
        
        for (ord = 1; ord <= maxOrder; ord++) {
            for (a=0; a<=ord; a++) {
                
                b = ord - a;
                Q_components[a][b][0] = NULL;
                
                for (c=1; c<=ord-a; c++) {
                    b = ord - (c+a);
                    
                    for (k=0;k<N3;k++)
                        Q_components[a][b][c][k] = 0;
                }
            }
            
            for (k=0;k<N3;k++)
                P_components[ord][k] = 0;
        }
        
        return 0;
    };

    //! Writes zeros into each channel of the spectral components 
    /*!
     \param correspondingData An object of type PPUModelData with which the P, Q components are associated (they will be handled according the number of channels in the object and its maxOrdToCalc if the order argument isn't used).  
     \param P_components Peak components to be zeroed. 
     \param Q_components Tail components to be zeroed.
     \param ordAllocated Order to allocated in the arrays when allocate_sized_components was called  (if absent, <tt> correspondingData.maxOrdToCalc </tt> is used)
     \sa zero_sized_components, free_sized_components
     */
    template<typename T>
    friend int free_sized_components( PPUModelData& correspondingData,  T **& P_components, T ****& Q_components, int ordAllocated=-1 ) {
        
        int ord,a,b,c,maxOrder;
        
        if ( correspondingData.readLevel < 2 ) return -1;
        
        maxOrder = correspondingData.maxOrdToCalc;
        
        if ( ordAllocated < 0 ) {
            maxOrder = correspondingData.maxOrdToCalc;
        } else {
            maxOrder = ordAllocated;
        }
        
        if ( P_components != NULL) Delete2D(P_components);
        
        if ( Q_components == NULL ) return 0;
        
        for (ord = 1; ord <= maxOrder; ord++) {
            for (a=0; a<=ord; a++) {
                
                for (c=1; c<=ord-a; c++) {
                    b = ord - (c+a);
                    
                    if ( Q_components[a][b][c] != NULL )
                        delete [] Q_components[a][b][c];
                }
            }  
        }
        
        Delete3D(Q_components);
        
        return 0;
    };
    
};

template<typename T>
inline int rebin_pdf_integ( T * inEdge, T * inAvgPDF, size_t numInChans, T * outEdge, T *& outHist, size_t numOutChans )
{
    T iXa, iXb, oXa, oXb, dX;
    
    size_t i, o;
    i=0;
    o=0;
    iXa = inEdge[i]; iXb = inEdge[i+1];
    oXa = outEdge[o]; oXb = outEdge[o+1];
    
    while (i < numInChans && o < numOutChans) {
        
        iXa = inEdge[i]; iXb = inEdge[i+1];
        oXa = outEdge[o]; oXb = outEdge[o+1];
        
        if ( oXa > iXb ) {
            i++;
            continue;
        }
        if ( iXa > oXb ) {
            o++;
            continue;
        }
        
        if ( oXa < iXa && oXb > iXa ) {
            
            if (i==1) {
                dX = iXa - oXa;
                outHist[o] += inAvgPDF[i-1]*dX;
                //cout << "1:  " << o << ", "<< i-1 <<  "  " << dX << endl;
            }
            
            while ( oXb >= iXb && i < numInChans ) {
                
                dX = iXb - iXa;
                outHist[o] += inAvgPDF[i]*dX;
                //cout << "2a:  " << o << ", "<< i << "  " << dX << endl;
                i++;
                iXa = inEdge[i]; iXb = inEdge[i+1];
            }
            
            if ( i > 1 && iXa < oXb ) {
                dX = oXb - iXa;
                outHist[o] += inAvgPDF[i]*dX;
                //cout << "2b:  " << o << ", "<< i << "  " << dX << endl;
            }
            
            o++;
        }
        else if ( oXa >= iXa && oXb <= iXb ) {
            outHist[o] += inAvgPDF[i]*(oXb - oXa);
            //cout << "3:  " << o << ", " << i << "  " << (oXb - oXa) << endl;
            o++;
        }
        else if ( oXa >= iXa && oXb >= iXb ) {
            outHist[o] += inAvgPDF[i]*(iXb - oXa);
            //cout << "4:  " << o << ", " << i << "  " << (iXb - oXa) << endl;
            i++;
        }
        
    }
    
    return 0;
}

template<typename T>
inline void interpolate1(double * oX, double * oY, size_t oN, T * nX, T * nY, size_t nN)
{
    size_t k;
    gsl_interp * gsl_interp_obj = gsl_interp_alloc(gsl_interp_linear, oN );
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp_init(gsl_interp_obj, oX, oY, oN);
    
    for (k=0; k<nN; k++)
    {
        nY[k] = gsl_interp_eval(gsl_interp_obj, oX, oY, nX[k], acc);
    }
    
    gsl_interp_free (gsl_interp_obj);
    gsl_interp_accel_free (acc);
    return;
}


template<typename T>
inline void interpolate0_integ(double * X_edges, double * Y_at_edges, size_t numEdges, T * nX_edges, T * nY, size_t numOutChans)
{
    size_t k;
    gsl_interp * gsl_interp_obj = gsl_interp_alloc(gsl_interp_linear, numEdges );
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp_init(gsl_interp_obj, X_edges, Y_at_edges, numEdges);
    
    double nXa, nXb;
    
    for (k=0; k<numOutChans; k++)
    {
        nXa = nX_edges[k];
        nXb = nX_edges[k+1];
        
        nY[k] = (T)gsl_interp_eval_integ(gsl_interp_obj, X_edges, Y_at_edges, nXa, nXb, acc);
        
    }
    
    gsl_interp_free (gsl_interp_obj);
    gsl_interp_accel_free (acc);
    return;
}



//These are also implemented in the class so the values of the TA,TB,TC do not have to handled by clients
inline double state_prob(int& a, int& b, int& c, double& cpus, double& TA, double& TB, double& TC) {
    return pow(cpus, a+b+c)*pow(TA, a)*pow(TB, b)*pow(TC, c)*exp( -cpus*(TA+TB+TC) ) / ( gsl_sf_fact(a)*gsl_sf_fact(b)*gsl_sf_fact(c) );
}
inline double overlap_prob(int c, double& cpus, double& TA, double& TC) {
    if ( c <=0 ) return 0;
    else
        return pow(cpus, c+c-1)*pow(TC, c)*pow(TA, c-1)*exp( -cpus*(TA+TC) ) / ( gsl_sf_fact(c)*gsl_sf_fact(c-1) );
}







#endif
