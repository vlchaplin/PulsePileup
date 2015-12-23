//
//  ppu_static_interface.h
//  C_ppu
//
//  Created by vchaplin on 3/14/13.
//

/** \file
 * PPU C-style interface to C++ routines.
 */

#ifndef C_ppu_ppu_static_interface_h
#define C_ppu_ppu_static_interface_h

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <map>
#include <time.h>

#include "ADC.h"
#include "AnalyticalPPU.h"
#include "ChannelEnergy.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef DECORRELATE_NUM
#define DECORRELATE_NUM 3
#endif

using namespace std;



/**********************************************************
 
 PPU analytical model interface (C-style interface routines are at the bottom)
 
 **********************************************************/



extern int ApproximationOrder;

#define PPU_MC_FIXED_NUMBER 1
#define PPU_MC_FIXED_EXPOSURE 2
#define PPU_MC_FIXED_DETECTION_NUMBER 3

struct FileData {

    PPUModelData * ppuData;
    int NumRefs;
};

struct ModelData {
    string h5file;
    int maxOrderAllocated;
    int currentOrder;
    
    FileData * h5_shared_data;
    
    
    float ** P_components;
    float **** Q_components;
    double * tempX;
    double * tempY;
    double * tempDE;
    
};

// this is just a very c++'ey kind of hash. maps filenames to memory objects that will hold their data.
// Fast access, ~O(log N).
typedef map<string, FileData *> map_type;
typedef map<string, ModelData *> model_map_type;
static map_type fileInstances; 
static model_map_type modelInstances; 



FileData * new_ppufile_instance( char * filename );
int get_ppufile_node( char * file, FileData ** node);
int get_instance_node( char * keyname, ModelData ** node, bool auto_create=true);
int get_instance( char * filename, PPUModelData **);
int ppu_compute_distorted_spectrum_modelbasis( char * keyname, double cpus, double *& outCounts, double *& outEdges, long * outN);

template<typename T>
inline long channel_search( T * edge_array, T value, long numElements )
{
    
    long median;
    long min_idx, max_idx;
	double mvalue;
    
    bool found = 0;
    
    min_idx = 0;
    max_idx = numElements-1;
	
    if ( value < edge_array[min_idx] || value > edge_array[max_idx] ) return -1;
    
    while (! found) {
        
        if ( min_idx > max_idx || min_idx < 0 ) return -1;
        //if ( min_idx == max_idx ) return min_idx;
        if ( max_idx - min_idx <= 2) return min_idx;
        
        median = ( min_idx+max_idx ) / 2;
        
        mvalue = edge_array[median];
        
        if (mvalue < value ) {
            min_idx = median;
            continue;
        }
        else if (mvalue > value) {
            max_idx = median+1;
            continue;
        }
        else {
            
            
            found = 1;
            
            return min_idx;
        }
    }
    
    return -1;
};

//Only used for the MC simulation
extern ChannelEnergy * chan2en;

/**********************************************************

 ApproximationOrder is a global variable, set in the .cpp file.
 
 If not using the default value, this routine should be called before an instance
 for the file has been created.  This variable controls how much data is read from the HDF5 file.
 
 If the ApproximationOrder is less than the max order contained in the file, only 'ApproximationOrder' orders are read.
 If it's greater,all orders are read.
 
 Thus, if the source is believed to be at a low rate (high orders are uneccesary), some efficiency can be obtained
 by using a small approximation order.
 
 A guide is:
 True rate does not exceed:     good approx. order:
 < 25kHz                        Use deadtime-only correction
 < 40kHz                        1
 <100kHz                        2
 <130kHz                        3
 <200kHz                        4
 <500kHz                        6
 < 1MHz                         8
 
 These are a bit generous, but the parameter space can not be fully explored without overshooting the true rate 
 (e.g., to define a minimum), so extra orders are needed.
 
**********************************************************/

int ppu_set_approximation_order(int order);

//  BEGIN FORTRAN (C-style) INTERFACE ROUTINES 

/********************

 Anayltical PPU routines (with Doxygen markup)
 
********************/

#endif

extern "C" {

    
    
    
    //! A function to be called before ppu_load_h5_file().  Sets the maximum correction order anticipated for the lifetime of the h5 file.  
    /*!
     Smaller orders can still be used for calculations by first creating an instances associated with the file, ppu_create_instance_from_file, and then using ppu_set_approximation_order.
     Can be reset by deleting all instances bound to the file using ppu_delete_instance_data, calling this function, then the load function again.
     */
    int ppu_set_maximum_order(int order);
    
    
    
    
    //!  Load the "pileup response" file, which is a specially formatted HDF5 file with numerous keywords and extensions.
    /*!
     The pileup response file contains probability functions for peak+peak, peak+tail, and dead+tail pulse-pileup interactions, describing how pulse-heights are restributed. These are the \f$\text{Pr}(\varepsilon | E_0, E_1)\f$ objects in the paper (Eqns. 31 and 38 are for peak+peak, 40 and 45 are for peak+tail, and 46 for the dead+tail). Code used to calculate these functions is mainly in MeasuredEnergies.cpp and PPUProb.cpp (the main program), however once calculated the files are portable to any system with HDF5 installed (HDF5 is a stable library in broad circulation).  Thus, it's not necessary to re-make the files.
     
     In the file are restribution functions for numerous pileup orders (7 by default). Depending on the intensity of the GBM source, one may wish to use a higher or lower order pileup approximation.  To safely do this, the <b> maximum pileup order </b> anticipated should be set prior to calling the load routine. This is done with the single-argument function ppu_set_maximum_order(). A lower order can be used later in actual calculations, but this ensures enough memory is allocated.  A good maximum order, even for TGFs, is 7 or 8.  By default, the maximum order is 8.  If the maximum order to calculate is greater than what's in the file, an approximation will be used.
     
     This routine uses static memory to store data from the file. The loaded data (internally a C++ object encapsulating several arrays, member variables, etc.) is referenced by the filename itself in a hash that is implemented internally.  Subsequent calls of ppu_load_h5_file() with the same filename are ignored and return a zero status (non-error). A file can be re-loaded by first calling ppu_delete_instance_data() on all tags associated with this file, then calling ppu_load_h5_file() again.
     
     Data remain in memory until all references are deleted or the program stops executing.
     
     \param filename Name of the pileup response file (an HDF5 file) to load
     \return 0 if succesfull, non-zero otherwise
     */
    int ppu_load_h5_file(char * filename);
    
    
    
    
    //! Set the approximation order used in the next calculation of the model.
    /*!
     With this routine the model instance referred to by \c keyname can be set to use a certain approximation order.  This should not exceed the maximum order initially loaded (see  ppu_set_maximum_order() ). For testing lower values of the rate, small approximation orders can be used.  This greatly increases the speed of routines ppu_initialize_zeroth_order_counts() and ppu_compute_distorted_spectrum().
     \param keyname Name of the model instance whose approximation order is being set
     \param order Integer approximation order >= 1 and <= the current maximum order (8 by default).
     \sq ppu_set_maximum_order()
    */
    int ppu_set_approximation_order(const char * keyname, int order);

    
    
    //! Create an instance of the model that uses the loaded data referenced by 'filename'. 
    /*!
     This routine enables calculating the PPU spectrum for multiple detectors or for multiple spectra, wherein the pileup response data can be shared.
     Since the PPU repsonse is detector independent (only the rate and input spectrum are detector dependent), sharing it improves code performance.
     
     Each new model instance must be given a unique string tag, such as "n9" or "BGO_00", followed by the (already loaded) HDF5 file name. \c keyname can be any alpha-numeric string. \c filename must have been succesfully loaded by ppu_load_h5_file.  Internal data, such as the initialized spectral shape and approximation order, are referenced by these tags and will persist until changed or deleted.  If the tagname already exists, this routine does nothing.
     
     This tag will should be used in subsequent calls to ppu_set_approximation_order(), ppu_initialize_zeroth_order_counts(), ppu_compute_distorted_spectrum(), and ppu_delete_instance_data(). 
     
     Finally, the shared data can only be deleted once all tags referening to it are deleted.
     
     \param keyname String containing the new tag by which this model will be referenced.
     \param filename The HDF5 filename (already loaded) that this model will use in its calculation.
     \return 0 if succesful, non-zero otherwise
     */
    int ppu_create_instance_from_file( char * keyname, char * filename );

    
    
    
    //! Return the number of channels of the pileup response function
    /*!
     This is a utility routine 
     \param filename Name of the HDF5 file.  If it's not already loaded, the file will first be loaded with ppu_load_h5_file().
     \param[out] numchannels  Pointer to long type which contains the number of channels.
     \return 0 if succesful, non-zero otherwise
     */
    int ppu_get_h5_numchannels( char * filename, long * numchannels );
    
    
    
    //! Return the channel edges of the pileup response function
    /*!
     This is a utility routine.  The outpur array \c channelEdges must already be allocated to hold \c numchannels+1 elements.
     \c numchannels is returned by ppu_get_h5_numchannels().
     \param filename   Name of the HDF5 file.  If it's not already loaded, the file will first be loaded with ppu_load_h5_file().
     \param[out] channelEdges Address of double array to hold the channel edges.
     \return 0 if succesful, non-zero otherwise
     */
    int ppu_get_h5_channels( char * filename, double * channelEdges );

    
    
    int ppu_initialize_zeroth_order_pdf( char * keyname, double * inEdges, double * inSpectrum, long numInChans );

    //! Set the input spectral shape, in terms of un-normalized counts.
    /*!
     This routine initializes the model instance \c keyname. The input spectrum is interpreted as an integral spectrum in the channels defined by \c inEdges. \c numInChans is the number of \b channels, not edges.  The spectrum \b does \b not need to be normalized--this is done internally.
     
     The input spectrum is re-binned to the intrinsic binning of the pileup response (which can be determined with ppu_get_h5_channels() ).  Thus \c inEdges must be in the same units as the PPU file (e.g., keV or volts, depending on the file).
     
     \param keyname Name of the model instance being initialized
     \param inEdges Array of channel boundaries of the input spectrum
     \param inputCountSpectrum  Input spectrum (array of channel values)
     \param numInChans  Number of channels in the input spectrum
     \return 0 if succesful, non-zero otherwise
     */
    int ppu_initialize_zeroth_order_counts( char * keyname, double * inEdges, double * inputCountSpectrum, long numInChans );

    
    
    //! Evaluate the model and return the PPU-distorted spectrum given the input rate
    /*!
     Once the model has been initialized, this routine can be called to evaulate the distortion at a given rate.
     The returned spectrum is in terms of probability per channel, instead of counts.  If deadtime and/or pulse-pileup losses occur, 
     the summed spectrum will be <1.0 .  
     
     The input rate in this routine is given in \b counts \b per \b microsecond. The input spectrum should have set already by ppu_initialize_zeroth_order_counts().
     
     To compute the counts spectrum, simply multiple each channel of the output by the expected number of input counts, which is <tt> Rate * Binsize </tt>, e.g.,
     <tt> outSpec[i] *= num_input_counts </tt>.
     
     The desired binning is passed \b to this routine, which interpolates the model result.
     
     Memory for \c outSpec must be allocated before calling this routine.
     
     \param keyname Name of the model instance being calculated
     \param cpus Input rate in counts per microsecond
     \param outEdges The desired channel boundaries for the output spectrum to be computed in. 
     \param numchans The number of channels in the output spectrum
     \param[out] outSpec A pre-allocated array to hold the output channel data (probability per channel)
     \return 0 if succesful, non-zero otherwise
     */
    int ppu_compute_distorted_spectrum( char * keyname, double cpus, double * outEdges, long numchans, double * outSpec );

    
    //! Free data associated with a model instance
    /*!
     This routine will free data associated with the model instance. If multiple models share the same HDF5 file data, the fiel data persists in memory until the last model reference is deleted.
     
     \param keyname Name of the model instance being deleted
     */
    int ppu_delete_instance_data( char * keyname );
    
    
    //! Free all data loaded from HDF5 ppu files via ppu_load_h5, and free data from all model instances.
    int ppu_delete_all();
    
};

/**********************************************************
 
 Monte carlo PPU simulation interface

**********************************************************/

extern "C" {

    //! Initialize the voltage-to-energy mapping before running the monte carlo.  
    /*!
     The first argument is the energy which corresponds to the LLT (e.g., ~4.0 keV for NAIs, ~150.0 for BGO 0, ~105.0 for BGO 1).
     The second argument corresponds to 5.0 volts, the start of the overflow channel. (e.g., ~1000.0 keV for NAIs, ~45000 for BGOs).
     
     \param minEnergy Energy of the LLT
     \param maxEnergy Energy of the overflow channel
     */
    void mc_initialize_chan2en_conversion( double minEnergy, double maxEnergy );

    //! Run a monte carlo simulation with the given parameters, and return the spectrum.
    /*!
     \c throwMode determines the meaning of several arguments. Set it to one of the predifined integer symbols: \c PPU_MC_FIXED_EXPOSURE, \c PPU_MC_FIXED_NUMBER, or \c PPU_MC_FIXED_DETECTION_NUMBER (defined in ppu_static_interface.h). 
     
     For <tt> throwMode = PPU_MC_FIXED_DETECTION_NUMBER </tt>, the simulation will run until the \c mc_size. Then \c exposure_us_ptr returns the elapsed exposure time.
     For <tt> throwMode = PPU_MC_FIXED_NUMBER </tt>, the simulation will generate \c mc_size events and then process them for pulse-pileup.
     For <tt> throwMode = PPU_MC_FIXED_EXPOSURE </tt>, the simulation will generate events until a certain amount of exposure time is reached, then process these events for PPU.
     
     
     
     \param throwMode 
     \param mc_size See above
     \param exposure_us_ptr See above
     */
    long mc_simulate_ppu_spectrum( int throwMode, long * mc_size, double * exposure_us_ptr, double rate_cpus, double * inSpectrum,  double * outSpectrum, double * inEdges, long numInChans, float ** counted_data=NULL, gsl_rng ** rng = NULL );

    
};

//  END OF FORTRAN (C-style) INTERFACE ROUTINES 
    
#endif
