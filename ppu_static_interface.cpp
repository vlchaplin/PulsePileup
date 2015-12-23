//
//  ppu_static_interface.cpp
//  C_ppu
//
//  Created by vchaplin on 3/14/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "ppu_static_interface.h"

int ApproximationOrder = 8;  //order 8 is sufficient up to 1 MHz, overkill for lower rates. 

//Only used for the MC simulation
ChannelEnergy * chan2en = NULL;

FileData * new_ppufile_instance( char * filename ) {
    
    if ( filename == NULL ) return NULL;
    
    string strFile=string(filename);
    
    map_type::iterator lookupIter;
    
    lookupIter = fileInstances.find(strFile);
    
    FileData * returnVal = NULL;
    
    if ( lookupIter == fileInstances.end() ) {
        
        //Create a new interface object.  The 'PPUModelData' class encapsulates reading the HDF5 file and handling
        //it's contents in memory
        PPUModelData * newPPUprobSet = new PPUModelData;
        
        
        //just defines the filename, doesn't read it yet
        newPPUprobSet->set_filename(strFile);
        
        
        
        //Try to read the 'MaxOrder' attribute from the /misc dataset. This represents the maximum order Pr(e|E0,E1) function recorded in the file.
        //This is needed so the reader knows how many extensions to try reading.
        int maxOrd = newPPUprobSet->read_and_set_maxorder();
        
        if ( maxOrd <= 0 ) {
            //failure. there is a problem reading the file
            delete newPPUprobSet;
            return NULL;
        } else {
            //success
            
            //newPPUprobSet->SetMaxOrdInFile(maxOrd);  ** uneccesary; the read_and_set_... function already did this
            
            //set the approximation order to use instead of maxOrd.. if it's greater than maxOrd (ie, the maximum order
            //  pulse-height redistribution function stored in the file), then an approximation will be used for 
            //  orders above maxOrd.  If it's less than maxOrd, only  orders up to ApproximationOrder are read from the file
            //  or used in computing the model... this can be useful for studying low-rate events as computation requirements 
            //  exponentially increase with max order.
            newPPUprobSet->SetMaxOrdToCalc(ApproximationOrder);
        }   
        
        
        
        
        // now that we know how many extensions to read, the object will attempt to allocate and read the file data.
        // the actual number read is the smaller of 'ApproximationOrder' and 'maxOrd'.
        newPPUprobSet->read_file_dims(1, PEAK_PPU_TENSOR);
        newPPUprobSet->read_file_ebounds();
        newPPUprobSet->read_file_pulse_width();
        
        
        if ( newPPUprobSet->load_prsps() == 0 ) {
            
            FileData * new_file_node = new FileData;
            
            new_file_node->ppuData = newPPUprobSet;
            new_file_node->NumRefs = 0;
            
            fileInstances.insert( map_type::value_type(strFile, new_file_node) );
            
            cout << "Created static instance for '" << strFile << "'"<< endl;
            
            returnVal = new_file_node;
            
        } else {
            //failure
            cerr << "Error loading " << strFile << endl;
            delete newPPUprobSet;
            return NULL;
        }
    } else {
        returnVal = lookupIter->second;
    }
    
    return returnVal;
}

int ppu_load_h5_file(char * filename)
{
    if ( new_ppufile_instance(filename) != NULL ) return 0;
    else return -1;
}

int ppu_create_instance_from_file( char * keyname, char * filename )
{
    
    if ( keyname == NULL || filename == NULL ) return -1;
    
    map_type::iterator lookupIter;
    model_map_type::iterator modelIter;
    
    string strFile=string(filename);
    string strKEY =string(keyname);
    
    modelIter = modelInstances.find(strKEY);
    
    if ( modelIter == modelInstances.end() ) {
        
        
        //An instance of PPUModelData for the keyname not found, so try to create a new instance.
        
        //check first for an HDF5 file of the matching name
        //These are quasi-singleton objects, only one instance per unique HDF5 file name
        lookupIter = fileInstances.find(strFile);
                
        PPUModelData * ppuObj;
        
        FileData * fileDataObject;
        
        fileDataObject = new_ppufile_instance( filename );
        
        if ( fileDataObject != NULL )
        {
            
            ppuObj = fileDataObject->ppuData;
            
            if ( ppuObj == NULL ) {
                //error... could try to reload
                return -1;
            }
            
            ModelData * new_instance_node = new ModelData;
            new_instance_node->h5_shared_data = fileDataObject;
            new_instance_node->P_components=NULL;
            new_instance_node->Q_components=NULL;
            new_instance_node->tempX=NULL;
            new_instance_node->tempY=NULL;
            new_instance_node->tempDE=NULL;
            
            new_instance_node->currentOrder = 0;
            new_instance_node->maxOrderAllocated = 0;
            
            size_t k, N3 = ppuObj->numChannels();
            new_instance_node->tempX = new double[N3+1];
            new_instance_node->tempY = new double[N3+1];
            new_instance_node->tempDE = new double[N3];
            
            for (k=0; k<N3; k++) {
                //node->tempX[k] = 0.5*( node->ppuData->edges[k] + node->ppuData->edges[k+1] );
                
                new_instance_node->tempX[k] = ppuObj->edges[k];
                new_instance_node->tempDE[k] =  ppuObj->edges[k+1] - ppuObj->edges[k];
            }
            new_instance_node->tempX[k] = ppuObj->edges[k];
            
            fileDataObject->NumRefs++;
            
            //insert the object into the instance map so it can be found later
            modelInstances.insert( model_map_type::value_type(strKEY, new_instance_node) );
            
            cout << "Created model instance with KEY '" << strKEY << "'"<< endl;
            
            return 0;
            
        } else {
            //failure
            cerr << "Error creating instance for key  " << strKEY << endl;
            
            return -1;
        }
        
    } 
    
    //if the instance was found, just return OK status (zero)
    return 0;
}

int ppu_delete_instance_data( char * KEY )
{
    ModelData * node;
    FileData * filenode;
    if ( get_instance_node(KEY, &node, false) == 0 )
    {
        free_sized_components( *(node->h5_shared_data->ppuData), node->P_components, node->Q_components, node->maxOrderAllocated );
        delete [] node->tempX;
        delete [] node->tempY;
        delete [] node->tempDE;
        
        node->h5_shared_data->NumRefs--;
        
        if ( node->h5_shared_data->NumRefs == 0 ) {
            string h5file = node->h5_shared_data->ppuData->get_filename();
            delete node->h5_shared_data->ppuData;
            delete node->h5_shared_data;
            
            fileInstances.erase(h5file);
        }
    
        delete node;
        modelInstances.erase(KEY);
        
        return 0;
    } else if ( get_ppufile_node(KEY, &filenode) == 0 )
    {
        if ( filenode->NumRefs != 0 ) {
            cout << "Unable to delete data for file '"<<KEY<<"'"<< endl;
            cout << "Model instances using it still exist." <<endl;
            return 0;
        }
        
        if ( filenode->ppuData != NULL ) {
            string h5file = filenode->ppuData->get_filename();
            fileInstances.erase(h5file);
            
            delete filenode->ppuData;
        }
    }
    
    return -1;
}

int ppu_delete_all()
{
    map_type::iterator fileIter;
    model_map_type::iterator modelIter;
    
    modelIter = modelInstances.begin();
    
    vector<string> keys;
    
    while ( modelIter != modelInstances.end() )
    {
        keys.push_back( modelIter->first );
        modelIter++;
    }
    
    vector<string>::iterator keyItr;
    
    keyItr = keys.begin();
    while (keyItr != keys.end() ) {
        ppu_delete_instance_data((char*)keyItr->c_str());
        keyItr++;
    }
    
    keys.clear();
    
    fileIter = fileInstances.begin();
    
    while (fileIter != fileInstances.end()) {
        keys.push_back( fileIter->first );
        fileIter++;
    }
    
    
    keyItr = keys.begin();
    while (keyItr != keys.end() ) {
        ppu_delete_instance_data((char*)keyItr->c_str());
        keyItr++;
    }
    
    
    return 0;
}

int get_ppufile_node( char * file, FileData ** node)
{
    FileData * stackNode = new_ppufile_instance(file);
    
    *node = stackNode;
    
    if (stackNode==NULL) return -1;
    else {
        return 0;
    }
}

int get_instance_node( char * KEY, ModelData ** node,  bool auto_create)
{
    if ( KEY == NULL ) return -1;
    
    model_map_type::iterator lookupIter;
    
    string strFile=string(KEY);
    
    lookupIter = modelInstances.find(KEY);
    
    if ( lookupIter == modelInstances.end() ) {
        
        //an instance doesn't currently exist
        
        return -1;
        
        /*
        if ( auto_create && ppu_create_instance_from_file(filename) == 0) {
            //worked
            lookupIter = modelInstances.find(KEY);
            *node = lookupIter->second;
        } else {
            //didn't work
            return -1;
        }
         */
        
    } else {
        
        //an instance was found
        
        *node = lookupIter->second;
        
    }
    
    return 0;
}


int get_instance( char * KEY,  PPUModelData ** theInstance )
{
    FileData * fileNode;
    ModelData * hashNode;
    int status = get_instance_node(KEY, &hashNode);
    
    if ( status == 0 ) {
        
        *theInstance = hashNode->h5_shared_data->ppuData;
        return 0;
        
    } else {
        status = get_ppufile_node(KEY, &fileNode);
    }
    
    if ( status != 0 ) return status;
    
    *theInstance = fileNode->ppuData;
    
    return 0;
}

int ppu_get_h5_numchannels( char * filename, long * numchannels )
{
    
    PPUModelData * h5ppuData;
    
    int status = get_instance(filename, &h5ppuData);
    
    if ( status != 0 ) return status;
    
    size_t N3 = h5ppuData->numChannels();
    
    *numchannels = (long)N3;
    
    return 0;
}

int ppu_get_h5_channels( char * filename, double * channelEdges )
{
    
    PPUModelData * h5ppuData;
    
    int status = get_instance(filename, &h5ppuData);
    
    if ( status != 0 ) return status;
    
    size_t N3 = h5ppuData->numChannels();
    
    for (size_t i=0; i<=N3; i++) {
        channelEdges[i] = h5ppuData->edges[i];
    }
    
    return 0;
}

int ppu_set_approximation_order(int order)
{
    ApproximationOrder = order; //default is 8
    return 0;
}

int ppu_set_maximum_order(int order)
{
    return ppu_set_approximation_order(order);
}

int ppu_set_approximation_order(const char * KEY, int order)
{
    ModelData * node;
    int status = get_instance_node((char*)KEY, &node);
    
    if ( status != 0) {
        cerr << "Error: instance not found" << endl;
        return status;
    }
    
    if ( order <= 0 ) {
        cerr << "Error: order must be greater than 0" << endl;
        return -1;
    }
    
    if ( node->maxOrderAllocated > 0 && order > node->maxOrderAllocated ) {
        
        free_sized_components( *(node->h5_shared_data->ppuData), node->P_components, node->Q_components, node->maxOrderAllocated );
        
        allocate_sized_components( *(node->h5_shared_data->ppuData), node->P_components, node->Q_components, order );
        
        node->maxOrderAllocated = order;
        node->currentOrder = order;
        return 0;
    } else {
        node->currentOrder = order;
    }
    
    //node->h5_shared_data->ppuData->SetMaxOrdToCalc(order);
    
    return 0;
}


int ppu_initialize_zeroth_order_pdf( char * KEY, double * inEdges, double * inSpectrum, long numInChans )
{
    
    ModelData * node;
    PPUModelData * ppuObj;
    size_t k, N3;
    int status = get_instance_node(KEY, &node);
    
    if ( status != 0) return status;
    
    ppuObj = node->h5_shared_data->ppuData;
    
    if ( node->P_components == NULL ) {
        //Allocates space to hold the intermediate pulse-state spectral components.
        //Must be freed with free_sized_components(...).  Both in 'AnalyticalPPU.h'
        if ( node->currentOrder <= 0 ) {
            node->currentOrder = ApproximationOrder;
        }
         
        status = allocate_sized_components( *(ppuObj), node->P_components, node->Q_components, node->currentOrder );
    
        if ( status == 0 ) node->maxOrderAllocated = node->currentOrder;
        
    } else {
        //re-initialize to zeros
        status = zero_sized_components( *(ppuObj), node->P_components, node->Q_components, node->maxOrderAllocated );
    }
    
    if ( status != 0) return status;  
    
    
    
    //The 'zeroth-order' term represents the input spectrum.
    //Contrary to DRMs, the energy edges of the PPU operators are generic, and do not match
    //those of the input spectrum.  Thus, 0th order is the input spectrum interpolated to proper channel bins.
    
    //The function below takes an input differential spectrum (or PDF) and produces an array of counts (or probabilities)
    //interpolated into the model channel definitions.  The output is stored as the 0th order spectrum in 'node.P_components'.
    
    N3 = ppuObj->numChannels();
    
    for (k=0; k<=N3; k++) {        
        node->tempY[k] = 0;
    }
    
    status = rebin_pdf_integ(inEdges, inSpectrum, numInChans, node->tempX, node->tempY, N3);
    
    if ( status != 0) return status;
    
    for (k=0; k < N3; k++) node->P_components[0][k] = (float)node->tempY[k];
    
    //for (k=0; k < N3; k++) node->P_components[0][k] = (float)inSpectrum[k];
    
    int prior_state = ppuObj->maxOrdToCalc;
    ppuObj->maxOrdToCalc = node->currentOrder;
    
    status = ppuObj->CalcOrderComponents(node->P_components, node->Q_components);
    
    ppuObj->maxOrdToCalc = prior_state;
    
    return status; //0 if succesfull
}

int ppu_initialize_zeroth_order_counts( char * keyname, double * inEdges, double * inputCountSpectrum, long numInChans )
{
    
    if ( numInChans <= 0 ) return -1;
    
    long i;
    double total=0.0;
    
    double * copy = new double[numInChans];
    
    for (i=0; i<numInChans; i++) {
        total += inputCountSpectrum[i]; 
    }
    
    for (i=0; i<numInChans; i++) {
        copy[i] = inputCountSpectrum[i] / (total*(inEdges[i+1] - inEdges[i])); 
    }
    
    int status = ppu_initialize_zeroth_order_pdf(keyname, inEdges, copy, numInChans);
    
    delete [] copy;
    return status;
}

/**********************************************************
 
 int ppu_compute_distorted_spectrum( char * filename, double cpus, float * answer )
 
 cpus - the input rate in counts per microsecond (cps *10e-06).
 
 answer - array of channel values (already allocated)
 
**********************************************************/
int ppu_compute_distorted_spectrum( char * filename, double cpus, double * outEdges, long outN, double * outCounts)
{
    ModelData * node;
    PPUModelData * ppuObj;
    
    int status = get_instance_node(filename, &node);
    
    if ( status != 0) {
        cerr << "Error: instance not found" << endl;
        return status;
    }
    
    ppuObj = node->h5_shared_data->ppuData;
    
    //The instance corresponding to filename needs to be re-initialized
    if ( node->P_components == NULL ){
        cerr << "Error: instance components are NULL" << endl;
        return -1;
    }
    
    //This returns the number of counts (or probability) per bin into tempY.
    int prior_state = ppuObj->maxOrdToCalc;
    ppuObj->maxOrdToCalc = node->currentOrder;
    
    status = ppuObj->CalcSpectrumExpansion(cpus, node->P_components, node->Q_components, node->tempY );
    
    ppuObj->maxOrdToCalc = prior_state;
    
    if ( status != 0) {
        return status;
    }
    
    //This converts tempY to a differential spectrum or PDF
    size_t N3 = ppuObj->numChannels();
    double * tempPDF = new double[N3+1];
    
    for (size_t k=0; k<N3; k++) {
        
        if ( node->tempY[k] < 0 ) 
            tempPDF[k] = 0;
        else
            tempPDF[k] = node->tempY[k] / node->tempDE[k];
    }
    tempPDF[N3] = tempPDF[N3-1];

    //Now tempPDF[k] represents the average value of the pdf or dn/dE in the interval tempX[k] to tempX[k+1] (i.e., channel 'k').
    
    //This routine integrates the discrete function { tempX, tempPDF } over each channel, given by outEdges.
    //'outCounts' contains the counts or probability per 'outEdges' channel
    
    for (long j=0; j<outN; j++) {
        (outCounts)[j]=0;
    }
    
    rebin_pdf_integ(node->tempX, tempPDF, N3, outEdges, outCounts, outN);
    
    //interpolate0_integ(node->tempX, node->tempY, N3, outEdges, outCounts, outN);
    
    delete [] tempPDF;
    
    return 0;
}



int ppu_compute_distorted_spectrum_modelbasis( char * filename, double cpus, double *& outCounts, double *& outEdges, long * outN)
{
    ModelData * node;
    PPUModelData * ppuObj;
    int status = get_instance_node(filename, &node);
    
    if ( status != 0) {
        cerr << "Error: instance not found" << endl;
        return status;
    }
    
    ppuObj = node->h5_shared_data->ppuData;
    
    //The instance corresponding to filename needs to be re-initialized
    if ( node->P_components == NULL ){
        cerr << "Error: instance components are NULL" << endl;
        return -1;
    }
    
    int prior_state = ppuObj->maxOrdToCalc;
    ppuObj->maxOrdToCalc = node->currentOrder;
    
    //This returns the number of counts (or probability) per bin into tempY.
    status = ppuObj->CalcSpectrumExpansion(cpus, node->P_components, node->Q_components, node->tempY );
    
    ppuObj->maxOrdToCalc = prior_state;
    
    if ( status != 0) return status;
    
    *outN = ppuObj->numChannels();
    
    outEdges = node->tempX;
    outCounts = node->tempY;
    
    /*
    long k;
    for (k=0; k < *outN; k++) {
        outCounts[k] = node->tempY[k];
        outEdges[k] = node->tempX[k];
    }
    outEdges[k] = node->tempX[k];
     */
    return 0;
}


void mc_initialize_chan2en_conversion(double min_keV, double max_keV)
{
    if ( chan2en != NULL ) delete chan2en;
    
    chan2en = new ChannelEnergy(min_keV, max_keV);
    
};


long mc_simulate_ppu_spectrum( int throwMode, long * mc_size_ptr, double * exposure_us_ptr, double rate_cpus, double * inSpectrum,  double * outSpectrum, double * inEdges, long numInChans, float ** counted_data, gsl_rng ** existingRNG )
{
    double exposure_us;
    long mc_size;
    const gsl_rng_type *generatorType;
    gsl_rng * gslRNG;
    
    if ( existingRNG == NULL ) {
     
        generatorType = gsl_rng_ran2;
        gslRNG = gsl_rng_alloc(generatorType);
        gsl_rng_set (gslRNG, (unsigned)time(0)); //state
        
    } else if ( *existingRNG == NULL ) {
        
        generatorType = gsl_rng_ran2;
        gslRNG = gsl_rng_alloc(generatorType);
        gsl_rng_set (gslRNG, (unsigned)time(0)); //state
        
        *existingRNG = gslRNG;
        
    } else {
        
        gslRNG = *existingRNG;
        
    }
    
    
    double energy;
    double erMin = inEdges[0];
    double erMax = inEdges[numInChans];
    double erRange = erMax - erMin;
    
    long j,k, trial_number, maxNumEvents=1;
    long num_counted;
    long detectionQuota;
    
    double trial_exposure, tsample;
    
    if ( throwMode==PPU_MC_FIXED_EXPOSURE ) 
    {
        exposure_us = *exposure_us_ptr;
        maxNumEvents = long(exposure_us*rate_cpus + 6.0 * sqrt(exposure_us*rate_cpus)); //allocate enough for 6-sigma event
    } else if ( throwMode==PPU_MC_FIXED_NUMBER )
    {
        mc_size = *mc_size_ptr;
        maxNumEvents = mc_size;
    }
    else if ( throwMode==PPU_MC_FIXED_DETECTION_NUMBER ) {
        
        double inputFactor;
        
        if ( rate_cpus < 0.3 ) inputFactor = 5.0 / ( 1.0 - rate_cpus * 2.6 );
        else {
            inputFactor = 6.0 / ( 1.0 - 0.3 * 2.6 );
        }
        
        detectionQuota =  *mc_size_ptr;
        
        maxNumEvents = inputFactor*detectionQuota;
        
    }
    else return -1;
    
    cout << "MaxNumEvents = " << maxNumEvents << endl;
    
    
    double specMax=0;
    for (k=0;k<numInChans;k++) {
        if ( inSpectrum[k] > specMax ) specMax = 1.01*inSpectrum[k];
        outSpectrum[k] = 0.0;
    }
    
    float * inputTime =  new float[maxNumEvents];
    float * inputVolts = new float[maxNumEvents];
    
    float * countData = new float[2*maxNumEvents];
    
    float * countTimes = &(countData[0]);
    float * countVolts = &(countData[maxNumEvents]);
    
    trial_number=0;
    
    while ((throwMode==PPU_MC_FIXED_EXPOSURE && trial_exposure < exposure_us && trial_number < maxNumEvents) || 
           (throwMode==PPU_MC_FIXED_NUMBER && trial_number < mc_size) ||
           (throwMode==PPU_MC_FIXED_DETECTION_NUMBER && trial_number < maxNumEvents)) {
        
        tsample = gsl_ran_exponential( gslRNG, 1.0 / rate_cpus );
        
        inputTime[trial_number] = 1e-06*(tsample + trial_exposure);
        
        trial_exposure += tsample;
        
        trial_number++;
        
#ifdef DECORRELATE_NUM
        for(j=0;j<DECORRELATE_NUM;j++) gsl_rng_uniform(gslRNG);
#endif
    }
    
    exposure_us = trial_exposure;
    
    cout << "exposure_us = " << exposure_us << endl;
    
    long chanSample;
    k=0;
    while (k < trial_number) {
        
        chanSample = floor( (double)numInChans * gsl_rng_uniform(gslRNG) );
#ifdef DECORRELATE_NUM
        for(j=0;j<DECORRELATE_NUM;j++) gsl_rng_uniform(gslRNG);
#endif
        tsample = gsl_rng_uniform(gslRNG);
        
        while (1) {
            
            if ( tsample <= ( inSpectrum[chanSample] /  specMax ) )
            {
                erMin = inEdges[chanSample];
                erRange = inEdges[chanSample+1] - erMin;
                
                energy = erMin + erRange * gsl_rng_uniform(gslRNG);
                
                inputVolts[k] = chan2en->EtoV(energy);
                
                //cout << energy << " - > " << inputVolts[k] << endl;
                break;
            } else {
                chanSample = floor( (double)numInChans * gsl_rng_uniform(gslRNG) );
#ifdef DECORRELATE_NUM
                for(j=0;j<DECORRELATE_NUM;j++) gsl_rng_uniform(gslRNG);
#endif
                tsample = gsl_rng_uniform(gslRNG);
            }
            
        }
        
        k++;
    }
    
    mc_size = k;
    
    if ( throwMode == PPU_MC_FIXED_DETECTION_NUMBER ) {
        
        size_t numNeeded;
        
        num_counted = (long)ProcessPPUtoQuota(detectionQuota, numNeeded, inputTime, inputVolts, (size_t)mc_size, countTimes, countVolts );
        
        mc_size = numNeeded;
        
        if ( num_counted > 0 ) {
            exposure_us = 1e06*countTimes[num_counted-1];
        }
        
    } else {
        num_counted = (long)ProcessPulsePileup( inputTime, inputVolts, (size_t)mc_size, countTimes, countVolts );
    }
    
    cout << "num_counted = " << num_counted << endl;
    
    j=-1;
    for (k=0; k<num_counted; k++) {
        
        energy = chan2en->VtoE(countVolts[k]);
        
        j = channel_search(inEdges, energy, numInChans);
        
        if ( j >= 0 ) {
            (outSpectrum)[j]++;
        }
        
    }
    
    
    
    delete [] inputTime;
    delete [] inputVolts;
    
    if ( counted_data != NULL ) {
        
        *counted_data = new float [(num_counted+1)*2];
        
        for (k=0; k<num_counted; k++) {
            (*counted_data)[k] = countTimes[k];
            (*counted_data)[k+num_counted] = countVolts[k];
            
            
            
        }
        
    }
    else delete [] countData;
    
    *exposure_us_ptr = exposure_us;
    *mc_size_ptr = mc_size;
    
    if ( existingRNG == NULL ) {
        gsl_rng_free(gslRNG);
    }
    
    return num_counted;
    
}


