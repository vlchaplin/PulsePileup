//
//  ppu_IDL_interface.cpp
//  C_ppu
//
//  Created by vchaplin on 7/30/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "ppu_static_interface.h"
#include "stddef.h"

//Included from $IDL_DIR/external/
#include "idl_export.h"

extern "C" {
    
    IDL_VPTR PPUSetMaxOrder(int argc, IDL_VPTR argv[] );
    IDL_VPTR PPULoadH5File(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUCreateInstance(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUSetApproxOrder(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUGetH5NumChannels(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUGetH5Channels(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUInitZerothOrderCounts(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUComputeSpectrum(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUDeleteInstance(int argc, IDL_VPTR argv[]);
    IDL_VPTR PPUDeleteAll(int argc, IDL_VPTR argv[]);
    
};

char * resolve_filename(IDL_STRING * idl_filename)
{
    if ( idl_filename->slen == 0 || idl_filename->s == NULL ) {
        return NULL;
    }
    
    //static allocation to hold the contents of the IDL string
    static char filename_c[600];
    size_t slen = (size_t)(idl_filename->slen);
    
    strncpy(filename_c, idl_filename->s, slen );
    filename_c[slen] = '\0';
    
    return filename_c;
}

char * resolve_modelname(IDL_STRING * idl_modelname)
{
    if ( idl_modelname->slen == 0 || idl_modelname->s == NULL ) {
        return NULL;
    }
    
    //static allocation to hold the contents of the IDL string
    static char modelname_c[600];
    size_t slen = (size_t)(idl_modelname->slen);
    
    strncpy(modelname_c, idl_modelname->s, slen );
    modelname_c[slen] = '\0';
    
    return modelname_c;
}



//Set the default maximum order
IDL_VPTR PPUSetMaxOrder(int argc, IDL_VPTR argv[])
{
    int maximumOrder;
    if ( argv[0]->type == IDL_TYP_INT ) {
        maximumOrder = (int)argv[0]->value.i;
    }
    else if ( argv[0]->type == IDL_TYP_LONG ) {
        maximumOrder = (int)argv[0]->value.l;
    }
    else {
        maximumOrder=-1;
    }
    
    //std::cout << maximumOrder << std::endl;
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    
    if ( maximumOrder <= 0 ) 
        var->value.l = -1;
    else 
        var->value.l = ppu_set_maximum_order(maximumOrder);
    
    return var;
};



//Load up to MaxOrder of the hdf5 pileup response file.
//If MaxOrder is greater than the number of orders in the file, 
//the whole file is read.
IDL_VPTR PPULoadH5File(int argc, IDL_VPTR argv[])
{
    IDL_STRING stringArg = argv[0]->value.str;
    
    char * filename_c = resolve_filename( &stringArg );    
    
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    var->value.l = ppu_load_h5_file(filename_c);
    
    return var;
};

IDL_VPTR PPUCreateInstance(int argc, IDL_VPTR argv[])
{
    char * instance_c = resolve_modelname( &(argv[0]->value.str) );
    char * filename_c = resolve_filename(  &(argv[1]->value.str) );   
    
    int status = ppu_create_instance_from_file(instance_c, filename_c);
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    var->value.l = (IDL_LONG)status;
    
    return var;
};

IDL_VPTR PPUSetApproxOrder(int argc, IDL_VPTR argv[])
{
    char * instance_c = resolve_modelname( &(argv[0]->value.str) );
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    
    if ( instance_c == NULL ) {
        cerr << "Argument 1 doesn't match any model instances.  Use PPUCreateInstance(name, h5file) to define a new instance."<<endl;
        var->value.l=-1;
        return var;
    }
    
    int order;
    if ( argv[1]->type == IDL_TYP_INT ) {
        order = (int)argv[1]->value.i;
    }
    else if ( argv[1]->type == IDL_TYP_LONG ) {
        order = (int)argv[1]->value.l;
    }
    else {
        order=-1;
        cerr << "Argument 2 must be an INT or LONG > 0" << endl;
        var->value.l=-1;
    }
    
    int status = ppu_set_approximation_order(instance_c, order);
    
    

    var->value.l = (IDL_LONG)status;
    
    return var;
}

IDL_VPTR PPUGetH5NumChannels(int argc, IDL_VPTR argv[])
{
    char * filename_c = resolve_filename( &(argv[0]->value.str) ); 
    
    long nchan=0;
    
    int status = ppu_get_h5_numchannels(filename_c, &nchan);
    
    if ( argv[1]->type == IDL_TYP_INT ) {
        argv[1]->value.i = (IDL_INT)nchan;
    }
    else if ( argv[1]->type == IDL_TYP_LONG ) {
        argv[1]->value.l = (IDL_LONG)nchan;
    }
    else {
        cerr << "Must be type INT or LONG (signed)" << endl;
        status = -1;
    }
    
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    var->value.l = (IDL_LONG)status;
    
    return var;
};


IDL_VPTR PPUGetH5Channels(int argc, IDL_VPTR argv[])
{
    char * filename_c = resolve_filename( &(argv[0]->value.str) );
    
    IDL_VARIABLE * returnData;
    long nchan=0;
    double * edgearray;
    
    int status = ppu_get_h5_numchannels(filename_c, &nchan);
    
    if ( status == 0 ){
        
        IDL_LONG64 dim[1];
        
        dim[0] = (IDL_LONG64)(nchan+1);
        
        edgearray = (double*)IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, dim, IDL_BARR_INI_ZERO, &returnData);
        
        status = ppu_get_h5_channels(filename_c, edgearray);
        
        IDL_VarCopy(returnData, argv[1]);
    }
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    var->value.l = (IDL_LONG)status;
    
    return var;
    
};

IDL_VPTR PPUInitZerothOrderCounts(int argc, IDL_VPTR argv[])
{
    char * instance_c = resolve_modelname( &(argv[0]->value.str) );
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    
    if ( instance_c == NULL ) {
        cerr << "Argument 1 doesn't match any model instances.  Use PPUCreateInstance(name, h5file) to define a new instance."<<endl;
        var->value.l=-1;
        return var;
    }
    
    if ( argv[1]->type != IDL_TYP_DOUBLE ) {
        cerr << "Argument 2 must be an array of type DOUBLE giving the channel edges of the input spectrum" << endl;
        var->value.l = -1;
        return var;
    }
    if ( argv[2]->type != IDL_TYP_DOUBLE ) {
        cerr << "Argument 3 must be an array of type DOUBLE giving the input count spectrum" << endl;
        var->value.l = -1;
        return var;
    }
    
    long nchan;
    
    /*
    if ( argv[3]->type == IDL_TYP_INT ) {
        nchan = (int)argv[3]->value.i;
    }
    else if ( argv[3]->type == IDL_TYP_LONG ) {
        nchan = (int)argv[3]->value.l;
    }
    else {
        nchan=-1;
        cerr << "Argument 4 must be an INT or LONG giving the number of input channels" << endl;
        var->value.l=-1;
    }
     */
    
    double * edgearray = (double*)(argv[1]->value.arr->data);
    double * countsarray = (double*)(argv[2]->value.arr->data);
    
    IDL_LONG64 * edgeDim = argv[1]->value.arr->dim;
    IDL_LONG64 * chanDim = argv[2]->value.arr->dim;
        
    if ( edgeDim[0] < 2 || chanDim[0] >= edgeDim[0] ) {
        cerr << "Error, edge array (arg 2) should have exactly one more element than the counts array (arg 3)" <<endl;
        var->value.l = -1;
        return var;
    }
    
    nchan = chanDim[0];
    argv[3]->type = IDL_TYP_LONG;
    argv[3]->value.l = nchan;
    
    int status = ppu_initialize_zeroth_order_counts(instance_c, edgearray, countsarray, nchan);
    
    var->value.l = (IDL_LONG)status;
    
    return var;
};



IDL_VPTR PPUComputeSpectrum(int argc, IDL_VPTR argv[])
{
    char * instance_c = resolve_modelname( &(argv[0]->value.str) );
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    
    if ( instance_c == NULL ) {
        cerr << "Argument 1 doesn't match any model instances.  Use PPUCreateInstance(name, h5file) to define a new instance."<<endl;
        var->value.l=-1;
        return var;
    }
    if ( argv[1]->type != IDL_TYP_DOUBLE ) {
        cerr << "Argument 2 must of type DOUBLE and give the input rate in counts per microsecond" << endl;
        var->value.l = -1;
        return var;
    }
    if ( argv[2]->type != IDL_TYP_DOUBLE ) {
        cerr << "Argument 3 must be an array of type double giving the desired channel edges of the output spectrum" << endl;
        var->value.l = -1;
        return var;
    }
    
    long nchan;
    if ( argv[3]->type == IDL_TYP_INT ) {
        nchan = (int)argv[3]->value.i;
    }
    else if ( argv[3]->type == IDL_TYP_LONG ) {
        nchan = (int)argv[3]->value.l;
    }
    else {
        nchan=-1;
        cerr << "Argument 4 must be an INT or LONG giving the number of output channels" << endl;
        var->value.l=-1;
    }
    
    IDL_VARIABLE * returnData;
    double cpus = argv[1]->value.d;
    double * outEdges = (double*)argv[2]->value.arr->data;
    double * outSpec;
    
    IDL_LONG64 * dim = argv[2]->value.arr->dim;
    IDL_LONG64 chdim[1];
    
    nchan = dim[0] - 1;
    chdim[0] = dim[0] - 1;
    
    outSpec = (double*)IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, chdim, IDL_BARR_INI_ZERO, &returnData);
    
    int status = ppu_compute_distorted_spectrum(instance_c, cpus, outEdges, nchan, outSpec);
    
    outSpec[0] = 1.0;
    
    IDL_VarCopy(returnData, argv[4]);
    
    var->value.l = (IDL_LONG)status;
    
    return var;
    
};

IDL_VPTR PPUDeleteInstance(int argc, IDL_VPTR argv[])
{
    char * instance_c = resolve_modelname( &(argv[0]->value.str) );
    int status = ppu_delete_instance_data(instance_c);
    
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    var->value.l = (IDL_LONG)status;
    
    return var;
}

IDL_VPTR PPUDeleteAll(int argc, IDL_VPTR argv[])
{
    int status = ppu_delete_all();
    
    
    IDL_VARIABLE *var;
    var = (IDL_VARIABLE*) IDL_Gettmp();
    var->type = IDL_TYP_LONG;
    var->value.l = (IDL_LONG)status;
    
    return var;
}







