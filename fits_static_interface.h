//
//  fits_static_interface.h
//  C_ppu
//
//  Created by vchaplin on 5/7/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef C_ppu_fits_static_interface_h
#define C_ppu_fits_static_interface_h

#include "PHA_IO.hh"
#include <string>
#include <stdlib.h>
#include <iostream>
#include <map>

using namespace std;



// this is just a very c++'ey kind of hash. maps filenames to memory objects that will hold their data.
// Fast access, ~O(log N).
typedef map<string, FxbReadStat *> map_type;
static map_type fileInstances; 

inline int fits_create_instance_from_file( const char * filename, int fits_mode = READWRITE )
{
    map_type::iterator lookupIter;
    
    string strFile=string(filename);
    
    lookupIter = fileInstances.find(strFile);
    
    if ( lookupIter == fileInstances.end() ) {
        
        FxbReadStat * newData = new FxbReadStat;
        
        if ( newData->open( (char*)filename, 0, fits_mode ) != 0 ) {
            //failure. there is a problem reading the file
            int rc = newData->status;
            delete newData;
            return rc;
        } else {
            //success

            //insert the object into the instance map so it can be found later
            fileInstances.insert( map_type::value_type(strFile, newData) );
            
            cout << "Created static instance for '" << strFile << "'"<< endl;
            return 0;
        }   
    }
    
    return 0;
}

inline int fits_get_instance_node( const char * filename, FxbReadStat ** node, bool auto_create=true) 
{
    map_type::iterator lookupIter;
    
    string strFile=string(filename);
    
    lookupIter = fileInstances.find(strFile);
    
    if ( lookupIter == fileInstances.end() ) {
        
        if ( auto_create && fits_create_instance_from_file( filename, READWRITE ) == 0 ) {
            
            lookupIter = fileInstances.find(strFile);
            
            *node = lookupIter->second;
            
            return 0;
            
        } else {
            
            *node = NULL;
            
            return -1;
        }
        
    }
    
    *node = lookupIter->second;
    
    return 0;

};


inline int fits_delete_instance_data( const char * filename ) {

    FxbReadStat * node=NULL;
    if ( fits_get_instance_node(filename, &node, false) == 0 && node != NULL )
    {
        node->flush();
        node->close();
        delete node;
        
        fileInstances.erase(filename);
        
        return 0;
    }
    
    return -1;
    
    
};


inline int fits_move2hdu( const char * filename, const char * hduname ) {
    FxbReadStat * node=NULL;
    if ( fits_get_instance_node(filename, &node, false) == 0 && node != NULL)
    {
        int rc= node->move2hdu((char*)hduname);
        
        node->status=0;
        
        return rc;
        
    }
    
    return -1;
};

inline int fits_move2hdu( const char * filename, int hdunum ) {
    FxbReadStat * node=NULL;
    if ( fits_get_instance_node(filename, &node, false) == 0 && node != NULL)
    {
        int rc =node->move2hdu(hdunum);
        
        node->status=0;
        
        return rc;
    }
    
    return -1;
};


#endif
