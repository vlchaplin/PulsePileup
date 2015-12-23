//
//  ppu_pha_spec.cpp
//  C_ppu
//
//  Created by vchaplin on 7/12/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iomanip>
#include <iostream>
#include <vector>

#include "PHA_IO.hh"
#include "ppu_static_interface.h"
#include "spoccExeUtilities.h"

#define LOWER_VOLTS 0.005
#define UPPER_VOLTS 5.0

//DEFAULT_H5_PPU is defined in Makefile.am

void _print_usage() {
    
    cout << endl;
    cout << "  Usage: ppu_pha_spec -i input_spec.pha [-o output_spec.pha] [-Rs rate_cps -Rus rate_cpus] [-h5 PPUfile.h5] " << endl;
    cout << endl;
    cout << NewArgument("", "Calculate the PPU distortion for the input spectra, which are passed as PHA files.  The rate is passed as an argument.  The output is written as PHA file. Input PHAs should have EBOUNDS in units of keV.", 10, 100) << endl << endl;
    cout << NewArgument("-i", "input_spec.pha ....   Input the name(s) of single-pha files (Type 1) containing the input spectral shape. Multiple input shapes will be processed with the same rate.", 20, 100) << endl;
    cout << NewArgument("-o", "out_spec.pha ....   Input the name(s) files to contain the output shape. If Multiple input shapes are given, the same number of outputs should be iven.  If the same number is given the user will be prompted for the filename.", 20, 100) << endl;
    cout << NewArgument("-Rs", "Input rate in counts per second.  Only one of -Rs or -Rus should be used.", 20, 100) << endl;
    cout << NewArgument("-Rus", "Input rate in counts per microsecond.  Only one of -Rs or -Rus should be used.", 20, 100) << endl;
    cout << NewArgument("-h5", "Name of the HDF5 'pulse-pileup response file' to use.  If abset the program will use the file named by the environment variable H5_PPU_FILE. If still that is not defined, the default file tried is: '"+string(DEFAULT_H5_PPU)+"'", 20, 100) << endl;
    cout << NewArgument("-fakenorm", "Use the input rate as given, but renormalize the output spectrum to 'match' the normalization of the input spectrum.  Normally, only the shape is taken from the input file, and the -R argument determines the output normalization. Using this argument, the output PPU probability-spectrum is multiplied by the total of the input spectrum, allowing a direct shape comparison without having to re-normalize the input as well.  Note that this results in a 'non-physical' output rate, because it doesn't correspond to the assumed input rate (-R..) which generated it.", 20, 100) << endl;
    cout << endl;
    
}

int main (int argc, char * const argv[]) {

    char ** argitr, ** argend;
    bool slurp;
    string arg;
    
    double rate_microsec = -1.0;
    
    vector<string> inputFiles;
    vector<string> outputFiles;
    
    char h5_ppu_file[500];
    bool h5_is_in_volts=false;
    bool use_fake_normalization=false;
    
    char * env = getenv("H5_PPU_FILE");
    if ( env != NULL )
        strcpy(env, h5_ppu_file);
    else
        strcpy(h5_ppu_file, DEFAULT_H5_PPU);
    
    argitr = (char **)argv;
    argend = (char **)(argv + argc);
    while ( argitr != argend )
    {
        if ( (**argitr) == '-' ) {
            
            slurp=0;
            arg=string(*argitr);
        } 
        
        if ( arg == "-i" ) {
            if (slurp)
                inputFiles.push_back( *argitr );
            else 
                slurp=1;
        }
        else if ( arg == "-o" ) {
            if (slurp)
                outputFiles.push_back( *argitr );
            else 
                slurp=1;
        } 
        else if ( arg == "--help" ) {
            _print_usage();
            return 0;
        }
        else if ( arg == "-Rs" ) {
            if (slurp)
                rate_microsec = 1e-06 * atof( *argitr );
            else 
                slurp=1;
        }
        else if ( arg == "-Rus" ) {
            if (slurp)
                rate_microsec = atof( *argitr );
            else 
                slurp=1;
        }
        else if ( arg == "-h5" ) {
            if (slurp)
                strcpy(h5_ppu_file,  *argitr );
            else 
                slurp=1;
        }
        else if ( arg == "--fakenorm" ) {
            use_fake_normalization=true;
        }
        
        
        argitr++;
    }
    
    if ( inputFiles.size() == 0 ) {
        cerr << "Input a PHA spectrum using -i" << endl;
        return 1;
    }
    if ( rate_microsec <= 0 ) {
        cerr << "Specify the rate (-Rs for count per sec, -Rus for counts per microsecond)" << endl;
        return 1;
    }
    
    size_t n;
    
    //do a quick calculation to determine approximation order
    double mu = rate_microsec* (6.0);
    int order = 0;
    double P = exp(-mu);
    
    while ( P < 0.99  ) {
        order++;
        P += pow(mu, order)*exp(-mu) / gsl_sf_fact(order);
    }
    
    if ( order == 0 ) order = 1;
    
    cout << endl << "Approximation order will be: " << order << endl << endl;
    
    if ( order > 15 ) {
        cout << "Warning: the order for 99th percentile (99% of all events have PPU order <= this) is very high. The model has been validated only up to rates ~ 1 MHz." << endl << endl;;
        order = 15;
    }
        
    ppu_set_maximum_order(order);
    
    if ( ppu_load_h5_file(h5_ppu_file) ) {
        cerr << "Error loading the pulse-pileup response (the HDF5 file)" << endl;
        return 1;
    }
    
    double * h5_basis_edges = NULL;
    long h5_basis_size=0;
    
    if ( ppu_get_h5_numchannels(h5_ppu_file, &h5_basis_size) || h5_basis_size <= 0 )
    {
        cerr << "Error reading the channel edges of the pulse-pileup response (the HDF5 file)" << endl;
        return 1;
    }
    
    h5_basis_edges = new double[h5_basis_size+1];
    
    if ( ppu_get_h5_channels(h5_ppu_file, h5_basis_edges) ) {
        cerr << "Error reading the channel edges of the pulse-pileup response (the HDF5 file)" << endl;
        return 1;
    }
    
    if ( h5_basis_edges[0] <= 0.03 && h5_basis_edges[h5_basis_size] <= 7.0 ) {
        cout << "Voltage basis automatically detected. Input PHA channel energies will be mapped linearly into "<<LOWER_VOLTS<<" -> "<<UPPER_VOLTS<<" volts" << endl << endl;;
        h5_is_in_volts=true;
    }
    
    
    if ( inputFiles.size() > 1 ) {
        
        //creates a dummy instance.  This is the h5 data is not deleted at the end of the first loop. 
        //Normally, when a named instance is deleted using ppu_delete_instance_data(), 
        //the number of references to the shared file data is decremented. If the ref. count
        //to this file is 1, deleting the instance also deletes the shared data.
        
        ppu_create_instance_from_file((char*)"dummy", h5_ppu_file);
    }
    
    for (n=0; n<inputFiles.size(); n++) {
        
        //a custom object for reading FITS with binary tables
        FxbReadStat fxb;
        
        fxb.open( (char*) inputFiles[n].c_str() );
        
        fxb.move2hdu(2);
        
        if ( fxb.status ) {
            fxb.errclose();
            continue;
        }
        
        long i, nedges = fxb.n_rows+1;
        double * edges = new double[fxb.n_rows+1];
        double * edges_volts;
        
        //read the edges
        fxb.read_col(3, 1, fxb.n_rows, 1+edges);
        fxb.read_col(2, 1, 1, edges);
        
        fxb.move2hdu(3);
        
        if ( fxb.status ) {
            fxb.errclose();
            delete [] edges;
            continue;
        }
        
        edges_volts = new double[fxb.n_rows+1];
        
        // ASSUMPTION regarding ENERGY-TO-VOLTS  
        // edges[0] --> LOWER_VOLTS V.  edges[n-2] --> UPPER_VOLTS V 
        // Intermediate energies are mapped linearly into voltage
        for (i=0; i<nedges; i++) {
            edges_volts[i] = (UPPER_VOLTS - LOWER_VOLTS)*( (edges[i] - edges[0]) / (edges[nedges-2]- edges[0]) ) + LOWER_VOLTS;
        }
        
        
        double * spectrum = new double[fxb.n_rows];
        double * out_spectrum = new double[fxb.n_rows];
        
        //read the spectrum
        fxb.read_col(2, 1, fxb.n_rows, spectrum);
        
        if ( fxb.status ) {
            fxb.errclose();
            delete [] spectrum;
            delete [] edges_volts;
            delete [] out_spectrum;
            continue;
        }
        
        double total=0.0;
        
        long nchan = fxb.n_rows;
        
        //calculate the input normalization
        for (i=0; i<nchan; i++) {
            total += spectrum[i];
        }
         
        const char * tagname = inputFiles[n].c_str();
        
        double * edges_we_use;
        
        if ( h5_is_in_volts ) {
            edges_we_use = edges_volts;
        } else {
            edges_we_use = edges;
        }
        
        
        ppu_create_instance_from_file( (char*)tagname, h5_ppu_file);
        
        ppu_set_approximation_order((char*)tagname, order);
        
        ppu_initialize_zeroth_order_counts( (char*)tagname, edges_we_use, spectrum, nchan);
                
        //get the output spectrum
        ppu_compute_distorted_spectrum((char*)tagname, rate_microsec, edges_we_use, nchan, out_spectrum);
        
        
        
        cout << "Change in spectrum by channel (%): " << endl;
        for (i=0; i<nchan; i++) {
            
            if ( use_fake_normalization ) {
                out_spectrum[i] *= total;
            } else {
                //convert to cps (a standard unit for PHA files)
                out_spectrum[i] *= rate_microsec*1e06;
                
                spectrum[i] *= (rate_microsec*1e06/total);
            }

            
            if ( spectrum[i] > 0 )
                cout << setw(8)<<i<<"    " << setprecision(3) << 100.0*((out_spectrum[i]-spectrum[i]) / spectrum[i]  ) << endl;
            else {
                cout << setw(8)<<i<<"    "<< 0 << endl;;
            }
        }
        
        cout << endl;
        
        string outfile;
        
        
        
        if ( n < outputFiles.size() ) {
            outfile = outputFiles[n];
        } else {
            bool ok = false;
            string userInput;
            while (! ok) {
                cout << "Type the output PHA file (relative or absolute path, 'q' cancels):" << endl;
                getline( cin, userInput );
                
                if ( userInput == "q" ) {
                    outfile = "";
                    ok=true;
                } else if (userInput.size()>0) {
                    outfile = userInput;
                    ok=true;
                }
                
            }
        }
        
        if ( outfile.length() > 0 ) {
         
            cout << "Writing output spectrum: " << outfile << endl;
            
            PHAWriter phawriter;
            
            SinglePHA<double, double, double> outphaObj;
            
            EdgeSet<double> * ebounds = new EdgeSet<double>;
            
            ebounds->newFromArray(edges, nchan+1);
            
            outphaObj.UseBinning(ebounds);
            
            outphaObj.IsRates(1);
            outphaObj.UsePoissonErrors(1);
            
            outphaObj.copySpectrum(out_spectrum, nchan);
            
            phawriter.setFile(outfile);
            phawriter.setTemplate(inputFiles[n]);
            
            phawriter.WriteDataFile(outphaObj);
            
            if ( phawriter.status() == 0 ) {
             
                char * in_spec = (char*)inputFiles[n].c_str();
                double cps = rate_microsec*1e06;
                phawriter.open(0,READWRITE);
                phawriter.write_key( (char*)"PPUH5", TSTRING, h5_ppu_file, -1, (char*)"HDF5 pileup response file");
                phawriter.write_key( (char*)"PPUFILE", TSTRING, in_spec, -1, (char*)"File with input spectrum");
                phawriter.write_key( (char*)"PPUCPS", TDOUBLE, &cps, -1, (char*)"Input rate to PPU model in counts per sec");
                phawriter.write_key( (char*)"PPUORDER", TINT, &order,-1, (char*)"Approximation order of PPU model");
                phawriter.up_checksum_date();
                
                phawriter.move2hdu(2);
                phawriter.write_key( (char*)"PPUH5", TSTRING, h5_ppu_file, -1, (char*)"HDF5 pileup response file");
                phawriter.write_key( (char*)"PPUFILE", TSTRING, in_spec, -1, (char*)"File with input spectrum");
                phawriter.write_key( (char*)"PPUCPS", TDOUBLE, &cps, -1, (char*)"Input rate to PPU model in counts per sec");
                phawriter.write_key( (char*)"PPUORDER", TINT, &order,-1, (char*)"Approximation order of PPU model");
                phawriter.up_checksum_date();
                
                
                phawriter.move2hdu(3);
                phawriter.write_key( (char*)"PPUH5", TSTRING, h5_ppu_file, -1, (char*)"HDF5 pileup response file");
                phawriter.write_key( (char*)"PPUFILE", TSTRING, in_spec, -1, (char*)"File with input spectrum");
                phawriter.write_key( (char*)"PPUCPS", TDOUBLE, &cps, -1, (char*)"Input rate to PPU model in counts per sec");
                phawriter.write_key( (char*)"PPUORDER", TINT, &order,-1, (char*)"Approximation order of PPU model");
                phawriter.up_checksum_date();
                
                phawriter.close();
                
            }
            
        }
        
        
        
        ppu_delete_instance_data((char*)tagname);

        
        delete [] edges;
        delete [] edges_volts;
        delete [] spectrum;
        delete [] out_spectrum;
        
    }
    
    
    

}