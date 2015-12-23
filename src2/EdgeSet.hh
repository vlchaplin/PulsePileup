/*
 *  EdgeSet.hh
 *  
 *
 *  Created by Vandiver L. Chapin on 4/5/10.
 *  Copyright 2010 The University of Alabama in Huntsville. All rights reserved.
 *
 */
 
#ifndef EDGESET_HH
#define EDGESET_HH

#ifndef EDGEDEBUG
 #define EDGEDEBUG 0
#endif

#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "DynMatrix.h"

using namespace std;



template<typename T>
class EdgeSet : public DynMatrix<T> {
	private:
	int lookup_method;
    T bounds[2];
    double binFactor;

	public:
	
	typedef typename DynMatrix<T>::size_type size_type;
	
	EdgeSet() : DynMatrix<T>(2,1) {
		lookup_method = 0;
	};
	EdgeSet(int channels, T nvalue) : DynMatrix<T>(2,1) {
		reinit( channels, nvalue);
        //MakeLinearBins( (T)10.0, (T)1000.0, channels );
	};
	EdgeSet(T min, T max, int channels) : DynMatrix<T>(2,1) {
        MakeLinearBins( min, max, channels );
	};
	~EdgeSet() {
	
	};
    
    void newFromArray(T * edge_array, size_type nedges)
    {
        this->reinit(nedges-1, 0);
        
        for (size_type i=0; i<nedges-1; i++) {
            (*this)[0][i] = edge_array[i];
            (*this)[1][i] = edge_array[i+1];
        }
    }
	
	void reinit(long nchan, T blankval) {
		this->DynMatrix<T>::reinit(2, nchan,blankval);
	};
	
	void appendBin(T lowerEdge, T upperEdge) {
		
		if ( this->size() < 1 ) {
			this->DynMatrix<T>::resize(2,1);
			(*this)[0][0] = lowerEdge;
			(*this)[1][0] = upperEdge;
		}
		else {
			vector<T> vec;
			vec.push_back(lowerEdge);
			vec.push_back(upperEdge);
			
			this->addColumn(vec);
		}
		
		this->calcBounds();
	};
	
	void calcBounds() {
		bounds[0] = (*this)[0][0];
		bounds[1] = (*this)[1][this->size()-1];
	};
	
	virtual void getBounds(T * minmax) {
		minmax[0] = bounds[0];
		minmax[1] = bounds[1];
	};
    
    void MakeLinearBins (T min, T max, size_t channels) {
        this->DynMatrix<T>::reinit(2, channels,(T) 0);
		
        T t_channelWidth = (max - min) / T(channels);
        T edge_iter = min;
		
        typename DynMatrix<T>::size_type b;
        b = 0;
		
        if (t_channelWidth == 0) {
			cout << "ERROR EdgeSet<>::MakLinBins() : Channel-width is zero with given parameters\n";
			return;
		}
        
        bounds[0] = min;
        while (b < channels) {
            (*this)[0][b] = edge_iter;
			//cout << edge_iter << ", " << t_channelWidth << ", ";
			edge_iter += t_channelWidth;
			//cout << edge_iter << "\n";
            (*this)[1][b] = edge_iter;
			//cout << (*this)[0][0] << "," << (*this)[1][0] << "\n";
			b++;
        }
        bounds[1] = edge_iter;
        //cout << "---1c---\n"<< *this;
        lookup_method = 1;
    };
    void MakeLogBins (T min, T max, size_type channels) {
		this->DynMatrix<T>::reinit(2, channels,(T) 0);
        if (min == 0.0) min += 1;
        double f = pow( (double)(max / min), 1.0f / double( channels ) );
        T edge_iter = min;
        
		if (f == 0) {
			cout << "ERROR EdgeSet<>::MakLogBins() : Channel-width is zero with given parameters\n";
			return;
		}
		
		typename DynMatrix<T>::size_type b;
        b = 0;

        bounds[0] = min;
        while (b < channels) {
            (*this)[0][b] = edge_iter;
            edge_iter *= f;
            (*this)[1][b] = edge_iter;
            b++;
        }
        bounds[1] = edge_iter;
        
        lookup_method = 2;
        this->binFactor = f;
    };
	
	static inline EdgeSet<T> logbingen(T min, T max, size_type nbins) {
		
		EdgeSet<T> temp;
		temp.MakeLogBins( min, max, nbins );
		return temp;
		
	 }; 
    
    virtual T binLower(long bin) {
        return (*this)[0][bin];
	};
	virtual T binUpper(long bin) {        
        return (*this)[1][bin];
	};
	virtual T binMid(long bin) {
		return 0.5*( binLower(bin) + binUpper(bin) );
	};
	virtual T binSize(long bin) {
		return binUpper(bin) - binLower(bin);
	};
    
    virtual size_t numEdges() { return this->size()+1; };
    
	virtual size_type size() {
		size_type j;
		this->DynMatrix<T>::size(NULL,&j);
		return j;
	};
	int findbin (T searchval, T min, T max, int imax) {
		int bin;
		switch ( lookup_method ) {
            case 1:
                bin = floor( (searchval - bounds[0]) / (bounds[1] - bounds[0] ) );
                break;
			default: 
				int guess = floor( (searchval - bounds[0]) / (bounds[1] - bounds[0] ) );
				bin = binsearch(searchval, min, max, guess, imax);
				break;
		};
		
		#if EDGEDEBUG 
		cout << "Bin search: "<<searchval<<", "<< bin << "\n";
		#endif
		return bin;
		
	};
	long findbin (T searchval) {
		long imax = this->size()-1;
        if (imax == -1) return -1;
		
		return findbin(searchval, bounds[0], bounds[1], imax);
	};
    virtual long lookup (T searchval) {
		return this->findbin(searchval);
	};
	long binsearch (T& searchval, T& minval, T& maxval, int i, int& imax) {
		
		if ( searchval >= maxval ) return imax;
		else if ( searchval < minval ) return -1;
        else if ( i < 0 ) return -1;
		
		if ( i > imax ) i = binsearch(searchval, minval, maxval, i-1, imax);
		if ( searchval >= (*this)[0][i] && searchval < (*this)[1][i] ) return i;
		
		if ( searchval < (*this)[0][i] ) i = binsearch(searchval, minval, maxval, i-1, imax);
		else if ( searchval > (*this)[1][i] ) i = binsearch(searchval, minval, maxval, i+1, imax);
		
		return i;
	};
	
	EdgeSet<T>& operator=(EdgeSet<T>& that) {
		if (this->size() != that.size()) { this->reinit(that.size(), (T) 0); }
		
		that.getBounds(this->bounds);
		this->lookup_method = that.lookup_method;
		
		this->DynMatrix<T>::operator=(that);

		return *this;
	};
	
	bool operator==(EdgeSet<T>& that) {
		if (this->size() != that.size()) return 0;
	
		if (this->bounds[0] != that.bounds[0]) return 0;
		if (this->bounds[1] != that.bounds[1]) return 0;
		
		if (this->lookup_method != that.lookup_method) return 0;
		
        //Possibly add code to check each edge value
		return 1;
	};
	
	friend ostream& operator<<(ostream &output, EdgeSet& edges) {
		int bins = edges.size();
		output << "EdgeSet Bins: "<<bins<<"\n";		
		for (int k=0;k<bins;k++) {
			//sprintf(line, "
			output << "[" ;
			output << setw(7) << edges[0][k];
			output << " : ";
			output << setw(7) << edges[1][k];
			output << "]\n";
		}
		output << "\n";
		
		return output;
	};

};



/***********************
 
 Class: EdgeFunction
 
 A class which encapsulates a function which converts a
 value of arbitrary precision (type T, specified by client)
 to an integer position.  The function parameters, such as min,max,
 and number of bins, are stored in the object. The base class
 EdgeFunction simply implements linear or log binning (SetLookupMethod =0,1 for linear, 2 for Log).
 
 It is designed to be sub-classed for more complicated combinatoric expressions; only the parameters
 and not the functional form can vary at run-time.
 
 The second class, EdgeSet, encapsulates an arbitrary list of edges and does a simple recursive searh
 on them.  For large sets of non-uniform binning, such as lightcurves, it is probably best to write
 a container of edge functions under which the time range is divided into intervals of homogenous binning.
 The container would select the proper function at lower resolution, and the funciton provide high resolution, like a B-tree.
 However, tree structures could also be devised to use EdgeSet type lookups instead of functional lookups.
 
 ************************/
template<typename T>
class EdgeFunction : public EdgeSet<T> {
    
private:
	int lookup_method;
    T bounds[2];
    T binFactor;
    T linStep;
	size_t edges;
	int useLog;
    
    
public:
	EdgeFunction() {
		useLog=0;
		this->lookup_method = 0;
		this->SetParameters( (T) 0, (T) 0, 1);
	};
    ~EdgeFunction() {};
    
	virtual void UseLog(int flog) { useLog = flog; lookup_method=1; };
	virtual size_t numEdges() { return this->edges; };
    virtual size_t size() { return this->edges-1; };
	virtual void SetLookupMethod(int i) {
		this->lookup_method = i;
	};
	virtual void getBounds(T * minmax) {
		minmax[0] = this->bounds[0];
		minmax[1] = this->bounds[1];
	};
	virtual void SetParameters(T min, T max, size_t channels) {
		if (min > max) {
			T temp = max;
			max = min;
			min = temp;
		}
		this->bounds[0] = min;
		this->bounds[1] = max;
		this->edges = channels+1;
		
		//if (min == 0.0) min = 1.0f / (max + 1.0f) ;
        double f = pow( (max / min), 1.0f / double( channels ) );
		this->binFactor = (T) f;
		this->linStep = (max - min) / channels;
        
	};
	
    T multiplyRule(long bin) {
        if (useLog) return pow( this->binFactor, bin );
        else return this->linStep * bin;
    };
    
	virtual T binLower(long bin) {
        if (useLog)
            return this->bounds[0] * pow( this->binFactor, bin );
        else
            return this->bounds[0] + this->linStep * bin;
	};
	virtual T binUpper(long bin) {        
        if (useLog)
            return this->bounds[0] * pow( this->binFactor, bin+1 );
        else
            return this->bounds[0] + this->linStep * (bin+1);
	};
	virtual T binMid(long bin) {
		return 0.5*( binLower(bin) + binUpper(bin) );
	};
	virtual T binSize(long bin) {
		return binUpper(bin) - binLower(bin);
	};
	
	virtual T operator[](long i) { return binLower(i); };
	virtual T operator()(long i) { return binMid(i); };
	
	virtual long int lookup(T searchval) {
		long int bin = 0;
		T low;
		T high;
		switch ( lookup_method ) {
            case 1:
				if (useLog) 
					bin = floor ( (log10( searchval ) - log10( bounds[0] )) / log10( binFactor ) );
				else 
					bin = floor( (searchval - bounds[0]) / this->linStep );
                break;
				
			case 2:
				low = bounds[0];
				high = low*binFactor;
				while ( bin < this->edges && (low <= searchval && high > searchval) ) {
					low = high;
					high *= binFactor;
					bin++;
				}
				break;
				
			default: 
				if (useLog) 
					bin = floor ( (log10( searchval ) - log10( bounds[0] )) / log10( binFactor ) );
				else 
					bin = floor( (searchval - bounds[0]) / this->linStep );
				break;
		};
		return bin;
	};
	
	
	//friend class EdgeSet;
};



#endif