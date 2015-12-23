/*
 *  PHAStructures.cpp
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef PHAICPP 
#define PHAICPP

#include "PHAStructures.hh"
#include <iostream>

using namespace std;

template<typename Cnts, typename Etype, typename T>
void SinglePHA<Cnts,Etype,T>::UseBinning(EdgeSet<Etype> * edges, bool delete_old) {

	if (delete_old && (! sharing[1]) && Ebounds != NULL) {
		if ( Ebounds != edges ) delete Ebounds;
		Ebounds = edges;
	} else if ( Ebounds != edges ) Ebounds = edges;
	
	if (edges == NULL) {
		this->reInit(0);
		return;
	}
	
	if ((size_t)nchan != Ebounds->size()) {
		reInit( Ebounds->size() );
	}
};

template<typename Cnts, typename Etype, typename T>
void SinglePHA<Cnts,Etype,T>::UseBinning(EdgeSet<Etype> * edges) { 
	this->UseBinning(edges, 0); 
};

template<typename Cnts, typename Etype, typename T>
int SinglePHA<Cnts,Etype,T>::lookup(Etype domain) {
	if (Ebounds == NULL) return -1;
	return Ebounds->findbin(domain);
};

template<typename Cnts, typename Etype, typename T>
void SinglePHA<Cnts,Etype,T>::LinBin(Etype min, Etype max, int channels) {
	if (Ebounds == NULL) this->Ebounds = new EdgeSet<Etype>;
	Ebounds->MakeLinearBins(min, max, channels);
	this->UseBinning(Ebounds,0);
};

template<typename Cnts, typename Etype, typename T>
void SinglePHA<Cnts,Etype,T>::LogBin(Etype min, Etype max, int channels) {
	if (Ebounds == NULL) this->Ebounds = new EdgeSet<Etype>;
	Ebounds->MakeLogBins(min, max, channels);
	this->UseBinning(Ebounds,0);
};

template<typename Cnts, typename Etype, typename T>
bool SinglePHA<Cnts,Etype,T>::AddCount(int chan) {
	if ((size_type)chan < 0 || (size_type)chan >= nchan) {
		return 0;
	}
	data[chan].cval++;
	return 1;
};

template<typename Cnts, typename Etype, typename T>
bool SinglePHA<Cnts,Etype,T>::AddCounts(int chan, Cnts counts) {
	if (chan < 0 || chan >= nchan) {
		return 0;
	}
	data[chan].cval+=counts;
	return 1;
};

template<typename Cnts, typename Etype, typename T>
bool SinglePHA <Cnts,Etype,T>::AddEvent(Etype energy) {
	int bin = this->lookup(energy);
	if ( bin == -1 ) return 0;
        
	return this->AddCount(bin);
};



#endif

/*
bool GBMResponse::loadMatrix() {

	if (file.length() <= 4) return 0;

	//ffopen

	return 0;
};
*/

