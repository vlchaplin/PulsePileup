/*
 *  spoccExeUtilities.h
 *  
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef SP_EXE_UTITILS_H
#define SP_EXE_UTITILS_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <list>
#include <algorithm>

#include <xlocale.h>

#include "DBStringUtilities.hh"
#include "VectorSpan.h"

#define PARSER_SEPERATOR_MATCH (const char*)","
#define PARSER_RANGE_TOKEN (const char*)"-"
#define PARSER_MAX_SEPARATRIX_ELEMENTS 100

#ifndef UTC_FMT_DEFAULT
#define UTC_FMT_DEFAULT "%Y-%m-%d"
#endif

using namespace std;



char** getCmdOption(char ** begin, char ** end, const std::string & option,int nfields=0, int max=0);

void printSingleTimeHelp();
void printTimeArgHelp();

bool getTimeArgMET(char ** argsBegin, char ** argsEnd, double& met, string& timesys, char ** nextArg=NULL, int maxArgs=0);
bool getTimeRangeArgs( char ** argsBegin, char ** argsEnd, double& metStart, double& metEnd, string& timesys );

struct ArgumentInfo {
    string optionName;
    string optionDesc;
    int col1w;
    int col2w;
    string sep;
    ArgumentInfo() { col1w=15; col2w=60; sep="  "; };
    
    ArgumentInfo(string name, string desc, int w1=15, int w2=60, string fsep="  ") {
        optionName=name;
        optionDesc=desc;
        col1w=w1;
        col2w=w2;
        sep=fsep;
    };
};

inline ArgumentInfo NewArgument( string name, string desc, int w1=15, int w2=60, string fsep="  ")
{
    ArgumentInfo temp(name,desc,w1,w2,fsep);
    return temp;
}

ostream& operator<<( ostream& outs, ArgumentInfo line );

inline vector<string> split( string& input, const char * delim = PARSER_SEPERATOR_MATCH, size_t guessAvgFieldLength = 3) 
{
	vector<string> substrings;
	
	substrings.reserve( input.size() / (guessAvgFieldLength+1) );
	
	size_t p0,p1;
	size_t spx = 0;
	size_t nchars;

	p0 = 0;
	p1 = input.find_first_of( delim );

	while ( p1 != string::npos ) {
	
		if (spx == 0) {
			nchars = p1;
			substrings.push_back( input.substr( p0, nchars ) );
		} else {
			nchars = p1 - p0 - 1;
			substrings.push_back( input.substr( p0+1, nchars ) );
		}
		
		p0 = p1;
		p1 = input.find_first_of( delim, p0+1 );
		spx++;
		
	}
	
	if ( p0+1 < input.size() ) substrings.push_back( input.substr( p0+1 ) );
	
	if ( spx == 0 ) {
		substrings.clear();
	}

	return substrings;
}

inline void parseRecordSelction(string& input, VectorSpan< size_t, size_t >& items) {
	
	//cout <<  "Input = " << input << endl;
	
	vector<string> ells = split( input, PARSER_SEPERATOR_MATCH );
	vector<string> range;
	vector<string>::iterator rit1;
	vector<string>::iterator rit2;
	size_t subrangeB;
	size_t subrangeE;

	size_t t = 0;
	size_t sz = ells.size();
	size_t thisItem;
	
	if ( sz == 0 ) {
		ells.push_back( input );
		sz++;
	}
	
	while (t < sz) {
		//cout << "Element: " << ells[t] << endl;
		range = split( ells[t], PARSER_RANGE_TOKEN );
		
		if (range.size() == 1) {
			t++;
		}
		else if (range.size() >= 2) {
			subrangeB = atoi( range[0].c_str() );
			subrangeE = atoi( range[1].c_str() );
			//cout << subrangeB << " <-> " << subrangeE << endl;
			if (subrangeE > subrangeB && subrangeB > 0 )
				items.Insert(subrangeB, subrangeE+1);
	
		} else {
			thisItem = (size_t)atoi( ells[t].c_str() );
			if ( thisItem > 0 ) items.Insert(thisItem, thisItem+1);
		}
		t++;
	}
}


inline void parseRecordSelction(string& input, list<size_t>& items) {
	
	//cout << input << endl;
	
	vector<string> ells = split( input, PARSER_SEPERATOR_MATCH );
	vector<string> range;
	vector<string>::iterator rit1;
	vector<string>::iterator rit2;
	size_t subrangeB;
	size_t subrangeE;

	size_t t = 0;
	size_t sz = ells.size();
	size_t thisItem;
	
	if ( sz == 0 ) {
		ells.push_back( input );
		sz++;
	}
	
	while (t < sz) {
		//cout << "Element: " << ells[t] << endl;
		range = split( ells[t], PARSER_RANGE_TOKEN );
		if (range.size() != 0) {
			subrangeB = atoi( range[0].c_str() );
			subrangeE = atoi( range[1].c_str() );
			//cout << subrangeB << " <-> " << subrangeE << endl;
			if (subrangeB > 0 && subrangeE > 0)
				while (subrangeB <= subrangeE) 
					items.push_back( subrangeB++ );
	
		} else {
			//No range token found, so it must be a single number
			thisItem = (size_t)atoi( ells[t].c_str() );
			if ( thisItem > 0 ) items.push_back( thisItem );
		}
		t++;
	}
}



#endif