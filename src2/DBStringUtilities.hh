/*
 *  DBStringUtilities.hh
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef DBSTRINGUTILS_HH
#define DBSTRINGUTILS_HH

#include <ctime>
#include <cmath>
#include <sys/stat.h>
#include <stddef.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>

#define GLAST_EPOCH_SHIFT time_t(978328800)
#define MJD_AT_EPOCH 51910.0
#define NLEAPSECONDS 2
// *** Leap second boundaries are defined in DBStringUtilities.cpp

using namespace std;

#ifdef _MSC_VER
#include <direct.h>
inline char * path_sep() {
	return (char*)"\\";
};
inline int makeDir(const char * dir) {
	return _mkdir(dir);
};
#else
inline char * path_sep() {
	return (char*)"/";
};
#endif


inline bool FileExists( const char * file ) {
	struct stat stFileInfo; 
	return ( stat(file,&stFileInfo) == 0 );
};
inline bool FileExists( string file ) {
	struct stat stFileInfo; 
	return ( stat(file.c_str(),&stFileInfo) == 0 );
};

inline int makeDir(const char * dir) {
	char cmd[300];
	strcpy(cmd, (char*)"mkdir -p ");
	strcat(cmd, dir);
	return system(cmd);
};
inline int makeDir(string dir) {
	char cmd[300];
	strcpy(cmd, "mkdir -p ");
	strcat(cmd, dir.c_str() );
	return system(cmd);
};

inline int makeSymLink(const char * sym, const char * actual) {
	char cmd[300];
	strcpy(cmd, "ln -s ");
	strcat(cmd, actual);
	strcat(cmd, " ");
	strcat(cmd, sym);
	cout << cmd << endl;
	return system(cmd);
};

inline int makeSymLink(string sym, string actual) {
	char cmd[300];
	strcpy(cmd, "ln -s ");
	strcat(cmd, actual.c_str() );
	strcat(cmd, " ");
	strcat(cmd, sym.c_str() );
	cout << cmd << endl;
	return system(cmd);
};

inline void test_FailsIfNosprintf() {
	char test[20];
	char fmt[] = "%d";
	sprintf( test, fmt, 1234 );
};

inline void test_FailsIfTooManyArgs() {
	char test[20];
	char fmt[] = "%d";
	sprintf( test, fmt, 1, 2, 3, 4, 5, 6 );
};

inline char * get_basename( char * file, char * dirname=NULL ) {
	
	char * splitStarts[100];
	char ** posItr;
	char ** endItr;
	
	splitStarts[0] = file;
	posItr = splitStarts;
	endItr = splitStarts+100;

	char * pathsep = path_sep();

	while ( *posItr != NULL && posItr != endItr ) 
	{
		*(posItr+1) = strpbrk( (*posItr), pathsep);
		//cout << (*posItr) << endl;
		posItr++;
		
		//Drop the leading matched character from the trailing substring
		if ( *posItr != NULL ) (*posItr)++;
	}
	
	if (dirname!=NULL) {
		size_t nchars = (size_t)( *(posItr-1) - file );
		strncpy( dirname, file, nchars );
		
		if (sizeof(dirname) >= nchars) dirname[nchars] = '\0';
	}
	
	return *(posItr-1);
};


/**** class Glast_Epoch *****
	Create a singleton which has the time conversion from GLAST MET's epoch
	on January 1, 2001 to the system epoch ( UNIX is January 1 1970 ).
****/

class Glast_Epoch {
	
	private:
	//static Glast_Epoch * the_instance;
	static bool exists;
	
	static double leap_times[NLEAPSECONDS];
	
	time_t relative_epoch;
	Glast_Epoch();
	
	public:
	void SetEpoch( time_t epoch ) {
		this->relative_epoch = epoch;
		//cout << "set epoch: "<< this->relative_epoch << endl;
	};
	//Adapted from met_leapsecs() in GbmTime.pm by Bill Cleveland.
	inline long double met_leapsecs(long double& met) {
		
		long double leaps = 0;
		
		int i=0;
		while (i<NLEAPSECONDS) {
			if ( leap_times[i++] < met ) leaps++;
			else return leaps;
		}
		return leaps;
	};
	inline double met_leapsecs(double& met) {
		
		double leaps = 0;
		
		int i=0;
		while (i<NLEAPSECONDS) {
			if ( leap_times[i++] < met ) leaps++;
			else return leaps;
		}
		return leaps;
	}
	
	inline time_t seconds() {
		return this->relative_epoch;
	};
	inline double mjd() {
		return MJD_AT_EPOCH;
	};
	static Glast_Epoch* get_instance();
	
};

//void sod_to_hhmmss( double sod, char * time);

struct tm * glast_met2local( long double met, long double * fraction = NULL );


//Adapted from met2mjd() in GbmTime.pm by Bill Cleveland.
inline double glast_met2mjd( long double met ) {
	Glast_Epoch * epoch = Glast_Epoch::get_instance();
	
	met -= epoch->met_leapsecs(met);
	
	return epoch->mjd() + (met / 86400.0);
};
inline long double glast_mjd2met( double mjd ) {
	Glast_Epoch * epoch = Glast_Epoch::get_instance();
	
	long double met = ( mjd - epoch->mjd() )*86400.0;
	met += epoch->met_leapsecs(met);
	
	return met;
};

inline void glast_met_to_utc( long double met, char * datetime){

	struct tm * gmet = glast_met2local( met );
	
	strftime( datetime, 80, (char*)"%Y-%m-%d %H:%M:%S", gmet );
	
	//cout << datetime << endl;
};
inline string glast_met_to_utc( long double met ){

	struct tm * gmet = glast_met2local( met );
	char datetime[20];
	strftime( datetime, 80, (char*)"%Y-%m-%d %H:%M:%S", gmet );
	
	return string(datetime);
};

inline void glast_met_strftime( long double met, char * datetime, char * formatStr){

	struct tm * gmet = glast_met2local( met );
	
	strftime( datetime, 80, formatStr, gmet );
	
};

inline string glast_met_strftime( long double met, char * formatStr){
    
	struct tm * gmet = glast_met2local( met );
	static char datetime[1000];
	strftime( datetime, 80, formatStr, gmet );
	
	return string(datetime);
	//cout << datetime << endl;
};

inline double glast_met_to_startofday( long double met ) {
	
	struct tm * gmet = glast_met2local( met );
	
	gmet->tm_sec = 0;
	gmet->tm_min = 0;
	gmet->tm_hour = 0;
	
	Glast_Epoch * epoch = Glast_Epoch::get_instance();
	
	double st = double( timegm( gmet ) - epoch->seconds() );
	
	st += epoch->met_leapsecs(st);
	
	return st;
};

double glast_tm_to_met( struct tm * gmet, double secfraction=0.0  );
double glast_utc_to_met( char * datetime, char fmt[]=(char*)"%Y-%m-%d %H:%M:%S" );

inline double glast_met_to_nearest_utc_unit( long double met, int tm_sec=-1, int tm_min=-1, int tm_hour=-1, int tm_mday=-1, int tm_mon=-1  ) {
	
	struct tm * gmet = glast_met2local( met );
	
	if ( tm_sec >= 0) gmet->tm_sec = tm_sec;
	if ( tm_min >= 0) gmet->tm_min = tm_min;
	if ( tm_hour >= 0) gmet->tm_hour = tm_hour;
	if ( tm_mday >= 0) gmet->tm_mday = tm_mday;
	if ( tm_mon >= 0) gmet->tm_mon = tm_mon;
	
	return glast_tm_to_met( gmet );
};

inline int glast_met_integer_dayfraction( long double met, int precision=3 ) {
    
    double sod = met - glast_met_to_startofday(met);
    
    //cout << glast_met_to_startofday(met) << endl;
    //cout << sod << endl;
    
    return round( (sod/86400.0)*pow( 10.0, precision ) );
}

inline void vectorizeFormat( vector<string>& chans, unsigned int nchan, const char * format, unsigned int start=0) {
	
	char chstr[100];
	vector<string>::iterator it;

	chans.resize(nchan);
	it = chans.begin();
	for(unsigned int ch=start; ch < nchan; ch++) {
		//sprintf( chstr, format, ch );
		sprintf( chstr, format, ch, ch, ch, ch, ch, ch ); //Allow as many as six %d in the format str
		*it = string(chstr);
		it++;
	}
};

inline void vectorizeColumns( vector<string>& chans, unsigned int nchan, string pre, string post, unsigned int start=0) {
	
	stringstream chstr;
	vector<string>::iterator it;

	chstr.str("");
	chans.resize(nchan);
	it = chans.begin();
	for(unsigned int ch=start; ch < nchan; ch++) {
		chstr << ch;
		string entry = pre+chstr.str()+post;
		*it = entry;
		it++;
		chstr.str("");
	}

};

template<typename X>
inline string str_join(const vector< X >& vec,const string sep)
{
		string tmp;
		stringstream strstr;
        if(vec.size()==0)
                return "";
/*				
        string::size_type size=sep.length()*vec.size();
		
        for(unsigned int i=0;i<vec.size();i++)
        {
                size += vec[i].size();
        }

        tmp.reserve(size);
	*/	
		strstr << vec[0];
        tmp = strstr.str();
        for(unsigned int i=1;i<vec.size();i++)
        {
			strstr.str("");
			strstr << vec[i];
			tmp = tmp + sep + strstr.str();
        }
        return tmp;
}

template<typename X, typename Y>
inline string str_thread_join(const vector< X >& xvec, const string xsepy, const vector< Y >& yvec, const string xy_sep)
{
	return str_thread_join( xvec, xsepy, yvec, xy_sep, 0,0);
}
template<typename X, typename Y>
inline string str_thread_join(const vector< X >& xvec, const string xsepy, const vector< Y >& yvec, const string xy_sep, size_t xoff, size_t yoff)
{
		string tmp;
		stringstream strstr;
        size_t i,size;
		
		if(  xvec.size() > yvec.size() ) size = yvec.size();
		else size = xvec.size();
		
        if (size == 0 || xoff >= size || yoff >= size ) return "";
				
		strstr << xvec[xoff] << xsepy << yvec[yoff];
        tmp = strstr.str();
		
        for(i=1; i<size; i++)
        {
			strstr.str("");
			strstr << xvec[i+xoff] << xsepy << yvec[i+yoff];
			tmp = tmp + xy_sep + strstr.str();
        }
        return tmp;
}


namespace sql_stmts {

	template<typename T>
	inline string where_TimeIntervalEnclosure( string timeColName, T& lowerBound, T& upperBound )
	{
		stringstream tempTimeCondition;
		
		tempTimeCondition << "(" << timeColName << " >= " << lowerBound << " AND " << timeColName << " <= " << upperBound << ")";
		return tempTimeCondition.str();
	}






}



//#include "DBStringUtilities.cpp"

#endif
