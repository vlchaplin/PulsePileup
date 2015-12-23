/*
 *  DBStringUtilities.cpp
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef DBSTR_CPP
#define DBSTR_CPP

//#define epoch_shift 978328800

#include "DBStringUtilities.hh"

//Glast_Epoch * Glast_Epoch::the_instance = NULL;

bool Glast_Epoch::exists = 0;
double Glast_Epoch::leap_times[] = { 157766399.0, 252460800.0 };

Glast_Epoch::Glast_Epoch() {

	struct tm epoch_tm = { 0 };
	epoch_tm.tm_sec = 0;
	epoch_tm.tm_min = 0;
	epoch_tm.tm_hour = 0;
	epoch_tm.tm_mday = 1;
	epoch_tm.tm_year = 101; //years since 1900
	epoch_tm.tm_yday = 0;
	epoch_tm.tm_isdst = 0;

/*	time_t stattime;
	time(&stattime);
	struct tm * currTime = localtime( &stattime );
*/

	this->relative_epoch = timegm( &epoch_tm );
	
	//Since we are using the standard mktime() function,
	//we are local time units, and thus must account for daylight
	//savings.
	//this->relative_epoch = GLAST_EPOCH_SHIFT;
	//if (currTime->tm_isdst) this->relative_epoch -= 3600;
}

Glast_Epoch* Glast_Epoch::get_instance() {
	
	static Glast_Epoch * the_instance;
	
	//cout << "get_instance()" << endl;
	if (! Glast_Epoch::exists) {
		Glast_Epoch::exists = 1;
		the_instance = new Glast_Epoch;	
	}
	//cout << "self epoch: " << the_instance->seconds() << endl;
	return the_instance;
};

struct tm * glast_met2local( long double met, long double * fraction ) {
	
	Glast_Epoch * epoch = Glast_Epoch::get_instance();
	met -= epoch->met_leapsecs(met);
	
	time_t met_int = floor(met);

	if ( fraction != NULL ) *fraction = met - met_int;
		
	time_t met_epoch = epoch->seconds() + met_int;
	
	struct tm * gmet = gmtime( &met_epoch );
	//cout << "Is DST " << gmet->tm_isdst << endl;
	//gmet->tm_isdst=-1;
	return gmet;
};

double glast_tm_to_met( struct tm * gmet, double secfraction  ) {
	
	Glast_Epoch * epoch = Glast_Epoch::get_instance();
	//cout << "Is DST " << gmet->tm_isdst << endl;

	//gmet->tm_isdst=-1;

	double st = double( timegm( gmet ) - epoch->seconds() ) + secfraction;
	
	st += epoch->met_leapsecs(st);
	
	return st;
};

double glast_utc_to_met( char * datetime, char * fmt ) {

	struct tm t0 = { 0 };
	
	double secfraction=0.0;
	char * decimalPntptr = strpbrk( fmt, "." );
	char * lastchar;
	
	if ( decimalPntptr != NULL ) {
		
		char * decimalPosition = decimalPntptr;
	
		int n_precision=strspn( decimalPntptr+1, "n" );
		
		decimalPntptr = strpbrk( datetime, "." );
		
        //cout << "here: " << " nprec: " << n_precision << endl;
        
		if ( n_precision > 0 && decimalPntptr!=NULL ) {

			static int digit;			
			static int k=1;

			while ( decimalPntptr[0] != '\0' && k <= n_precision ) {
								
				digit = decimalPntptr[0] - '0';
				
				if (digit <= 9 && digit >= 0) secfraction += (digit / pow( 10.0, k++ ) );
				
				decimalPntptr++;
			}
			
			//cout << "here2: " << " fsec: " << secfraction << endl;
		}

		size_t dpos = decimalPosition - fmt;
		
		//decimalPosition[0] = '\0';
        
		
		static char tempFmt[1000];
		strncpy( tempFmt, fmt, dpos );
		tempFmt[dpos] = '\0';
		lastchar = strptime( datetime, tempFmt, &t0 );
        
        //cout << "here3" << endl;
        
	} else {
		lastchar = strptime( datetime, fmt, &t0 );
	}

	if (lastchar == NULL ) {
		cout << "Failed to parse utc time" << endl;
		return -1;
	}
    
    //cout << "final: " << " fsec: " << secfraction << endl;
    
	return glast_tm_to_met( &t0, secfraction );

}

/*
//Convert seconds of day to hh:mm:ss string
void sod_to_hhmmss( double sod, char * time) {
	
	int hours,minutes;

	hours = floor( sod / 3600.0 );
	sod -= hours * 3600.0;
	
	minutes = floor( sod / 60.0 );
	sod -= minutes * 60.0;

	sprintf( time, "%2d:%2d:%2lf", hours, minutes, sod );
}
*/
#endif
