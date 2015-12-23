/*
 *  spoccExeUtilities.cpp
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#include "spoccExeUtilities.h"

using namespace std;

ostream& operator<<( ostream& outs, ArgumentInfo line )
{
    if ( line.optionDesc.size() > (size_t)line.col2w ) {
        size_t textLength = line.optionDesc.size();
        size_t first=0;
        size_t numchars=(size_t)line.col2w;
        size_t space = line.sep.size() + line.col1w;
        //bool last=0;
        while ( first < textLength) 
        {
            if (first+numchars+1 > textLength) {
                numchars = textLength - first;
            
                if (first == 0) {
                    outs << right;
                    outs << setw(line.col1w) << line.optionName << line.sep << flush;
                    outs << left << setw(line.col2w) << line.optionDesc.substr(first, numchars);
                }
                else {
                    outs << right << setw(space) << " " << flush;
                    outs << left << setw(line.col2w) << line.optionDesc.substr(first, numchars);
                }
                first+=numchars;
                
            } else{
                
                if (first == 0) {
                    outs << right << setw(line.col1w) << line.optionName << line.sep << flush;
                    outs << left << setw(line.col2w-1) << line.optionDesc.substr(first, numchars-1);
                }
                else {
                    outs << right << setw(space) << " " << flush;
                    outs << left << setw(line.col2w-1) << line.optionDesc.substr(first, numchars-1);
                }
                
                //cout << "|" << line.optionDesc.substr(first+numchars-2,3) << "|" << endl;
                
                if ( line.optionDesc.substr(first+numchars-2,2) == "  " || 
                    (line.optionDesc.substr(first+numchars-1,1) == " " && line.optionDesc.substr(first+numchars,1) != " ") ||
                     (line.optionDesc.substr(first+numchars-2,1) == " " && line.optionDesc.substr(first+numchars-1,1) != " ")
                    ) {
                    outs << " " << endl;
                    first+=(numchars-1);
                }
                else {
                    outs << "-" << endl;
                    first+=(numchars-1);
                }
            }
            
            
            
        }
    } else {
        outs << right << setw(line.col1w) << line.optionName  << flush;
        outs << left << line.sep << setw(line.col2w) << line.optionDesc;
    }
    
    return outs;
};

char** getCmdOption(char ** begin, char ** end, const std::string & option, int nfields, int max)
{

	if ( max > 0 && ( max < end - begin ) ) {
		end = begin + max;
	}

    char ** itr = find(begin, end, option);
    if (itr != end && (itr+nfields) != end)
    {
        return itr;
    }
    return NULL;
}

void printSingleTimeHelp() {
	cout << endl;
	cout << "Time argument format:" << endl;
	cout << "-met nnn.nnn "<<endl<<"   Specify the time in Fermi MET" << endl;
	cout << "-mjd nnn.nnn "<<endl<<"   Specify the time in MJD" << endl;
	cout << "-utc 'timestring' [-fmt 'fmtstr']"<<endl;
	cout << "  Specifies the time in UTC." << endl;
	cout << "  `timestring` is interpreted as a date-time using POSIX strftime() formatting rules and should be quoted."<<endl;
	cout << "  `fmtstr` is the strftime() format string and should also be quoted if present."<<endl;
	cout << "     e.g., a time that looks like 'yyyy-mm-dd hh:mm:ss' is formatted '%Y-%m-%d %H:%M:%S' in POSIX notation." << endl;
	cout << "  See POSIX strftime() documentation for more info." << endl ;
	cout << endl;
	cout << "  As an extension to UTC time formats, this program allows '.nnnn'... at the end of the format code " << endl;
	cout << "  to stand for the fraction of a second with n digits precision. Any digits (0-9) trailing a '.' in " << endl;
	cout << "  the timestring will be interpreted as the fractional part, and will be parsed until the string of" << endl;
	cout << "  digits is broken or until the input precision is reached. Note that the last digit is not rounded." << endl;
	cout << "  '.nnnn' always terminates the format code, regardless of what follows it." << endl;
	cout << "  The timestring may have zero or fewer decimal places than requested in `fmstr`." << endl;
	cout << "     E.g., a format string requesting a time to microsecond precision is '%H:%M:%S.nnnnnn'." << endl;
	cout << "  In this case a properly formatted timestring would be:" << endl;
	cout << "      '08:26:15.1234' (i.e.,seconds interpreted as '15.123400')" << endl;
	cout << endl << "  The default UTC format is '" << UTC_FMT_DEFAULT << "'"<< endl << endl;;
}

void printTimeArgHelp() {
	cout << endl;
	cout << "Time arguments format:" << endl;
	cout << "-all"<<endl<<"   Specifies the entire range of available data" << endl;
	cout << "-met nnn.nnn nnn.nnn"<<endl<<"   Specifies the begin and end time in Fermi MET" << endl;
	cout << "-mjd nnn.nnn nnn.nnn"<<endl<<"   Specifies the begin and end time in MJD" << endl;
	cout << "-utc 'string1' 'string2' [-fmt 'fmtstr']"<<endl;
	cout << "  Specifies the begin and end time in UTC." << endl;
	cout << "  string1 and string2 are interpreted as date-time strings using POSIX strftime() formatting rules."<<endl;
	cout << "     e.g., 'yyyy-mm-dd hh:mm:ss' is '%Y-%m-%d %H:%M:%S' in POSIX notation." << endl;
	cout << "  The default format is '" << UTC_FMT_DEFAULT << "'"<< endl;
	cout << "  See POSIX strftime() documentation for more info." << endl << endl;
}

bool getTimeArgMET(char ** argsBegin, char ** argsEnd, double& met, string& timesys, char ** nextArg, int maxArgs) 
{
	
    
    char ** opt;
	
	bool found=0;
	
	timesys = "utc";
	
	if (! found) {
		opt = getCmdOption( argsBegin, argsEnd, "-met", 1, maxArgs);
		if (opt != NULL) {
			met = atof( *(opt+1) );
			found=1;
			timesys="met";
			//if (nextArg != NULL) nextArg = opt+2;
		} 
	}
	if (! found) {
		opt = getCmdOption( argsBegin, argsEnd, "-mjd", 1, maxArgs);
		if (opt != NULL) {
			
			double mjd1 = atof( *(opt+1) );
			
			met = glast_mjd2met( mjd1 );
			found = 1;
			timesys="mjd";
			
			//if (nextArg != NULL) nextArg = opt+2;

		}
	}
	if (! found) {
        
        //cout << "here" << endl;
        
		opt = getCmdOption( argsBegin, argsEnd, "-utc", 1, maxArgs);
       
		if (opt != NULL) {
			char formatStr[1000];//way overkill
			char ** frm = getCmdOption( argsBegin, argsEnd, "-fmt", 1);
			//if (nextArg != NULL) nextArg = opt+2;

			if (frm != NULL) {
				strcpy(formatStr, *(frm+1) );
				//if (nextArg != NULL) nextArg = opt+3;

			} else {
				strcpy(formatStr, (char*)UTC_FMT_DEFAULT);
			}
        
			met = glast_utc_to_met( *(opt+1), formatStr );
			
			if ( met < 0 ) {
				found=0;
				cout << "error parsing time" << endl;
			} else {
				found = 1;
				timesys="utc";
			}	
		}
        
        
        
	}
    
    
    
//	if ( (! found )|| (getCmdOption( argsBegin, argsEnd, "--help", 0) !=NULL ) ) printSingleTimeHelp();
	
	if (!found) return 0;
	return 1;
    
}

bool getTimeRangeArgs(char ** argsBegin, char ** argsEnd, double& metStart, double& metEnd, string& timesys ) 
{
	char ** opt;
	
	bool found=0;
	
	timesys = "utc";
	opt = getCmdOption( argsBegin, argsEnd, "-all", 0);
	if (opt != NULL) {
		metStart = -1.0;
		unsigned long temp = 0;
		metEnd = (double)~temp;
		
		return 1;
	}
	
	if (! found) {
		opt = getCmdOption( argsBegin, argsEnd, "-met", 2);
		if (opt != NULL) {
			metStart = atof( *(opt+1) );
			metEnd = atof( *(opt+2) );
			found=1;
			timesys="met";
		} 
	}
	if (! found) {
		opt = getCmdOption( argsBegin, argsEnd, "-mjd", 2);
		if (opt != NULL) {
			
			double mjd1 = atof( *(opt+1) );
			double mjd2 = atof( *(opt+2) );
			
			metStart = glast_mjd2met( mjd1 );
			metEnd = glast_mjd2met( mjd2 );
			found = 1;
			timesys="mjd";
		}
	}
	if (! found) {
		opt = getCmdOption( argsBegin, argsEnd, "-utc", 2);
		if (opt != NULL) {
			char formatStr[1000];
			char ** frm = getCmdOption( argsBegin, argsEnd, "-fmt", 1);
			if (frm != NULL) {
				strcpy(formatStr, *(frm+1) );
			} else {
				strcpy(formatStr, (char*)UTC_FMT_DEFAULT);
			}
			struct tm t1 = { 0 };
			struct tm t2 = { 0 };
			char * lastchar1 = strptime( *(opt+1), formatStr, &t1 );
			char * lastchar2 = strptime( *(opt+2), formatStr, &t2 );
			
			//cout << *opt << " : " << *(opt+1) << ", "  << *(opt+2) << endl;
			if (lastchar1 == NULL || lastchar2 == NULL ) {
				cout << "Failed to parse utc time" << endl;
				return 0;
			}
			
			metStart = glast_tm_to_met( &t1 );
			metEnd = glast_tm_to_met( &t2 );
			
			found = 1;
			timesys="utc";
		} 
	}


	if ( (! found )|| (getCmdOption( argsBegin, argsEnd, "--help", 0) !=NULL ) ) printTimeArgHelp();
	
	if (!found) return 0;

	cout << *opt << " => " << setprecision(16) << metStart << ", " << metEnd << endl;
	return 1;

}