/*
 *  PHA_IO.cpp
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef PHAIODEFINED
#define PHAIODEFINED

//#define PHA1TEMPLATE "/gbm/Analysis/LightcurveSimulator/fits_templates/pha1template.pha"
//#define TESTOUTFILE "/gbm/Analysis/LightcurveSimulator/data/goc_test_pha1.pha"

#define COPYMOD_KEY_OR_CARD 1

#include "PHA_IO.hh"

using namespace std;

template<typename T>
int PHAReader::ReadPHAKeys( OGIP_PHA_misc<T> * keys )
{
	if (keys == NULL) return -1;
	
	fitsfile *fptr;
	int status=0;
	
	ffopen(&fptr, (char*)this->file.c_str(), READONLY, &status);
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		return status;
	}
	
	this->bindKeysToMem( *keys );
//	ffuky(fptr, TFLOAT, (char *)"AREASCAL", &(set.ogip_fields->areasc), NULL, &status);



	vector< FitsKeyLocation >::iterator k1 = fileToMemLink.begin();
	while (k1 != fileToMemLink.end() ) 
	{
		ffmahd( fptr, k1->hdunum, NULL, &status);
		ffgky( fptr, k1->ftype,  (char*)(k1->key.c_str()), k1->ptr, NULL, &status);
		k1++;
	}

	
	
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		return status;
	}
	
	ffclos(fptr, &status);
	return status;
}

template<typename C,typename E,typename T>
int PHAReader::ReadDataFile(SinglePHA<C,E,T> * newset) {

	if (newset == NULL) return -1;

	int namelength = this->file.size();
	int status = 0;
	
	char * filename;
	filename = new char[namelength+1];
	
	strcpy(filename, this->file.c_str() );
	
	int ratecol, errcol;
	long nchan;
	
	
	
	fitsfile *fptr;
	
	ffopen(&fptr, filename, READONLY, &status); 
	
	delete [] filename;
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		return status;
	}
	char detname[72];
	ffgky(fptr, TSTRING, (char *)"DETNAM", detname, NULL, &status);
	if ( status==0 ) {
		string temp( detname );
		newset->setDetName( temp );
	}
	
	ffmnhd(fptr, ANY_HDU, (char *)"EBOUNDS", 0, &status);
	ffgnrw(fptr,&nchan,&status);
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		return status;
	}
	int b;
	float * chdata;
	chdata = new float[nchan];

	EdgeSet<E> * ebounds;
	ebounds = new EdgeSet<E>;
	ebounds->reinit(nchan, (E) 0.0);
	
	//Read Emin, Emax.  Use one float array to minimize heap allocations.
	ffgcv(fptr,TFLOAT,2,1,1,nchan,0,chdata,NULL,&status);
	for (b=0;b<nchan;b++) { (*ebounds)[0][b] = chdata[b]; }
	
	ffgcv(fptr,TFLOAT,3,1,1,nchan,0,chdata,NULL,&status);
	for (b=0;b<nchan;b++) { (*ebounds)[1][b] = chdata[b]; }
	
    ebounds->calcBounds();
	
	//cout << *ebounds;
	
	ffmnhd(fptr, ANY_HDU, (char *)"SPECTRUM", 0, &status);
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		delete [] chdata;
		delete ebounds;
		return status;
	}
	
	ffgcno(fptr, CASEINSEN, (char *)"STAT_ERR", &errcol, &status);
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		delete [] chdata;
		delete ebounds;
		return status;
	}
	
	bool isRates = 1;
	
	ffgcno(fptr, CASEINSEN, (char *)"RATE", &ratecol, &status);
	if (status != 0) ffgcno(fptr, CASEINSEN, (char *)"COUNTS", &ratecol, &status);
	if (status != 0) {  
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		delete [] chdata;
		delete ebounds;
		return status;
	} else { isRates = 0; }
	
	EdgeSet<E> * temp = newset->getEdgeSet();
	bool temp_is_shared = newset->sharededges();
	newset->UseBinning(ebounds);
	
	ffgcv(fptr,TFLOAT,ratecol,1,1,nchan,0,chdata,NULL,&status);
	for (b=0;b<nchan;b++) { (*newset)[b].cval = chdata[b]; }
	
	ffgcv(fptr,TFLOAT,errcol,1,1,nchan,0,chdata,NULL,&status);
	for (b=0;b<nchan;b++) { (*newset)[b].stat_err = chdata[b]; }
	
	float exp = -1;
	ffgky(fptr, TFLOAT, (char *)"EXPOSURE", &exp, NULL, &status);
	if ( status==0 ) {
		for (b=0;b<nchan;b++) { (*newset)[b].cval *= exp; }
		newset->SetExposure( (T)exp );
	}
	
	if (status != 0) { 
		newset->UseBinning(temp);
		ffrprt( stderr, status);
		ffclos(fptr, &status);
		delete [] chdata;
		delete ebounds;
		return status;
	}
	
	ffclos(fptr, &status);
	
	newset->setParentFile( this->file );
	
	if (! temp_is_shared ) delete temp;
	
	delete [] chdata;
	return status;
};





template<typename E>
int PHAWriter::WriteEboundsHDU(fitsfile ** fptr, int * status, EdgeSet<E> *& temp) {
		
		if (temp == NULL) return -2;
		
		int nchan = temp->size();
		if (nchan == 0) return -1;
		
		int * chanindex;
		chanindex = new int[nchan];
		
		E *emin,*emax;
		emin = new E[nchan];
		emax = new E[nchan];
		
		for (int i=0;i<nchan;i++) {
			emin[i] = (*temp)[0][i];
			emax[i] = (*temp)[1][i];
			chanindex[i] = i;
		}
		
		int eb_fitstype = type_to_fitstype< E >();
		
		ffuky(*fptr, TINT, (char *)"DETCHANS", &nchan, NULL, status);
		
		ffpcl(*fptr, TINT, 1,1,1, nchan, chanindex, status);
		ffpcl(*fptr, eb_fitstype, 2,1,1, nchan, emin, status);
		ffpcl(*fptr, eb_fitstype, 3,1,1, nchan, emax, status);
#ifndef NO_DATE_STAMP    
        ffpdat(*fptr, status);
#endif
		ffpcks(*fptr, status);
		
		delete [] chanindex;
		delete [] emin;
		delete [] emax;
		
		return *status;
	
	return 0;
}

template<typename C,typename E,typename T>
int PHAWriter::formatSpecHDUCLAS(fitsfile** fptr, int * status, SinglePHA<C,E,T> * data) {
	
	bool ChAreRates;
	bool poissonErr;
	typename SinglePHA<C,E,T>::size_type nchan;
	
	ChAreRates = data->IsDifferential();
	poissonErr = data->UsePoissonErrors();
	nchan = data->nchannels();
	
	ffuky(*fptr, TSTRING, (char *)"HDUCLAS", (char*)"OGIP", NULL, status);
	ffuky(*fptr, TSTRING, (char *)"HDUVERS", (char*)"1.2.1", NULL, status);
	ffuky(*fptr, TSTRING, (char *)"HDUVERS1", (char*)"1.1.1", NULL, status);
	ffuky(*fptr, TSTRING, (char *)"HDUCLAS1", (char*)"SPECTRUM", NULL, status);
	
	if (this->specCompositionLevel == 0) 
		ffuky(*fptr, TSTRING, (char *)"HDUCLAS2", (char*)"BKG", (char*)"Data is background only", status);
	else if (this->specCompositionLevel == 1) 
		ffuky(*fptr, TSTRING, (char *)"HDUCLAS2", (char*)"NET", (char*)"Data is source - background", status);
	else 
		ffuky(*fptr, TSTRING, (char *)"HDUCLAS2", (char*)"TOTAL", (char*)"Data is source + background", status);

	char tcol[20];
	if ( ChAreRates ) {
		//If the spectrum is in rates, write a stat_err column
		
		ffuky(*fptr, TSTRING, (char *)"HDUCLAS3", (char*)"RATE", (char *)"Data is counts/s", status);
		
		ttype_to_tform(TFLOAT, tcol);
		ffdcol( *fptr, pha_col, status );
		fficol( *fptr, pha_col, (char*)"RATE", tcol, status);
		
		sprintf(tcol, "TUNIT%d", pha_col); 
		ffuky(*fptr, TSTRING, tcol, (char*)"counts/s", NULL, status);
		
		//Add the stat_err column
		ttype_to_tform(TFLOAT, tcol);
		int temp;
		ffgcno( *fptr, CASEINSEN, (char *)"STAT_ERR", &temp, status);
		if (*status==0) ffdcol( *fptr, temp, status ); else *status=0;
		
		sprintf(tcol, "TUNIT%d", err_col); 
		fficol( *fptr, err_col, (char*)"STAT_ERR", (char*)"E", status);
		ffuky( *fptr, TSTRING, tcol, (char*)"counts/s", NULL, status);
	} else {
		ffuky(*fptr, TSTRING, (char *)"HDUCLAS3", (char*)"COUNT", (char *)"Data is counts", status);
		//If not, write counts.  Update the column type to be a long integer.
		ttype_to_tform(TLONG, tcol);
		ffdcol( *fptr, pha_col, status );
		fficol( *fptr, pha_col, (char*)"COUNTS", tcol, status);
		
		sprintf(tcol, "TUNIT%d", pha_col);
		ffuky(*fptr, TSTRING, tcol, (char*)"count", NULL, status);
		
		int temp;
		ffgcno( *fptr, CASEINSEN, (char *)"STAT_ERR", &temp, status);
		if (*status==0) ffdcol( *fptr, temp, status ); else *status=0;
		fficol( *fptr, err_col, (char*)"STAT_ERR", (char*)"E", status);
		sprintf(tcol, "TUNIT%d", err_col); 
		ffuky( *fptr, TSTRING, tcol, (char*)"count", NULL, status);
		
	}
	if ( poissonErr ) {
		int tf = 1;
		ffuky(*fptr, TLOGICAL, (char *)"POISSERR", &tf, NULL, status);
		//int temp;
		//ffgcno( *fptr, CASEINSEN, (char *)"STAT_ERR", &temp, status);
		//if (*status==0) ffdcol( *fptr, temp, status ); else *status=0;
		//temp=0;
		//ffuky(*fptr, TINT, (char *)"STAT_ERR", &temp, NULL, status);
	} else {
		int tf = 0;
		ffuky(*fptr, TLOGICAL, (char *)"POISSERR", &tf, NULL, status);
	}

	if (*status != 0) {
		cout << "Error formating SPECTRUM HDU" << endl;
		ffrprt( stderr, *status);
	}

	return *status;
};

//*fptr Must already point to the Spectrum HDU.
template<typename C,typename E,typename T>
int PHAWriter::WriteSpectrumHDU(fitsfile ** fptr, int * status, SinglePHA<C,E,T>& set) {

	typedef typename SinglePHA<C,E,T>::err_type err_type;

	int nchan = set.size();
	if (nchan == 0) return -1;
	
	int tlmin = 0;
	int tlmax = nchan-1;
	
	int chtype,errtypef,t_typef;
	
	this->formatSpecHDUCLAS( fptr, status, &(set) );
	
	int * chanindex;
	chanindex = new int[nchan];
	
	C * channeldata;
	channeldata = new C[nchan];
	
	err_type * stat_err;
	stat_err = new err_type[nchan];
	
	for (int i=0;i<nchan;i++) {
		channeldata[i] = set[i].cval;
		stat_err[i] = set[i].stat_err;
		chanindex[i] = i;
	}
	
	chtype = type_to_fitstype< C >();
	t_typef = type_to_fitstype< T >();
	errtypef = type_to_fitstype< err_type >();
	
	T exp = set.GetExposure();
	
	ffuky(*fptr, TINT, (char *)"DETCHANS", &nchan, NULL, status);
	ffuky(*fptr, TINT, (char *)"TLMIN1", &tlmin, NULL, status);
	ffuky(*fptr, TINT, (char *)"TLMAX1", &tlmax, NULL, status);
	ffuky(*fptr, t_typef, (char *)"EXPOSURE", &exp, NULL, status);
		
	ffpcl(*fptr, TINT, 1/*colnum*/, 1/*startrow*/, 1/*firstell*/, nchan/*nelems*/, chanindex/*data*/, status);
	ffpcl(*fptr, chtype, pha_col,1,1, nchan, channeldata, status);
	ffpcl(*fptr, errtypef, err_col,1,1, nchan, stat_err, status);
#ifndef NO_DATE_STAMP    
	ffpdat(*fptr, status);
#endif
	ffpcks(*fptr, status);
			
	delete [] chanindex;
	delete [] channeldata;
	delete [] stat_err;
	
	return *status;
};

template<typename Oclass>
int PHAWriter::WriteStandardKeys(fitsfile ** fptr, int * status, Oclass *& fields) {

	if (*status != 0) return *status;

	string copyout;

	if ( fields->getDetName(copyout) ) {
		char detname[80];
		strcpy(detname, copyout.c_str() );
		ffuky(*fptr, TSTRING, (char *)"DETNAM", detname, NULL, status);
	} else {
		ffukyu(*fptr, (char *)"DETNAM", NULL, status);
	}
	
	if (fields->object.length() > 0) {
		ffuky(*fptr, TSTRING, (char *)"OBJECT", (char*)fields->object.c_str(), NULL, status);
	}
	
	char date_obs[30];
	char date_end[30];
	
	glast_met_strftime(  fields->tmin, date_obs, (char*)"%Y-%m-%dT%H:%M:%S" );
	glast_met_strftime(  fields->tmax, date_end, (char*)"%Y-%m-%dT%H:%M:%S" );
	
	int T_fits = type_to_fitstype(fields->tzero);
	ffuky(*fptr, T_fits, (char *)"TRIGTIME", &(fields->tzero), NULL, status);
	ffuky(*fptr, T_fits, (char *)"TSTART", &(fields->tmin), NULL, status);
	ffuky(*fptr, T_fits, (char *)"TSTOP", &(fields->tmax), NULL, status);
	ffuky(*fptr, TSTRING, (char *)"DATE-OBS", date_obs, NULL, status);
	ffuky(*fptr, TSTRING, (char *)"DATE-END", date_end, NULL, status);

	return *status;
}

template<typename C,typename E,typename T>
int PHAWriter::WriteDataFile(SinglePHA<C,E,T>& set) {

	string name;
	
    int status = 0;
	int nchan = set.size();
	if (nchan == 0) return -1;
	
	set.getParentFile(name);
	int namelength = name.size();
	
	
	char * templatefl;
	if (namelength == 0) {
        
        name = this->getTemplateFile();
        namelength = name.size();
        
        if ( namelength == 0 ) {  
            templatefl = (char *)PHA1TEMPLATE;
            if ( templatefl == NULL ) {
                cout << "Error, PHA1TEMPLATE not defined" << endl;
                return -1;
            }
        } 
    }
	
    if ( namelength > 0) {
		templatefl = new char[namelength+1];
		strcpy(templatefl, name.c_str() );
	}
    
    //cout << "file to make = " << this->file << endl;
    //cout << "template = " << templatefl << endl;
		
	fitsfile *fptr;
	fitsfile *tplt;

	ffinit(&fptr, this->file.c_str(), &status);
	
	if (status == 104 || status == 105) {
		remove ( this->file.c_str() );
		status = 0;
		ffinit(&fptr, this->file.c_str(), &status);
	}
	
	if (status != 0) goto WritePHADataErrorClosure;

	ffopen(&tplt, templatefl , READONLY, &status);
	if (status != 0) goto WritePHADataErrorClosure;
	
	
	//Copy the Primary HDU to fptr
	ffcpfl( tplt, fptr, 1,1,0, &status);
	ffuky( fptr, TSTRING, (char *)"CREATOR", (char *)"PHA_IO.cpp", NULL, &status);
	ffuky( fptr, TSTRING, (char *)"FILENAME", (char *)this->file.c_str(), NULL, &status);
	
	
	typedef typename SinglePHA<C,E,T>::fields_class fields_class;
	fields_class * fields;
	fields = set.getMiscFields();
	
	this->WriteStandardKeys(&fptr, &status, fields);
	
#ifndef NO_DATE_STAMP    
	ffpdat(fptr, &status);
#endif
	ffpcks(fptr, &status);
	
	
	long n_rows;
	char* fitsTableName;
	fitsTableName = (char*)"EBOUNDS";
	ffmnhd(tplt, ANY_HDU, fitsTableName, 0, &status);
	ffcphd(tplt, fptr, &status);
	
	ffmnhd( fptr, ANY_HDU, fitsTableName, 0, &status);
	ffgnrw( fptr, &n_rows, &status);
	ffdrow( fptr, 1, n_rows, &status);
	if (status != 0) goto WritePHADataErrorClosure;

	EdgeSet< E > * edgetable;
	edgetable = set.getEdgeSet();
	this->WriteEboundsHDU(&fptr, &status, edgetable);
	this->WriteStandardKeys(&fptr, &status, fields);
	
	if (status != 0) goto WritePHADataErrorClosure;
#ifndef NO_DATE_STAMP    
	ffpdat(fptr, &status);
#endif
	ffpcks(fptr, &status);
	
	
	fitsTableName = (char*)"SPECTRUM";
	//Copy the extension from the template file
	ffmnhd( tplt, ANY_HDU, fitsTableName, 0, &status);
	ffcphd( tplt, fptr, &status);
	
	//Now delete the rows written from the template file
	ffmnhd( fptr, ANY_HDU, fitsTableName, 0, &status);
	ffgnrw( fptr, &n_rows, &status);
	ffdrow( fptr, 1, n_rows, &status);
	if (status != 0) goto WritePHADataErrorClosure;
	
	this->WriteSpectrumHDU(&fptr, &status, set);
	this->WriteStandardKeys(&fptr, &status, fields);
	
	ffuky(fptr, TFLOAT, (char *)"AREASCAL", &(set.ogip_fields->areasc), NULL, &status);
	ffuky(fptr, TFLOAT, (char *)"BACKSCAL", &(set.ogip_fields->backsc), NULL, &status);
	ffuky(fptr, TUSHORT,(char *)"GROUPING", &(set.ogip_fields->grpg), NULL, &status);
	ffuky(fptr, TUSHORT,(char *)"QUALITY", &(set.ogip_fields->qual), NULL, &status);
	
	char rsp[100];
	set.CpyRmf_file(rsp);
	
	//cout << "RSP FILE= " << rsp << endl;
	
	if ( rsp == NULL ) ffuky(fptr, TSTRING, (char *)"RESPFILE", (char*)"None", NULL, &status);
	else ffuky(fptr, TSTRING, (char *)"RESPFILE", rsp, NULL, &status);
	
	if ( fields->bakfile.length() > 0 ) {
		ffuky(fptr, TSTRING, (char *)"BACKFILE", (char*)fields->bakfile.c_str(), NULL, &status);
	}

#ifndef NO_DATE_STAMP    
	ffpdat(fptr, &status);
#endif
	ffpcks(fptr, &status);
	if (status != 0) goto WritePHADataErrorClosure;
	
	cout << endl<< "----------------------"<<endl;
	cout << "Wrote single spectrum to:" << this->file << endl;
	cout << "----------------------"<<endl<<endl;
	
	WritePHADataNormalClosure:
	ffclos(fptr, &status);
	ffclos(tplt, &status);
	return status;
	
	WritePHADataErrorClosure:
	if (namelength > 0) delete [] templatefl;
	ffrprt( stderr, status);
	goto WritePHADataNormalClosure;
};

#endif
