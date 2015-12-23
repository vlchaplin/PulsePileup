/*
 *  PHA_IO.hh
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef PHAIOUNIT
#define PHAIOUNIT

#include "EdgeSet.hh"
#include "PHAStructures.hh"
//#include "TTEventTable.hh"
#include "DBStringUtilities.hh"

#include <string>
#include <vector>
#include <iostream>
#include <typeinfo>

#ifndef PHA1TEMPLATE
    #define PHA1TEMPLATE getenv("PHA1TEMPLATE")
#endif

using namespace std;

struct FitsKeyLocation {
	string key;
	int ftype;
	int hdunum;
	void * ptr;
};


inline void ttype_to_tform(int ttype, char * colform) {
    
	string tform;
    
	switch (ttype) {
		case TBIT		: tform = "X"; break;
		case TBYTE		: tform = "B"; break;
		case TSTRING	: tform = "A"; break;
		case TSHORT		: tform = "I"; break;
		case TUSHORT	: tform = "U"; break;
		case TLONG		: tform = "I"; break;
		case TLONGLONG	: tform = "K"; break;
		case TULONG		: tform = "U"; break;
		case TFLOAT		: tform = "E"; break;
		case TDOUBLE	: tform = "D"; break;
		case TCOMPLEX	: tform = "C"; break;
		case TINT		: tform = "I"; break;
		default : tform = "\0";
	};
	
	strcpy(colform, tform.c_str());
	
};

inline size_t sizeof_ttype(int ttype) {
    
	size_t nbytes;
    
	switch (ttype) {
		case TBIT		: nbytes = 0; break;
		case TBYTE		: nbytes = 1; break;
		case TSTRING	: nbytes = sizeof(char); break;
		case TSHORT		: nbytes = sizeof(int); break;
		case TUSHORT	: nbytes = sizeof(unsigned int); break;
		case TLONG		: nbytes = sizeof(long int); break;
		case TLONGLONG	: nbytes = sizeof(long long int); break;
		case TULONG		: nbytes = sizeof(unsigned long int); break;
		case TFLOAT		: nbytes = sizeof(float); break;
		case TDOUBLE	: nbytes = sizeof(double); break;
		case TCOMPLEX	: nbytes = 0; break;
		case TINT		: nbytes = sizeof(int); break;
		default : nbytes = 0;
	};
	
	return nbytes;
};

template<typename memtype>
int type_to_fitstype( ) {
    
	int coltype;
	
	if ( typeid(memtype) == typeid(int) ) coltype = TINT;
	else if ( typeid(memtype) == typeid(short) ) coltype = TSHORT;
	else if ( typeid(memtype) == typeid(unsigned short) ) coltype = TUSHORT;
	else if ( typeid(memtype) == typeid(float) ) coltype = TFLOAT;
	else if ( typeid(memtype) == typeid(double) ) coltype = TDOUBLE;
	else if ( typeid(memtype) == typeid(long) ) coltype = TLONG;
	else if ( typeid(memtype) == typeid(char) ) coltype = TSTRING;
	else if ( typeid(memtype) == typeid(unsigned char) ) coltype = TBYTE;
	else if ( typeid(memtype) == typeid(long double) ) coltype = TDOUBLE;
	else if ( typeid(memtype) == typeid(long long) ) coltype = TLONGLONG;
	else coltype = TSTRING;
    
	return coltype;
};

template<typename memtype>
int type_to_fitstype(memtype& object) {
    
	int coltype;
	
	if ( typeid(memtype) == typeid(int) ) coltype = TINT;
	else if ( typeid(memtype) == typeid(short) ) coltype = TSHORT;
	else if ( typeid(memtype) == typeid(unsigned short) ) coltype = TUSHORT;
	else if ( typeid(memtype) == typeid(float) ) coltype = TFLOAT;
	else if ( typeid(memtype) == typeid(double) ) coltype = TDOUBLE;
	else if ( typeid(memtype) == typeid(long) ) coltype = TLONG;
	else if ( typeid(memtype) == typeid(char) ) coltype = TSTRING;
	else if ( typeid(memtype) == typeid(unsigned char) ) coltype = TBYTE;
	else if ( typeid(memtype) == typeid(long double) ) coltype = TDOUBLE;
	else if ( typeid(memtype) == typeid(long long) ) coltype = TLONGLONG;
	else coltype = TSTRING;
    
	return coltype;
};

/*
 
 check_for_archive - simple utility routine to check whether a string represents a possibly compressed file name
                        by checking suffixes '.tar', '.gz', '.tgz', and '.tar.gz'
 */

inline void check_for_archive(char * inputFile, bool * isTar, bool * isZip, char ** archive, char ** relativeFile) {


	//cout << "Checking for archive" << endl;

	char * temp;
	char * temp2;
	size_t arcbase_nchars;
	
	(*archive)[0]='\0';
	(*relativeFile)[0]='\0';
	*isTar=0;
	*isZip=0;

	temp = strstr(inputFile, (char*)".tar");
	
	if ( temp != NULL ) {
		
		//cout << ".tar";
		
		temp += 4;
		*isTar=1;
		
		temp2 = strstr(temp, (char*)".gz");
		
		if ( temp2 != NULL ) {
			
			//cout << ".gz";
			
			temp2 += 3;
			*isZip=1;
			
			arcbase_nchars = (size_t)(temp2 - inputFile);
			
			//copy from temp+1 since the path seperator is there...
			if ( temp2[0] != '\0' ) strcpy( *relativeFile, temp2+1 );
			
		} 
		else
		{
			arcbase_nchars = (size_t)(temp - inputFile);
			
			if ( temp[0] != '\0' ) strcpy( *relativeFile, temp+1 );
		}
		
		strncpy( *archive, inputFile, arcbase_nchars );
		
		(*archive)[arcbase_nchars]='\0';
		
		return;
	}
	
	
	temp = strstr(inputFile, (char*)".tgz");
	
	if ( temp != NULL ) {
		
		//cout << ".tgz";
		
		*isTar=1;
		*isZip=1;
		temp += 4;
		
		if ( temp[0] != '\0' ) strcpy( *relativeFile, temp+1 );
		
		arcbase_nchars = (size_t)(temp - inputFile);
		
		strncpy( *archive, inputFile, arcbase_nchars );
		
		(*archive)[arcbase_nchars]='\0';
		
		return;
	}
	
	temp = strstr(inputFile, (char*)".gz");
	
	if ( temp != NULL ) {
		
		*isTar=0;
		*isZip=1;
		temp += 3;
		
		if ( temp[0] != '\0' ) strcpy( *relativeFile, temp );
		
		arcbase_nchars = (size_t)(temp - inputFile);
		
		strncpy( *archive, inputFile, arcbase_nchars );
		
		(*archive)[arcbase_nchars]='\0';
		
		return;
	}
	
	
};

/***
 
 CLASS: FxbReadStat
 
 Interior class designed to be used inside the PHA IO objects, though could be used out right by clients
 ***/

class FxbReadStat {
	public:
	fitsfile * fptr;
	int status;
	long n_eff_rows;
	long n_rows;
	long n_cols;
	
	bool isTar;
	bool isZip;
	bool unlinkOnClose;
	
	bool CloseOnError;
	
	FxbReadStat() {
		fptr=NULL;
		status=0;
		n_eff_rows=0;
		n_rows=0;
		n_cols=0;
		isTar=0;
		isZip=0;
		unlinkOnClose=0;
		CloseOnError=1;


		//cout << "Fxb CTOR " << this << endl;

	};
	
	~FxbReadStat() {
		//cout << "Fxb DTOR fptr:" << fptr << endl;
        


	};
	
	inline bool okForOps() {
		return ( this->fptr != NULL ) && ( this->status == 0 );
	}
	
	inline int template_init(char * file, char * tempFile, bool overwrite=1) {
		
		string ovrw;
		
		if (file==NULL || tempFile==NULL || strlen(file)*strlen(tempFile)==0 ) return -1;
		
		if ( overwrite ) ovrw = "!"; else ovrw = "";
		
		string extendedFileName = ovrw + string(file) + "(" + string(tempFile) + ")";
		
		ffinit( &fptr, extendedFileName.c_str(), &status );
		
		unlinkOnClose=0;
		if (status) this->errmsg();
	
		return status;
	};
	
	inline int open(char * file, int extno=0, int mode=READONLY ) {
	
		this->status=0;
	
		ffopen( &(this->fptr), file, mode, &(this->status) );
		
		if (this->status != 0) {
			
			char basename[200];
			char relative[100];
			char command[310];
		
			char * b_ptr = basename;
			char * r_ptr = relative;
			char ** b_ptrptr = &b_ptr;
			char ** r_ptrptr = &r_ptr;
		
			check_for_archive( file, &isTar, &isZip, b_ptrptr, r_ptrptr );
			
			if ( isTar || isZip ) {
			
				if ( isTar && strlen(relative) == 0 ) {
					cout << "Unable to extract from a tar archive without a relative filename appended" << endl;
					return -1;
				}
			
				if ( FileExists((const char *)relative) ) unlinkOnClose=0;
				else if ( isTar && isZip ) {
					unlinkOnClose=1;
					sprintf( command, (char*)"tar -xzvf %s %s", basename, relative );
					cout << command << endl;
					
					system( command );
				} else if ( isTar ) {
					unlinkOnClose=1;
					sprintf( command, (char*)"tar -xvf %s %s", basename, relative );
					cout << command << endl;
					
					system( command );
				}
				
				this->status = 0;
				ffopen( &(this->fptr), relative, mode, &(this->status) );
				
			}
			
		
		}
		
		if (this->status != 0) {  
			if (this->status != 0) {
				if (CloseOnError)
					return this->errclose();
				else 
					return this->errmsg();
			}
		} else if (extno != 0) {
			ffmahd(this->fptr, extno, NULL, &(this->status) );
			if ( this->chdu_is_table() )
			{
				ffgnrw(this->fptr,&(this->n_rows),&(this->status) );
				ffgrsz(this->fptr,&(this->n_eff_rows),&(this->status) );
			
				if ( n_eff_rows > n_rows ) n_eff_rows = n_rows;
			}			
			if (this->status != 0) {
				if (CloseOnError)
					return this->errclose();
				else 
					return this->errmsg();
			}
		}
		return 0;
	};
	
	inline int move2hdu(char * extname) {
	
		if (this->fptr == NULL)  return -1;
	
		if (this->status != 0) {  
			if (this->status != 0) {
				if (CloseOnError)
					return this->errclose();
				else 
					return this->status;
			}
		} else {
			ffmnhd(this->fptr, ANY_HDU, extname, 0, &(this->status) );
			if ( this->chdu_is_table() )
			{
				ffgnrw(this->fptr,&(this->n_rows),&(this->status) );
				ffgrsz(this->fptr,&(this->n_eff_rows),&(this->status) );
			
				if ( n_eff_rows > n_rows ) n_eff_rows = n_rows;
			}
			
			if (this->status != 0) {
				if (CloseOnError)
					return this->errclose();
				else 
					return this->status;
			}
		}
		return 0;
	};
	inline int move2hdu(int extno) {
		if (this->status != 0 || this->fptr == NULL) {  
			if (this->status != 0) {
				if (CloseOnError)
					return this->errclose();
				else 
					return this->status;
			}
		} else {
			ffmahd(this->fptr, extno, NULL, &(this->status) );
			
			if ( this->chdu_is_table() )
			{
				ffgnrw(this->fptr,&(this->n_rows),&(this->status) );
				ffgrsz(this->fptr,&(this->n_eff_rows),&(this->status) );
			
				if ( n_eff_rows > n_rows ) n_eff_rows = n_rows;
			}
			
			if (this->status != 0) {
				if (CloseOnError)
					return this->errclose();
				else 
					return this->status;
			}
		}
		return 0;
	};
	
	inline bool chdu_is_table() {
		int chdu;
		ffghdt(this->fptr, &chdu, &(this->status) );
		
		if ( chdu == BINARY_TBL || chdu == ASCII_TBL ) return 1;
		else return 0;
	}
	
	inline int close() {
		unlinkOnClose=0;
		
		if ( unlinkOnClose ) {
			int st=0;
			char file[300];
			file[0]='\0';
			ffflnm(fptr, file, &st );
			
			if ( strlen(file) > 0 && FileExists(file) ) remove( file );
		}
		
		
		ffclos(fptr, &status );
		
		
		fptr = NULL;
		isTar=0;
		isZip=0;
		unlinkOnClose=0;
		
		return status;
	
	};
	inline int errmsg() {
		ffrprt( stderr, status );
		return status;
	}
	inline int errclose() {
		if (this->fptr == NULL)  return -1;
		ffrprt( stderr, status);		
		return this->close();
	};
    inline int up_checksum_date() {
        if (this->status != 0 || this->fptr == NULL) return this->status;
#ifndef NO_DATE_STAMP    
        ffpdat(fptr, &status);
#endif
        ffpcks(fptr, &status);
        return this->status;
    };
    inline int flush() {
        if (this->status != 0 || this->fptr == NULL) return this->status;

        ffflus(this->fptr, &status);
        
        return this->status;
    };
	
	inline int read_key( char * key, int t_type, void* data, int hdu_num=-1, char * comm=NULL )
	{	
		if (this->fptr == NULL)  return 104;
		if (this->status) return status;
	
		int curr_hdu;
		ffghdn(this->fptr, &curr_hdu);
		
		if ( hdu_num != curr_hdu && hdu_num > 0 ) {
			this->move2hdu( hdu_num );
			ffgky( this->fptr, t_type, key, data, comm, &status );
			this->move2hdu( curr_hdu );
		} else 
			ffgky( this->fptr, t_type, key, data, comm, &status );
		
		return status;
	};
    inline int write_key( char * key, int t_type, void* data, int hdu_num=-1, char * comm=NULL )
	{	
		if (this->fptr == NULL)  return 104;
		if (this->status) return status;
        
		int curr_hdu;
		ffghdn(this->fptr, &curr_hdu);
		
		if ( hdu_num != curr_hdu && hdu_num > 0 ) {
			this->move2hdu( hdu_num );
			ffuky( this->fptr, t_type, key, data, comm, &status );
			this->move2hdu( curr_hdu );
		} else 
			ffuky( this->fptr, t_type, key, data, comm, &status );
		
		return status;
	};
    
    inline int delete_key( char * key, int hdu_num=-1 )
	{	
		if (this->fptr == NULL)  return 104;
		if (this->status) return status;
        
		int curr_hdu;
		ffghdn(this->fptr, &curr_hdu);
		
		if ( hdu_num != curr_hdu && hdu_num > 0 ) {
			this->move2hdu( hdu_num );
			ffdkey(this->fptr, key, &status);
			this->move2hdu( curr_hdu );
		} else 
			ffdkey(this->fptr, key, &status);
		
		return status;
	};
	
    
    inline int modify_comment( char * key, char * comment, int hdu_num=-1 )
	{	
		if (this->fptr == NULL)  return 104;
		if (this->status) return status;
        
		int curr_hdu;
		ffghdn(this->fptr, &curr_hdu);
		
		if ( hdu_num != curr_hdu && hdu_num > 0 ) {
			this->move2hdu( hdu_num );
			ffmcom(this->fptr, key, comment, &status);
			this->move2hdu( curr_hdu );
		} else 
			ffmcom(this->fptr, key, comment, &status);
		
		return status;
	};
    
    
    inline long col_vector_size( int colnum )
    {
        if (this->fptr == NULL || this->status)  return -1;
        
        int type;
        long repeat; long numbytes;
        
        ffgtcl(this->fptr, colnum, &type, &repeat, &numbytes, &(this->status) );
        
        if ( this->status ) return -1;
        
        return repeat;
        
    };
    
    template <typename T>
    inline int read_col( int colnum, long rowA, long rowB, T * data, long colrepeat=1 )
    {
        if (this->fptr == NULL)  return -1;
		if (this->status) return this->status;
        
        if ( this->n_rows <= 0 ) return -1;
        
        long rowItr = rowA;
        
        if ( rowA <= 0 ) rowItr = 1;
        if ( rowB <= 0 ) rowB = this->n_rows;
        
        long readSize = this->n_eff_rows;
        
        int ttype = type_to_fitstype<T>();
        
        long x=0;
        
        while ( rowItr <= rowB && this->status == 0 )
        {
            if ( rowItr + readSize > rowB ) readSize = rowB - rowItr + 1;
            
            ffgcv(this->fptr, ttype, colnum, rowItr, 1, readSize*colrepeat, NULL, data + x*readSize*colrepeat, NULL, &(this->status) );
            
            x++;
            rowItr += readSize;
        }
        
        return this->status;
        
    }
    
};


class PHAIO_Base {

	public:
	
	string file;
	FxbReadStat stat;
	
    string templateFile;
    
	long lastOpRowSize;

	FxbReadStat * getFxb() { return &(this->stat); };
	
	PHAIO_Base() {};
	~PHAIO_Base() {};
		
	long operationRows() { return this->lastOpRowSize; };
	long incOpRows(long inc) { return (this->lastOpRowSize+=inc); };
	long setOpRows(long val) { return (this->lastOpRowSize=val); };
	
	int open( string filename, int extno=0, int mode=READWRITE ) {
	
		FxbReadStat * fStat = this->getFxb();
	
		this->setFile( filename );
	
		return fStat->open( (char*)file.c_str(), extno, mode );
	};
	
	int open( int extno=0, int mode=READWRITE ) {
	
		FxbReadStat * fStat = this->getFxb();
	
		return fStat->open( (char*)file.c_str(), extno, mode );
	};
	
	int reopen( PHAIO_Base * obj ) {
		FxbReadStat * fStat = this->getFxb();
		FxbReadStat * oStat = obj->getFxb();
		return ffreopen( (oStat->fptr), &(fStat->fptr), &(fStat->status) );
	};
	
	int reopen( fitsfile * ifptr ) {
		FxbReadStat * fStat = this->getFxb();
		return ffreopen( ifptr, &(fStat->fptr), &(fStat->status) );
	};
	
	int close() {

		FxbReadStat * fStat = this->getFxb();

		fStat->close();
		
		return 0;
	};
	
	int status() {
		FxbReadStat * fStat = this->getFxb();
		return fStat->status;
	};
	int errmsg() {
		FxbReadStat * fStat = this->getFxb();
		return fStat->errmsg();
	};
	int errclose() {
		FxbReadStat * fStat = this->getFxb();
		return fStat->errclose();
	};
	
	int set_status(int newstatus) {
	
		FxbReadStat * fStat = this->getFxb();
	
		fStat->status=newstatus; 
		return newstatus;
	};
	
	inline int move2hdu(int extno) {
		FxbReadStat * fStat = this->getFxb();
		return fStat->move2hdu(extno);
	};
	inline int move2hdu(char * extname) {
		FxbReadStat * fStat = this->getFxb();
		return fStat->move2hdu(extname);
	};
	int up_checksum_date() {
        FxbReadStat * fStat = this->getFxb();
        return fStat->up_checksum_date();
    };
	int read_key( char * key, int t_type, void* data, int hdu_num=-1, char * comm=NULL ) {
		FxbReadStat * fStat = this->getFxb();
		return fStat->read_key(key, t_type, data, hdu_num, comm);
	};
   
    int write_key( char * key, int t_type, void * data, int hdu_num=-1, char * comm=NULL ) {
		FxbReadStat * fStat = this->getFxb();
		return fStat->write_key(key, t_type, data, hdu_num, comm);
	};
    int delete_key( char * key, int hdu_num=-1) {
		FxbReadStat * fStat = this->getFxb();
		return fStat->delete_key(key,  hdu_num);
	};
    int modify_comment( char * key, char * comment, int hdu_num=-1 ) {
		FxbReadStat * fStat = this->getFxb();
		return fStat->modify_comment(key,  comment, hdu_num);
	};
    
	long n_rows() {
		FxbReadStat * fStat = this->getFxb();
		return fStat->n_rows;
	};
	long n_eff_rows() {
		FxbReadStat * fStat = this->getFxb();
		return fStat->n_eff_rows;
	};

    void setTemplate(char * filename) {
		this->templateFile.assign(filename);
	};
    void setTemplate(string filename) {
		this->templateFile.assign(filename);
	};
    
	void setFile(char * filename) {
		this->file.assign(filename);
	};
	void setFile(string filename) {
		this->file = filename;
	};
	string getFile() {
		return this->file;
	};
    string getTemplateFile() {
		return this->templateFile;
	};
    


};


class PHAReader : public PHAIO_Base {

	
	protected:
	vector< FitsKeyLocation > fileToMemLink;
	
	template<typename T>
	void bindKeysToMem( OGIP_PHA_misc<T>& keys ) {
		
		FitsKeyLocation links[] = { {"TRIGTIME",TDOUBLE,1,(void*)&keys.tzero},
									{"TSTART",TDOUBLE, 1,(void*)&keys.tmin}, 
									{"TSTOP",TDOUBLE,  1,(void*)&keys.tmax},
									{"DETNAM",TSTRING, 1,keys.charDetName},
									{"RESPFILE",TSTRING, 3,keys.charRmfFile} };
		
		this->fileToMemLink.assign( links, links+5 );
	};
	
	public:
	PHAReader(){};
	~PHAReader() {
	
	
	};
	
	template<typename T>
	int ReadPHAKeys( OGIP_PHA_misc<T> * keys );

	template<typename C,typename E,typename T>
	int ReadDataFile(SinglePHA<C,E,T> * newset);

};

class PHAWriter : public PHAIO_Base {

	private:
	
	template<typename C,typename E,typename T>
	int WriteSpectrumHDU(fitsfile** fptr, int * status, SinglePHA<C,E,T>& set);
	
	protected:
	string writephase;
	unsigned int specCompositionLevel;
	int pha_col;
	int err_col;
	
	
	template<typename E>
	int WriteEboundsHDU(fitsfile** fptr, int * status, EdgeSet<E> *& edges);
	template<typename Oclass>
	int WriteStandardKeys(fitsfile** fptr, int * status, Oclass *& fields);
	template<typename C,typename E,typename T>
	int formatSpecHDUCLAS(fitsfile** fptr, int * status, SinglePHA<C,E,T> * data);
	
	public:
	
	PHAWriter(){
		specCompositionLevel=2;
		pha_col=2;
		err_col=3;
	};
	
	bool setSpecComposition(unsigned int setting) {
		//Effects HDUCLAS2: 0= BGK, 1= NET, 2= TOTAL
		if (setting <= 2) {
			this->specCompositionLevel = setting;
			return 1;
		} else return 0;
		
	};
	
	int flush( ) {
		FxbReadStat * fxb = this->getFxb();
		
		return ffflus( fxb->fptr, &(fxb->status) );
	}
	
	template<typename C,typename E,typename T>
	int WriteDataFile(SinglePHA<C,E,T>& set);

};

#include "PHA_IO.cpp"
#endif