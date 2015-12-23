/*
 *  PHAStructures.hh
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef PHASTRUCTURES
#define PHASTRUCTURES

#define EDGEDEBUG 0
#define PHA2_ADDCNT_DEBUG 0
#define PHADEBUG 0

#include "string.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "fitsio.h"

#include "EdgeSet.hh"
#include "PHAElements.hh"

using namespace std;



template<typename T>
T uniform_deviate( int rand )
{
    return (T)rand*(1.0/(RAND_MAX+1.0));
};

template<typename T>
T randomNumber(T low, T hi)
{
    int r = rand();
    T rand = low + uniform_deviate<T>( r ) * (hi - low);
    return rand;
};

template<typename Cnts, typename Etype, typename T>
class SinglePHA {

	public:
	typedef PHAChannelData<Cnts> chan_type;
	typedef OGIP_PHA_misc<T> fields_class;
	typedef typename chan_type::err_type err_type;
	typedef Cnts cunit_type;
	typedef Etype energy_type;
	typedef T time_type;
	typedef typename vector< chan_type >::size_type size_type;
	
	unsigned long int spec_num;

	private:
	vector< chan_type > data;
	typename vector< chan_type >::iterator chvalue; //typename required because iterator type unknown until instantiation (since it's a template). 
	fields_class * ogip_fields;
	
	
	protected:
	/*The following are pointers to allow data sharing ("flyweighting"),but are allocated within the class
	by default.  The method ShareNull() deletes the protected member and sets the pointer to the shared object. 
	
	Clients depend on these pointers being NULL if the object isn't in a valid or defined state,
	thus constructors and accessor methods must set the pointers to NULL where appropriate, since
	some compilers do not initialize pointers to NULL at declaration.
	*/
	
	chan_type * NullEntry; //Null Channel entry
	EdgeSet<Etype> * Ebounds;
	string * parentFile;
	bool sharing[3];
	
	char * rmf_file;
	
	T exposure;
	
	size_t nchan;
	bool isDifferential;
	bool isPoisson;
	
	
	//Protected for encapsulation, except friend class
	virtual EdgeSet<Etype> * getEdgeSet() {
		return this->Ebounds;
	};
	

	public:
	SinglePHA() {
		#if PHADEBUG
		cout << "SinglePHA DEFAULT CONSTRUCTOR\n";
		#endif
		Ebounds = NULL;
		rmf_file = NULL;
		(this->NullEntry) = new chan_type( 0.0, 0.0, 0.0 );
		exposure = (T) 1;
		isDifferential = 0;
		this->reInit(0);
		parentFile = NULL;
		ogip_fields = new fields_class;
		this->isPoisson = 0;
		
		for (int i=0;i<2;i++) sharing[i] = 0;
	};
	SinglePHA(EdgeSet<Etype> * edges, chan_type * nullvall, char * rspfile) {
		#if PHADEBUG
		cout << "SinglePHA SHARED CONSTRUCTOR\n";
		#endif
		(this->ogip_fields) = new fields_class();
		this->rmf_file = rspfile;
		this->isDifferential = 0;
		this->isPoisson = 0;
		this->exposure = (T) 1;
		this->parentFile = NULL;
		
		for (int i=0;i<2;i++) sharing[i] = 1;
		
		this->ShareNull( nullvall );
		this->reInit(0);
		this->ShareBinning( edges );
		
	};
	virtual ~SinglePHA() {
		#if PHADEBUG
		cout << "SinglePHA DESTRUCTOR\n";
		#endif
		if ( (this->Ebounds != NULL) && (! sharing[0] ) ) {
			#if PHADEBUG
			cout << "Ebounds destroyed "<< sharing[0] << "\n";
			#endif
			delete this->Ebounds;
		}
		if ( (this->NullEntry != NULL) && (! sharing[1]) ) delete this->NullEntry;
		
		if (this->parentFile != NULL) delete this->parentFile;
		if (this->ogip_fields != NULL  && (! this->ogip_fields->isShared() ) ) delete this->ogip_fields;
	};
	
	virtual void reInit(int size) {
		chan_type temp;
		temp = *NullEntry;
		this->nchan = size;
		(this->data).clear();
		(this->data).resize( (size_t)size, temp );
	};
	//differential refers to time (rates) at the moment... need to change...
	virtual bool IsDifferential() {
		//return this->isDifferential;
        return this->IsRates();
	};
	virtual void IsDifferential(bool flag) {
		//this->isDifferential = flag;
        this->IsRates(flag);
	};
	virtual bool IsRates() {
		return this->isDifferential;
	};
	virtual void IsRates(bool flag) {
		this->isDifferential = flag;
	};
    
    
	virtual void UsePoissonErrors(bool flag) {
		this->isPoisson = flag;
		
		if ( flag ) {
			chvalue = data.begin();
			while (chvalue != data.end() ) {
				chvalue->stat_err = sqrt( chvalue->cval );
				++chvalue;
			}
		}
	};
	virtual bool UsePoissonErrors() {
		return this->isPoisson;
	};
	
	void getEdgeSet(EdgeSet<Etype>& ptr) {
		ptr = *(this->Ebounds);
	};
    
    int getEminEmax(Etype * emin, Etype * emax)
    {
        if ( this->Ebounds != NULL) {
            Etype temp[2];
            this->Ebounds->getBounds(temp);
            *emin = temp[0];
            *emax = temp[1];
            return 0;
        }
        return 1;
    }
    
	virtual fields_class * getMiscFields() {
		if (this->ogip_fields == NULL) this->ogip_fields = new fields_class();
		return this->ogip_fields;
	};
	chan_type copyNullChan() {
		return *this->NullEntry;
	};
	
	
	virtual void nullify() {
		chvalue = data.begin();
		while (chvalue != data.end() ) {
			*chvalue = *NullEntry;
			++chvalue;
		}
	};
	/*
	void nullify(size_type ch) {        
        if (ch>=0 && ch< this->nchan) {
            this->data[ch] = *NullEntry;
        };
	};*/
	void nullify(unsigned int ch) {        
        if (ch>=0 && ch< this->nchan) {
            this->data[ch] = *NullEntry;
        };
	};
    
    Cnts max() {
		
        Cnts max=0;
        chvalue = data.begin();
		while (chvalue != data.end() ) {
			if (max < chvalue->cval) max = chvalue->cval;
			++chvalue;
		}
		return max;
	};
	
	Cnts total() {
		//This function can be changed to handle the differential state.
		return sum();
	};
	
	Cnts sum() {
		Cnts chsum = 0;
		chvalue = data.begin();
		while (chvalue != data.end() ) {
			chsum += chvalue->cval;
			++chvalue;
		}
		return chsum;
	};
	void setDetName(char * newname) {
		string name = newname;
		fields_class * this_fields = this->getMiscFields();
		this_fields->setDetName(name);
	};
	void setDetName(string& newname) {
		fields_class * this_fields = this->getMiscFields();
		this_fields->setDetName(newname);
	};
	virtual bool getDetName(string& copy) {
		fields_class * this_fields = this->getMiscFields();
		return this_fields->getDetName(copy);
	};
	virtual void setParentFile(string& filename) {
		if (this->parentFile == NULL) this->parentFile = new string;
		this->parentFile->assign(filename);
	};
	virtual bool getParentFile(string& filename) {
		if (this->parentFile == NULL) return 0;
		filename.assign( *(this->parentFile) );
		return 1;
	};
		
	void SetExposure(T seconds) {
		exposure = seconds;
	};
	T GetExposure() {
		return exposure;
	};
	
	
	virtual void UseRmf_file(const char * filename) {
	
		fields_class * fields = getMiscFields();
	
		fields->response.assign(filename);
	};
	virtual void ShareRmf_file(char * filename) {
	
		sharing[2] = 1;
		fields_class * fields = getMiscFields();
	
		fields->response.assign(filename);
	};
	virtual void CpyRmf_file(char * filename) {
	
		fields_class * fields = getMiscFields();
		strcpy( filename, fields->response.c_str() );
		
	};
	
	virtual void ShareFields( fields_class * fields ) {
		if (this->ogip_fields != NULL && this->ogip_fields != fields && (! this->ogip_fields->isShared()) ) delete this->ogip_fields;
		this->ogip_fields = fields;
		this->ogip_fields->setShared();
	};
	
	void UseValueForNull(Cnts nullvall) {
		if (this->NullEntry == NULL) this->NullEntry = new chan_type ( nullvall, 0.0, 0.0 );
		else (this->NullEntry)->cval = nullvall;
	};
	void ShareNull(chan_type* flyw) {
		if (this->NullEntry != NULL && (! sharing[1] ) ) delete this->NullEntry;
		this->NullEntry = flyw;
		(this->sharing)[1] = 1;
	};
	void ShareBinning(EdgeSet<Etype> * edges, bool delete_old) {
		this->UseBinning(edges, delete_old);
		sharing[0] = 1;
	};
	void ShareBinning(EdgeSet<Etype> * edges) {
		this->ShareBinning(edges,0);
	};
	
	void ShareAllFromSpec( SinglePHA<Cnts,Etype,T>* otherSpectrum ) {
		this->ShareFields( otherSpectrum->getMiscFields() );
		this->ShareNull( otherSpectrum->NullEntry );
		this->ShareBinning( otherSpectrum->Ebounds, 1);	
	};
	
	
	bool sharededges() { return sharing[0]; };
	void setsharingflags(bool edges, bool nulls, bool rmf) {
		sharing[0] = edges;
		sharing[1] = nulls;
		sharing[2] = rmf;
	};
	void getsharingflags(bool& edges, bool& nulls, bool& rmf) {
		edges = sharing[0];
		nulls = sharing[1];
		rmf = sharing[2];
	};
	
	virtual size_type size() {
		return this->nchan;}
	;
	virtual size_type nchannels() {
		return this->nchan;}
	;
	
	virtual void UseBinning(EdgeSet<Etype> * edges, bool delete_old);
	virtual void UseBinning(EdgeSet<Etype> * edges);
	int lookup(Etype domain);
	void LinBin (Etype min, Etype max, int channels);
	void LogBin (Etype min, Etype max, int channels);
	bool AddEvent(Etype domainvalue);
	virtual bool AddCount(int channel);
	virtual bool AddCounts(int channel, Cnts counts);
    
    void copySpectrum( Cnts * spectrum, size_t numells, err_type * errors = NULL )
    {
        size_t i;
        
        if ( numells > this->nchannels() ) {
            numells = this->nchannels();
        }
        
        for (i=0; i<numells; i++) {
            data[i].cval = spectrum[i];
        }
        
        if ( errors != NULL) {
            for (i=0; i<numells; i++) {
                data[i].stat_err = errors[i];
            }
        }
        
    }
	
	void phaToVectors( vector<Cnts> * spectrum=NULL, vector<err_type> * errors=NULL, vector<Etype> * chedges=NULL) {
		
		size_type sp_c = this->size();
		
		if ( spectrum != NULL ) {
			spectrum->clear();
			spectrum->reserve(sp_c);
		}
		if ( errors != NULL ) {
			errors->clear();
			errors->reserve(sp_c);
		}
		if ( chedges != NULL ) {
			chedges->clear();
			if (this->Ebounds != NULL) {
				*chedges = (*this->Ebounds)[0];
				chedges->push_back( (*this->Ebounds)[1][sp_c-1] );
			}
		}
		
		chvalue = data.begin();
		while (chvalue != data.end() ) {

			if ( spectrum != NULL ) {
				spectrum->push_back( chvalue->cval );
			}
			if ( errors != NULL ) {
				errors->push_back( chvalue->stat_err );
			}
			
			++chvalue;
		}
	
	}

	//Several of the following operators operate channelwise, and they call the corresponding operator
	//for the class 'chan_type'. Design note: if the this PHA class is ever templated with a variable chan_type, the 
	//variable class must have these operators defined (e.g., chan_type::operator+=(SinglePHA<>& )  ).
	//Basically the PHA is object is a CONTAINER of channel objects, and forwards some operations to the contained channels.

    chan_type& operator [](size_type i) {
		return data[i];
	};
	//const operator in case the object is created as immutable (via const).
	const chan_type& operator [](size_type i) const { 
		return data[i];
	};
	
	SinglePHA<Cnts,Etype,T>& operator=(SinglePHA<Cnts,Etype,T>& that) {
		
		#if PHADEBUG
		cout << "SinglePHA = SinglePHA"<<endl;
		#endif
		
		if (this == &that) return *this;
		
		//Copy channel data (using std::vector assignment), and misc.
		this->data = that.data;
		this->nchan = that.nchan;
		this->IsDifferential( that.IsDifferential() );
		this->UsePoissonErrors( that.UsePoissonErrors() );
		this->exposure = that.exposure;
		
		//Copy null channel value (=> chan_type::operator= must be defined !)
		*(this->NullEntry) = *(that.NullEntry);
		
		//*Safely* Copy edge set into this.
		if ( that.Ebounds != NULL ) {
			EdgeSet<Etype> * edgecopy;
			edgecopy = new EdgeSet<Etype>;
			*(edgecopy) = *(that.Ebounds);
			this->UseBinning( edgecopy, 1);	
		}
			
		//Copy parentFile
		if ( that.parentFile == NULL && this->parentFile != NULL ) delete this->parentFile;
		else if (that.parentFile != NULL) {
			this->setParentFile( *that.parentFile );
		}
		
		if (that.ogip_fields != NULL) {
			*this->ogip_fields = *that.ogip_fields;
		};
		
		return (*this);
	};
    
    SinglePHA<Cnts,Etype,T>& operator+=(SinglePHA<Cnts,Etype,T>& that) {
		//Check identity
		if (this == &that) {
			(*this)*=2;
			return (*this);
		}
	
		if (this->nchan != that.nchan) {
            cerr << "SinglPHA::operator+= sizes do not agree\n";
            return *this;
        }; 
        size_type ch;
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] += that.data[ch];
        };
		
		if (that.ogip_fields != NULL) {
			*this->ogip_fields += *that.ogip_fields;
		};
        return *this;
	};
	SinglePHA<Cnts,Etype,T>& operator+=(Cnts scalar) {
        size_type ch;
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] += scalar;
        };
        return *this;
	};
	SinglePHA<Cnts,Etype,T>& operator-=(Cnts scalar) {
        size_type ch;
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] -= scalar;
        };
        return *this;
	};
	SinglePHA<Cnts,Etype,T>& operator*=(Cnts scalar) {
        size_type ch;
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] *= scalar;
        };
        return *this;
	};
	SinglePHA<Cnts,Etype,T>& operator/=(Cnts scalar) {
        size_type ch;
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] /= scalar;
        };
        return *this;
	};
	SinglePHA<Cnts,Etype,T>& operator-=(SinglePHA<Cnts,Etype,T>& that) {
		//Check identity
		if (this == &that) {
			(*this)*=0;
			return (*this);
		}
	
		if (this->nchan != that.nchan) {
            cerr << "SinglPHA::operator-= sizes do not agree\n";
            return *this;
        }; 
        size_type ch;
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] -= that.data[ch];
        };
		if (that.ogip_fields != NULL) {
			*this->ogip_fields -= *that.ogip_fields;
		};
        return *this;
	};
    
    SinglePHA<Cnts,Etype,T>& operator*=(SinglePHA<Cnts,Etype,T>& that) {
		if (this->nchan != that.nchan) {
            cerr << "SinglPHA::operator*= sizes do not agree\n";
            return *this;
        }; 
        
        size_type ch;
        
        for (ch=0;ch< this->nchan;ch++) {
            this->data[ch] *= that.data[ch];
        };
		if (that.ogip_fields != NULL) {
			*this->ogip_fields *= *that.ogip_fields;
		};
        return *this;
	};

	friend class SpectrumModel;
	friend class PHAReader;
	friend class PHAWriter;

	friend ostream& operator<<(ostream &output, SinglePHA& obj) {
		int bins = obj.size();
		output << "PHA Channels: "<<bins<<"\n";
		//char line[80];
		
		EdgeSet<Etype> * channels = obj.Ebounds;
		
		if (bins == 0 || channels == NULL) {
			output << "\n";
			return output;
		}
		if (obj.parentFile != NULL) {
			output << "Initial file for PHA state (excluding client changes):\n"<< *obj.parentFile << "\n\n";
		};
		
		//Cnts directsum = 0;
		Cnts sum = 0;
		
		float cerror = 0;
		
		Etype integral = 0;
				
		obj.chvalue = obj.data.begin();
		
		
		typename EdgeSet<Etype>::size_type k;
		
		
		if (! obj.sharededges() ) {
			for (k=0;k<bins;k++) {
				//sprintf(line, "
				output << "[" ;
				output << setw(7) << (*channels)[0][k];
				output << " : ";
				output << setw(7) << (*channels)[1][k];
				output << "]  ";
				output << setw(9) << obj.chvalue->cval;
				output << " ( +/- ";
				output << setw(5) << ( obj.chvalue->stat_err );
				output << ", +/- ";
				output << setw(5) << ( obj.chvalue->sys_err );
				output << ")\n";
				cerror += pow(obj.chvalue->stat_err,2);
				integral += ((*channels)[1][k] - (*channels)[0][k])*(obj.chvalue->cval);
				
				sum += obj.chvalue->cval;
				++obj.chvalue;
			}
			output << "Sum { n } = " << sum << " +/- " <<  sqrt(cerror) << endl;
			output << "Integral n*dX = " << integral << endl<<endl;
		} else {
			for (k=0;k<bins;k++) {
				sum += obj.chvalue->cval;
				++obj.chvalue;
			}
			output << "Sum { n } = " << sum << " +/- " <<  sqrt(cerror) << endl<<endl;
		}
		
		return output;
	};
};



/************************** PHA2 Class ********************************

	Container of SinglePHA objects.  Sub-class of SinglePHA.  Most
	base methods are virtualized, but this class forwards the call
	to each contained object.  See comments in class defintion.
	
	Access operators are extended to the time domain, so PHA2 access
	is like a 2-D array, with additional time searches with ().  Notes below.

***********************************************************************/


template<typename Cnts, typename Etype, typename Ttype>
class PHA2 : public SinglePHA<Cnts,Etype,Ttype> {


	private:
	
	//Run-time error catching for use with the () operator.
	class TimeAccessError : public std::runtime_error {
		public:
		TimeAccessError() : std::runtime_error("TimeAccessError") {};
	};



	public:
	typedef SinglePHA<Cnts,Etype,Ttype> pha_type;
	typedef pha_type* pha_ref;
	typedef OGIP_PHA_misc<Ttype> fields_class;
	typedef typename SinglePHA<Cnts,Etype,Ttype>::chan_type chan_type;

	//3/24: Defined pha_ref as pointer to pha_type, and data containers as vectors of pha_ref.
	//Realized that data containers must be pointers to PHAs rather than static 
	//PHA elements.  The static elements don't get proper internal iterators...not
	//sure how to fix this. But also the vector asignment would be expensive, and requires 
	//copying the defined object from the class scope into the vector. Some vector functions
	//take element aliases (pha_type&) instead of copying the data (pha_type), which go out of scope, 
	//and so heap corruption was occuring.
	//Note, the PHA2:: [] and () operators should return pha_type& and not pha_ref,
	//to give clients the [i][c] and (t)[c] syntax for accesing a
	//SinglePHA channel 'c' at a time bin i or time t.
	
	private:
	//Data container and access iterator
	vector< pha_ref > datasets;
	typename vector< pha_ref >::iterator phaitr;
	
	fields_class * ogip_fields;
	
	// The time binning
	EdgeSet<Ttype> * TEdges;
	Ttype tbounds[2];
	
	//PHAChannelData<Cnts> default_ch_val;

	int nbins;


	protected:
	virtual EdgeSet<Etype> * getEdgeSet() {
		return this->Ebounds;
	};

	public:
	//Short-hand type for vector accessing
	typedef typename vector< pha_ref >::size_type size_type;
	
	friend class PHA2Writer;
	
	//Constructors
	PHA2() {
		this->nbins = 0;
		this->Ebounds = NULL;
		this->rmf_file = NULL;
		this->NullEntry = new chan_type( 0.0, 0.0, 0.0 );
		this->ogip_fields = new fields_class();
		TEdges = NULL;
	};
	
	//Destructor
	virtual ~PHA2() {
	
		#if PHADEBUG
		cout << "PHA2 DESTRUCTOR\n";
		#endif
	
	//	if (this->TEdges != NULL) delete this->TEdges;
	//	if (this->ogip_fields != NULL) delete this->ogip_fields;
	/*
		Since destructor is not virtual
		if (this->Ebounds != NULL) delete this->Ebounds;
		if (this->rmf_file != NULL) delete this->rmf_file;
		if (this->NullEntry != NULL) delete this->NullEntry;
	*/	
	//	this->clearData();
	};
	
	virtual bool IsDifferential() {
		return this->isDifferential;
	};
	virtual void IsDifferential(bool flag) {
		this->isDifferential = flag;
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(*phaitr)->IsDifferential(flag);
			++phaitr;
		}
	};
	
	
	//*************************** Energy Edge Set Handling **********************//
	/*
	Encapsulate the EdgeSet for the composite object by copying the data into a new private EdgeSet.
	SinglePHAs do not enforce encapsulation to allow simple flyweighting.
	
	Three methods:
	
	(public, virtual)
	1. UseBinning( EdgeSet *)
	2. UseBinning( EdgeSet *, bool )
		-virtualized from SinglePHA.  #2 is identical to #1,the second argument is meaningless.
		Gives the EdgeSet to use by pointer reference.  The data will be copied into a private member
		(method #3).  The private pointer is checked to make sure it doesn't point to EdgeSet*. This
		really shouldn't occur, but in case it does this bit of checking ensures it doesn't delete the 
		object it's supposed to copy. A new pointer is created to restore encapsulation.
				
	3. (private) UseBinning (EdgeSet&, bool delete_current)
		-Copies the data in EdgeSet& to a new object, then shares it with the child PHAs and sets the sharing bit.
	
	*/
	private:
	void UseBinning(EdgeSet<Etype>& edges, bool delete_existing) {	
	
		if (delete_existing && this->Ebounds != NULL) delete this->Ebounds;
		
		this->Ebounds = new EdgeSet<Etype>;
		*(this->Ebounds) = edges;
		
		//(this->Ebounds)->getBounds(tbounds);
	
		for (int i=0;i<nbins;i++) {
			datasets[i]->ShareBinning( this->Ebounds, bool(1) );
		}
	};
	
	public:
	
	void ForceUseBinning(EdgeSet<Etype> ** edges) {
		this->Ebounds = *edges;
	};
	
	virtual void UseBinning(EdgeSet<Etype> * edges, bool delete_old) {
		this->UseBinning(edges);
	};
	virtual void UseBinning(EdgeSet<Etype> * edges) {
		if (this->Ebounds != NULL && this->Ebounds == edges ) this->UseBinning(*edges,bool(0) );
		else this->UseBinning(*edges, bool(1) );
	};
	
	//****************************************************************************//
	
	friend class SpectrumModel;
	
	//**** Data Methods ****//
	virtual void UsePoissonErrors(bool flag) {
		this->pha_type::UsePoissonErrors(flag);
		for (int i=0;i<nbins;i++) {
			datasets[i]->UsePoissonErrors(flag);
		}
	};
	virtual bool UsePoissonErrors() {
		return this->pha_type::UsePoissonErrors();
	};
	
	virtual bool AddCount(int channel) {
	
		bool success = 1;
		for (int i=0;i<nbins;i++) {
			bool added = datasets[i]->pha_type::AddCount( channel );
			success = added && success;
		}	
		return success;
	};
	
	virtual bool AddCounts(int channel, Cnts counts) {
	
		bool success = 1;
		for (int i=0;i<nbins;i++) {
			bool added = datasets[i]->pha_type::AddCounts( channel, counts );
			success = added && success;
		}	
		#if PHA2_ADDCNT_DEBUG
		cout << "AddCounts to spectra? " << success << ": "<< channel << ", "<< counts <<"\n";
		#endif
		return success;
	};
	
	//**** Access OPERATORS *********//
	virtual pha_type& operator [](size_type i) {
		return *(datasets[i]);
	};
	virtual const pha_type& operator [](size_type i) const { 
		return *(datasets[i]);
	};
	
	pha_type& operator() (Ttype t) {
		int tbin = TEdges->findbin(t);
		
		TimeAccessError err;
		try {
			if (tbin < 0 || tbin >= this->nbins) throw &err;
		}
		catch (TimeAccessError * p) {
			cerr << "PHA2 Error, Time t="<<t<<" out of range.\n";
		};
		
		return *(datasets[tbin]);
	};
	
	pha_type& operator ()(size_type i, Ttype * t) {
		
		if (i >= 0 || i < this->nbins) *t = (*TEdges)[0][i];
		//else t = NULL;
		
		return *(datasets[i]);
	};
	
	pha_type& operator() (Ttype t,int& bin) {
		int tbin = TEdges->findbin(t);
		
		TimeAccessError err;
		try {
			if (tbin < 0 || tbin >= this->nbins) throw &err;
		}
		catch (TimeAccessError * p) {
			cerr << "PHA2 Error, Time t="<<t<<" out of range.\n";
		};
		
		bin = tbin;
		return *(datasets[tbin]);
	};
	
	virtual fields_class * getMiscFields() {
		if (this->ogip_fields == NULL) this->ogip_fields = new fields_class();
		return this->ogip_fields;
	};
	
	
	//**** Utilities for clearing, initializing *********//
	void clearData () {
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			delete *phaitr;
			++phaitr;
		}
		datasets.clear();
	};
	
	virtual void nullify() {
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(*phaitr)->nullify();
			++phaitr;
		}
	};
	
	virtual void reInit(int tbins) {
		
		this->clearData();
		datasets.reserve(tbins);

		if (TEdges != NULL && this->nbins != tbins) {
			TEdges->resize(2,tbins);
			this->nbins = tbins;
		}
	};
	
	virtual size_type size() {
		return datasets.size();
	};
	virtual size_type nchannels() {
		return this->Ebounds->size();
	};
		
	//**** Utilities for data sharing *********//
	virtual void SetExposure(Ttype exp) {
		this->SinglePHA<Cnts,Etype,Ttype>::SetExposure(exp);
		
		for (int i=0;i<nbins;i++) {
			datasets[i]->SetExposure(exp);
			(*TEdges)[1][i] = (*TEdges)[0][i] + exp;
		}
	};
	
	virtual void UseRmf_file(const char * file) {
		this->SinglePHA<Cnts,Etype,Ttype>::UseRmf_file(file);
		
		/*Uneccessary since the misc fields is shared
		
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(*phaitr)->ShareRmf_file( this->rmf_file );
			++phaitr;
		}
		*/
	};
	virtual void ShareRmf_file(char * file) {
		this->SinglePHA<Cnts,Etype,Ttype>::ShareRmf_file(file);
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(*phaitr)->ShareRmf_file( this->rmf_file );
			++phaitr;
		}
	};
	

	void AppendSpectrum(pha_type *& pha) {
		
		Ttype tstart;
		Ttype tstop;
		
		if (TEdges == NULL) {
			TEdges = new EdgeSet<Ttype>(0,0);
			tstart = 0;
			ogip_fields->tmin = tstart;
		} 
		if (this->nbins < 1) {
			tstart = 0;
			ogip_fields->tmin = tstart;
		} else {
			tstart = (*TEdges)[1][nbins-1];
		}
		
		tstop = tstart + pha->GetExposure();
		
		datasets.push_back( pha );
		TEdges->appendBin(tstart,tstop);
		
		
		ogip_fields->tmax = tstop;
		
		this->nbins+=1;
	};
	
	//Add a single spectrum to the lightcurve
	void AppendSpectrum(pha_type& pha) {
		
		if (TEdges == NULL) {
			this->AppendSpectrum(pha, 0);
			return;
		}
		
		pha_type * newpha = new pha_type;
		*newpha = pha;
		
		datasets.push_back( newpha );
		
		Ttype tstart = (*TEdges)[1][nbins-1];
		Ttype tstop = tstart + newpha->GetExposure();
		
		TEdges->appendBin(tstart,tstop);
		
		nbins++;
	};
	
	void AppendSpectrum(pha_type& pha, Ttype tstart) {
		
		pha_type * newpha = new pha_type;
		*newpha = pha;
		
		datasets.push_back( newpha );
		Ttype tstop = tstart + newpha->GetExposure();
		
		if (TEdges == NULL) {
			TEdges = new EdgeSet<Ttype>;
			(*TEdges)[0][0] = tstart;
			(*TEdges)[1][0] = tstop;
			ogip_fields->tmin = tstart;
			ogip_fields->tmax = tstop;
		} else {
			TEdges->appendBin(tstart,tstop);
			ogip_fields->tmax = tstop;
		}
		
		nbins++;
	};
	
	//Construct SinglePHAs with data sharing and uniform time binning.
	void BuildLightCurve(Ttype min, Ttype max, long tbins) {
		#if PHADEBUG
		cout << "PHA2::BuildLightCurve Begin\n";
		#endif
		this->clearData();
		datasets.reserve(tbins);
		
		if (TEdges == NULL) { 
			TEdges = new EdgeSet<Ttype>;
			TEdges->MakeLinearBins(min,max,tbins);
		};
		
		fields_class * metaData = this->getMiscFields();
		
		for (int t=0;t<tbins;t++)
		{
			pha_ref newpha = new pha_type(this->Ebounds, this->NullEntry, this->rmf_file);
			newpha->SetExposure( (*TEdges)[1][t] - (*TEdges)[0][t] );
			newpha->ShareFields(metaData);
			datasets.push_back( newpha );
		}
		
		this->nbins = tbins;
		
		#if PHADEBUG
		cout << "PHA2::BuildLightCurve End\n";
		#endif
	};
	
	//Math operators
	PHA2<Cnts,Etype,Ttype>& operator+=(Cnts scalar) {
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(**phaitr)+=scalar;
			++phaitr;
		}
		return (*this);
	};
	PHA2<Cnts,Etype,Ttype>& operator-=(Cnts scalar) {
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(**phaitr)-=scalar;
			++phaitr;
		}
		return (*this);
	};
	PHA2<Cnts,Etype,Ttype>& operator*=(Cnts scalar) {
		phaitr = datasets.begin();
		typename vector< pha_ref >::iterator end = datasets.end();
		while ( phaitr != end ) {
			(**phaitr)*=scalar;
			++phaitr;
		}
		return (*this);
	};
		
	PHA2<Cnts,Etype,Ttype>& operator=(pha_type& that) {
		#if PHADEBUG
		cout << "PHA2 = SinglePHA"<<endl;
		#endif
		
		if (this == &that) return *this;
		
		*(this->NullEntry) = that.copyNullChan();
		
		EdgeSet<Etype> newEb;
		that.getEdgeSet(newEb);
		
		this->UseBinning( newEb,1);
		
		return (*this);
	};

	
	friend ostream& operator<<(ostream &output, PHA2& obj) {
		int bins = obj.size();
		output << "# PHA2 Begin #####\n";
		output << "N Spectra: "<<bins<<"\n";

		if (bins == 0) {
			output << "\n";
			return output;
		}
		
		obj.phaitr = obj.datasets.begin();
		
		
		typename EdgeSet<Ttype>::size_type k;
		
		int nchan;
		
		for (k=0;k<bins;k++) {
			//sprintf(line, "
			output << "(" ;
			output << setw(9) << (*obj.TEdges)[0][k];
			output << " : ";
			output << setw(9) << (*obj.TEdges)[1][k];
			output << ") |";
			
			nchan = (*obj.phaitr)->size();
			for (int c=0;c<nchan;c++) {
				output << setw(11) << (**obj.phaitr)[c].cval << "|";
			}
			output <<  endl;
			++obj.phaitr;
		}
		//output << "\n";
		output << "# End ############\n\n";
		return output;
	};

};


class GBMResponse {
	private:
	DynMatrix<float> response;
	
	string file;
	fitsfile *fptr;
	int status;
	
	
	public:
	GBMResponse() {};
	GBMResponse(string file) {this->file = file;};
	bool loadMatrix();


};

#include "PHAStructures.cpp"

#endif
