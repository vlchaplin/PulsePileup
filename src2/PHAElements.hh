/*
 *  PHAElements.hh
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef PHAELEMENTS
#define PHAELEMENTS

#define TRUNC_DETNAME_HISTORY_ATLENGTH 20
#define TRUNC_DETNAME_STRING "..."

//#define DEBUG_FIELD_OPERATION

#include <string>
#include <stdexcept>
#include <iostream>
using namespace std;

template<typename C>
class PHAChannelData {

	private:
	class ChannelAccesError : public std::runtime_error {
		public:
		ChannelAccesError() : std::runtime_error("ChannelAccessError") {
		
		};
	};

	public:
	
	typedef float err_type;
	
	C cval;
	err_type stat_err;
	err_type sys_err;
	
	PHAChannelData() {
		stat_err = 0.0;
		sys_err = 0.0;
	};
	PHAChannelData(C c,err_type err1, err_type err2) {
		cval = c;
		stat_err = err1;
		sys_err = err2;
	};
	~PHAChannelData() {};
	
	C& operator[](err_type channel_err) {
		stat_err = channel_err;
		return cval;
	};
	
	C& operator[](err_type * channel_err) {
	
		ChannelAccesError err;
		
		try {
			if(1) {
				*channel_err = stat_err;
				return cval;
			}
			throw &err;
		}
		catch (ChannelAccesError *p) {
			cerr << "ERROR: Channel value parameter access must be float *\n\n";
		}
	};
    
    PHAChannelData<C>& operator+=(PHAChannelData<C>& channel) {
		this->cval += channel.cval;
		this->stat_err = sqrt( pow(this->stat_err,2) + pow(channel.stat_err,2) );
		this->sys_err = sqrt( pow(this->sys_err,2) + pow(channel.sys_err,2) );
		return (*this);
	};
	PHAChannelData<C>& operator-=(PHAChannelData<C>& channel) {
		this->cval -= channel.cval;
		this->stat_err = sqrt( pow(this->stat_err,2) + pow(channel.stat_err,2) );
		this->sys_err = sqrt( pow(this->sys_err,2) + pow(channel.sys_err,2) );
		return (*this);
	};
    PHAChannelData<C>& operator*=(PHAChannelData<C>& channel) {
		this->cval *= channel.cval;
		this->stat_err = sqrt( pow(this->stat_err,2) + pow(channel.stat_err,2) );
		this->sys_err = sqrt( pow(this->sys_err,2) + pow(channel.sys_err,2) );
		return (*this);
	};
	PHAChannelData<C>& operator/=(PHAChannelData<C>& channel) {
		this->cval /= channel.cval;
		this->stat_err = sqrt( pow(this->stat_err,2) + pow(channel.stat_err,2) );
		this->sys_err = sqrt( pow(this->sys_err,2) + pow(channel.sys_err,2) );
		return (*this);
	};
	PHAChannelData<C>& operator+=(C& scalar) {
		this->cval += scalar;
		return *this;
    };
	PHAChannelData<C>& operator-=(C& scalar) {
		this->cval -= scalar;
		return *this;
    };
    PHAChannelData<C>& operator*=(C& scalar) {
		this->cval *= scalar;
		return *this;
    };
	PHAChannelData<C>& operator/=(C& scalar) {
		this->cval /= scalar;
		return *this;
    };
	PHAChannelData<C>& operator=(PHAChannelData<C>& channel) {
		this->cval = channel.cval;
		this->stat_err = channel.stat_err;
		this->sys_err = channel.sys_err;
		return (*this);
	};
	const PHAChannelData<C>& operator=(const PHAChannelData<C>& channel) {
		this->cval = channel.cval;
		this->stat_err = channel.stat_err;
		this->sys_err = channel.sys_err;
		return (*this);
	};

};

template<typename T>
class OGIP_PHA_misc {
	
	protected:
	string * detname;
	string keyterminator;
	size_t kyt_size;
	
	void NamedOperationsLogging(string& operation, OGIP_PHA_misc& operand ) {
		if (operand.detname != NULL) {
			if (this->detname != NULL) {
				size_t length = this->detname->length();
				#ifdef DEBUG_FIELD_OPERATION
				cout << "SUBSTR: "<<this->detname->substr(length - kyt_size, kyt_size) << ", " << this->keyterminator<<"\n";
				cout << "COMPARE: "<<this->detname->compare( length - kyt_size, kyt_size, this->keyterminator )<<"\n";
				#endif
				if ( length <= TRUNC_DETNAME_HISTORY_ATLENGTH) {
					(*this->detname) += (operation + *operand.detname);
				} else if ( this->detname->compare( length - kyt_size, kyt_size, this->keyterminator ) != 0 ) {
					(*this->detname) += this->keyterminator;
				}
			} 
			else this->detname = new string( operation + *operand.detname );
		};
		
		#ifdef DEBUG_FIELD_OPERATION
			cout << "DETNAME: "<<*this->detname<<"\n";
		#endif
	};
	
	bool shared;

	public:
	
	unsigned long int refNum;
	
	float areasc;
	float backsc;
	unsigned short qual;
	unsigned short grpg;
	string response;
	string bakfile;
	string object;

	char * charDetName;
	char * charRmfFile;

	T tzero;
	T tmin;
	T tmax;
	

	OGIP_PHA_misc() {
		areasc = 1.0;
		backsc = 1.0;
		tzero = 0.0;
		tmin = 0.0;
		tmax = 0.0;
		qual = 0;
		grpg = 0;
		detname = NULL;
		keyterminator = TRUNC_DETNAME_STRING;
		kyt_size = keyterminator.length();
		shared = 0;
	};
	~OGIP_PHA_misc(){
		if (this->detname != NULL) delete detname;
	};
	
	void setShared() { this->shared = 1; };
	bool isShared() { return this->shared; };
	
	inline void setDetName(string& detname) {
		if (this->detname == NULL) this->detname = new string;
		this->detname->assign(detname);
	};
	inline bool getDetName(string& copy) {
		if (this->detname == NULL) return 0;
		copy.assign( *(this->detname) );
		return 1;
	};
	inline bool getDetName(char * copy) {
		if (this->detname == NULL) return 0;
		strcpy( copy, this->detname->c_str() );
		return 1;
	};
	inline void setResponse(const char * rspFile) {
		this->response.assign(rspFile);
	};
	inline void setResponse(string& rspFile) {
		this->response.assign(rspFile);
	};
	inline string getResponse() {
		return this->response;
	};
	
	OGIP_PHA_misc& operator=(OGIP_PHA_misc& that) {
		this->areasc = that.areasc;
		this->backsc = that.backsc;
		this->qual = that.qual;
		this->grpg = that.grpg;
		if (that.detname != NULL) this->setDetName( *that.detname );
		this->response = that.response;
		
		return *this;
	};
	
	OGIP_PHA_misc& operator+=(OGIP_PHA_misc& that) {
		if (that.detname != NULL) {
			string operation("+");
			this->NamedOperationsLogging( operation, that );
		}
		return *this;
	};
	OGIP_PHA_misc& operator-=(OGIP_PHA_misc& that) {
		if (that.detname != NULL) {
			string operation("-");
			this->NamedOperationsLogging( operation, that );
		}
		return *this;
	};
	OGIP_PHA_misc& operator*=(OGIP_PHA_misc& that) {
		if (that.detname != NULL) {
			string operation("*");
			this->NamedOperationsLogging( operation, that );
		}
		return *this;
	};
};

#endif
