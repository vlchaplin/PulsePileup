/*
 *  VectorSpan.h
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */
 
 /* EXAMPLES: ****
 
 
 	int data[] = {0,1,2,3,4,5,6,7,8,9};
	VectorSpan< int, int* > intsp;
	intsp.Insert( data, data+5 );
	intsp.Insert( data+2, data+5 );
	VectorSpan< int, int* >::iterator iit = intsp.begin();
	while (iit != intsp.end()) {
		int * tstart = (*iit)->t_start;
		int * tend = (*iit)->t_end;
		while (tstart != tend) {
			cout << *(tstart++) << " ";
		}
		cout << endl;
		iit++;
	}
	
	
	char tempSubstr[30];
	char test[] = "AAbbCCddEEffGGhh";

	VectorSpan< char, char* > strsp;

	strsp.Insert( test, test+2 );
	strsp.Insert( test+2, test+4 );
	strsp.Insert( test+6, test+8 );
	
	cout << "Total chars: " << strsp.GetSpanSize() << endl;
	
	VectorSpan< char, char* >::iterator it = strsp.begin();
	while (it != strsp.end()) {
		size_t nchars = (*it)->t_end - (*it)->t_start;
		
		strncpy( tempSubstr, (*it)->t_start, nchars );
		tempSubstr[nchars]='\0';
		it++;
		
		cout << tempSubstr << endl;
	}
	size_t ncharTot = strsp.GetSpanSize();
	strsp.ExtractCopy( tempSubstr );
	tempSubstr[ncharTot] = '\0';
	cout << tempSubstr << endl;
 
 *******/


#ifndef VECTORSPAN_HH
#define VECTORSPAN_HH

#include <stddef.h>
#include <vector>
#include <list>

using namespace std;

template<typename vtype, typename t_iterator>
class VectorSpan {

	public:
	
	class range_slice {
		public:
		t_iterator t_start;
		t_iterator t_end;
		range_slice(){};
		range_slice(t_iterator& st, t_iterator& en) 
			: t_start(st), t_end(en) {};
		size_t diff() {
			size_t difference=0;
			t_iterator temp = t_start;
			while (temp != t_end) {
				difference++;
				temp++;
			}
			return difference;
		};
	};
	
	typedef vtype data_type;
	typedef t_iterator data_iterator_type;
	typedef typename vector<range_slice *>::iterator iterator;
	typedef typename vector<range_slice *>::size_type size_type;

	private:
	vector<range_slice *> span;
	typename vector<range_slice *>::iterator itr;
	
	public:


	VectorSpan() {};
	~VectorSpan() {
		this->clear();
	};
	
	void index(vector< vtype >& data,  VectorSpan< size_t, size_t >& spanIdx, long shift=0 ) {
	
		typename vector< vtype >::iterator dataP1;
		typename vector< vtype >::iterator dataP2;
		
		VectorSpan< size_t, size_t >::iterator idxtr = spanIdx.begin();
		size_t maxNum = data.size();
		while ( idxtr != spanIdx.end() ) {
		
			dataP1 = data.begin();
			dataP2 = data.begin();
			
			
			if ( (*idxtr)->t_start+shift < maxNum ) advance( dataP1, (*idxtr)->t_start+shift ); 
			else {
				idxtr++;
				continue;
			}
			
			if ( (*idxtr)->t_end+shift+1 < maxNum ) advance( dataP2, (*idxtr)->t_end+shift );
			else dataP2 = data.end();
			
			span.push_back( new range_slice( dataP1, dataP2 ) );
			idxtr++;
		}

	};
	
	void Insert( t_iterator intrvl_begin, t_iterator intrvl_end ) {
		span.push_back( new range_slice( intrvl_begin, intrvl_end ) );
	};
	
	list< t_iterator > List() {
		list< t_iterator > temp;
		t_iterator temp_t_itr;
		itr = span.begin();
		while (itr != span.end()) {
			
			temp_t_itr = (*itr)->t_start;
			
			while (temp_t_itr != (*itr)->t_end) temp.push_back( temp_t_itr++ );
			
			itr++; 
		}
		return temp;
	};
	void Enumerate(vector< t_iterator >& ref) {
		t_iterator temp_t_itr;
		itr = span.begin();
		while (itr != span.end() ) {
			temp_t_itr = (*itr)->t_start;
			while (temp_t_itr != (*itr)->t_end) {
				ref.push_back( temp_t_itr++ );
			}
			itr++; 
		}
	};
	vector< t_iterator > Enumerate() {
		vector< t_iterator > temp;
		t_iterator temp_t_itr;
		itr = span.begin();
		
		while (itr != span.end() ) {
			temp_t_itr = (*itr)->t_start;
			
			while (temp_t_itr != (*itr)->t_end) {
				temp.push_back( temp_t_itr++ );
			}
			itr++; 
		}

		return temp;
	};
	vector< vtype > ExtractCopy() {
		vector< vtype > temp;
		t_iterator temp_t_itr;
		itr = span.begin();
		while (itr != span.end()) {
			
			temp_t_itr = (*itr)->t_start;
			
			while (temp_t_itr != (*itr)->t_end) {
				temp.push_back( *temp_t_itr++ );
			}
			
			itr++; 
		}
		return temp;
	};
	void ExtractCopy(vtype * copy) {
		
		static size_t ellsz = sizeof(vtype);
		vtype * v_itr;
		
		t_iterator temp_t_itr;
		
		v_itr = copy;
		itr = span.begin();
		while (itr != span.end()) {
			
			temp_t_itr = (*itr)->t_start;
			
			while (temp_t_itr != (*itr)->t_end) {
				*v_itr = *temp_t_itr++;
				v_itr += ellsz;
			}
			
			itr++; 
		}
	};
	size_t GetSpanSize() {
		size_t countSpan=0;
		t_iterator temp_t_itr;
		itr = span.begin();
		while (itr != span.end()) {
			countSpan += (*itr)->diff();
			itr++; 
		}
		return countSpan;
	};
	
	iterator begin() {
		return span.begin();
	};
	iterator end() {
		return span.end();
	};
	
	size_t nslices() {
		return span.size();
	}

	void clear() {
		itr = span.begin();
		while (itr != span.end()) { 
			delete *itr;
			itr++;
		}
		span.clear();
	};


};




#endif