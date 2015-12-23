#ifndef DYNMATRIX
#define DYNMATRIX

#include <vector>


template <typename T>
class DynMatrix
{

public:
	typedef typename std::vector<T>::size_type size_type;

private:
	size_type d0size;
	size_type d1size;
	
	bool variableLength;

protected:
  std::vector< std::vector<T> > data;

public:
	DynMatrix(){
		d0size = 0;
		d1size = 0;
		variableLength = 0;
	};
	DynMatrix(int rows, int cols)
	{
		for(int i=0; i<rows; ++i)
		{
		 data.push_back( std::vector<T>(cols) );
		}
		d0size = rows;
		d1size = cols;
		variableLength = 0;
	};
	inline ~DynMatrix(){};

	inline std::vector<T> & operator[](size_type i) { return data[i]; };
	inline const std::vector<T> & operator[] (size_type i) const { return data[i]; };
	inline DynMatrix<T>& operator=(DynMatrix<T>& that) {
		
		variableLength = that.variableLength;
		
		data.resize(that.d0size); 
		for(size_type i = 0; i < that.d0size; ++i) {
			data[i].resize( that.data[i].size() );
			data[i] = that.data[i];
		}
		
		d0size = that.d0size;
		d1size = that.d1size;
		
		return *this;
	};
	
	inline void useVariableLength(bool flag) {
		variableLength = flag;
	};
	inline bool isVariableLength() {
		return variableLength;
	};
	
	void size(size_type* rows, size_type* cols=NULL) {
		if (rows != NULL) *rows = d0size;
		if (cols != NULL) *cols = d1size;
	};
	
	inline void addRow() {
		data.push_back( std::vector<T>(d1size) );
		d0size++;
	};
	
	inline void addColumn(std::vector<T>& values) {
		
		for(size_type i = 0; i < values.size(); ++i) {
			data[i].push_back( values[i] );
		}
		//std::cout << "Added col: "<< values[0] << ", " << values[1] << std::endl;
		d1size++;
	};
	
	inline std::vector< T > pop_back() {
		d0size--;
		std::vector< T > back = data.back();
		data.pop_back();
		return back;
	};
	
	inline size_type rowSize(size_type i) {
		return data[i].size();
	};
	inline void resizeRow(size_type row, size_type newsize) {
		data[row].resize(newsize);
		useVariableLength(1);
	};

	inline void resize(size_type rows, size_type cols)
	{
		data.resize(rows);
		for(size_type i = 0; i < rows; ++i)
			data[i].resize(cols);
			
		d0size = rows;
		d1size = cols;
	};
	
//Fixed error where resize led to recycled references.
	inline void reinit(size_type rows, size_type cols, T revalue)
	{
		data.clear();
		for(size_type i = 0; i < rows; ++i) {
			std::vector<T> colvec(cols,revalue);
			data.push_back( colvec );
		}
			
		d0size = rows;
		d1size = cols;
	};
/*	virtual inline void reinit(int rows, int cols, T revalue)
	{
		data.clear();
		for(size_type i = 0; i < rows; ++i) {
			std::vector<T> colvec(cols,revalue);
			data.push_back( colvec );
		}
			
		d0size = rows;
		d1size = cols;
	};
*/

};

#endif