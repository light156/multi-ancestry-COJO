#include "tcojo.h"


double median(const ArrayXd &eigen_vector)
{
    int size = eigen_vector.size();
    vector<double> b(eigen_vector.data(), eigen_vector.data() + size);
    double b_median; 

    sort(b.begin(), b.end());
    if (size%2==1)
        b_median = b[(size-1)/2];
    else 
        b_median = (b[size/2]+b[size/2-1])/2;

    vector<double>().swap(b);
    return b_median;
}


void to_upper(string &str)
{
	int i=0;
	for(i=0; i<str.size(); i++){
		if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
	}
}


int split_string(const string &str, vector<string> &vec_str, string separator)
{
	if(str.empty()) return 0;
	vec_str.clear();

	int i=0;
	bool look=false;
	string str_buf;
	string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
	string::size_type pos;

	for(i=0; i<separator.size(); i++){
		pos=symbol_pool.find(separator[i]);
		if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
	}

	for(i=0; i<str.size(); i++){
		if( symbol_pool.find(str[i])!=string::npos ){
			if(!look) look=true;
			str_buf += str[i];
		}
		else{
			if(look){
				look=false;
				vec_str.push_back(str_buf);
				str_buf.erase(str_buf.begin(), str_buf.end());
			}
		}
	}
	if(look) vec_str.push_back(str_buf);

	return vec_str.size();
}


void append_row(ArrayXXd &matrix, const ArrayXXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void append_row(MatrixXd &matrix, const MatrixXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void append_column(MatrixXd &matrix, const MatrixXd &vector)
{
    int numCols = matrix.cols();
    matrix.conservativeResize(NoChange, numCols+1);
    matrix.col(numCols) = vector;
}


void remove_row(ArrayXXd &matrix, int index)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void remove_row(MatrixXd &matrix, int index)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void remove_column(MatrixXd &matrix, int index)
{   
    // -1 indicates the last column
    int numRows = matrix.rows(), numCols = matrix.cols()-1;

    if (index != -1 && index < numCols)
        matrix.middleCols(index, numCols-index) = matrix.rightCols(numCols-index).eval();

    matrix.conservativeResize(NoChange, numCols);
}