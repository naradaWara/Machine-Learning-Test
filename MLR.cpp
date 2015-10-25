#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <fstream>
#include <math.h>
using namespace std;

//-------------------------------------
// This class handles Matrix operations
// including solving  a linear system
//-------------------------------------
class CMatrix
{	
private:
	int m_rows;
	int m_cols;
	char m_name[128];
	CMatrix();
public:
	double **m_pData;
	
	//----------------------
	//Constructor
	//-----------------------
	CMatrix(const char *name, int rows, int cols) : m_rows(rows), m_cols(cols)
	{
		m_pData = new double*[m_rows];
		for(int i = 0; i < m_rows; i++)
			m_pData[i] = new double[m_cols];
		for(int i = 0; i < m_rows; i++)
		{
			for(int j = 0; j < m_cols; j++)
			{
				m_pData[i][j] = 0.0;
			}
		}
	}
	
	//------------------
	//Copy constructor 
	//------------------
	CMatrix(const CMatrix &other)
	{
		//strcpy(m_name, other.m_name);
		m_rows = other.m_rows;
		m_cols = other.m_cols;
		m_pData = new double*[m_rows];
		for(int i = 0; i < m_rows; i++)
			m_pData[i] = new double[m_cols];
		for(int i = 0; i < m_rows; i++)
		{
			for(int j = 0; j < m_cols; j++)
			{
				m_pData[i][j] = other.m_pData[i][j];
			}
		}
	}
	
	//------------
	//Destructor
	//------------
	~CMatrix()
	{
		for(int i = 0; i < m_rows; i++)
			delete [] m_pData[i];
		delete [] m_pData;
		m_rows = m_cols = 0;
	}
	
	
	//--------------------
	//Overloaded operator
	//--------------------
	CMatrix& operator = (const CMatrix &other)
	{
		if( this->m_rows != other.m_rows ||
			this->m_cols != other.m_cols)
		{
			std::cout << "WARNING: Assignment is taking place with by changing the number of rows and columns of the matrix";
		}
		for(int i = 0; i < m_rows; i++)
			delete [] m_pData[i];
		delete [] m_pData;
		m_rows = m_cols = 0;
		//strcpy(m_name, other.m_name);
		m_rows = other.m_rows;
		m_cols = other.m_cols;
		m_pData = new double*[m_rows];
		for(int i = 0; i < m_rows; i++)
			m_pData[i] = new double[m_cols];
		for(int i = 0; i < m_rows; i++)
		{
			for(int j = 0; j < m_cols; j++)
			{
				m_pData[i][j] = other.m_pData[i][j];
			}
		}
		return *this;
	}
	
	//-----------------------
	//Matrix transpose
	//-----------------------
	CMatrix Transpose()
	{
		CMatrix trans("TR", m_cols, m_rows);
		for(int i = 0; i < m_rows; i++)
		{
			for(int j = 0; j < m_cols; j++)
			{
				trans.m_pData[j][i] = m_pData[i][j];
			}
		}
		return trans;
	}
	
	
	//-------------------------------
	// Multiplication operator
	//------------------------------
	CMatrix operator * (const CMatrix &other)
	{
		if( this->m_cols != other.m_rows)
		{
			std::cout << "Multiplication could not take place because number of columns of 1st Matrix and number of rows in 2nd Matrix are different";
			return *this;
		}
		CMatrix result("", this->m_rows, other.m_cols);
		for(int i = 0; i < this->m_rows; i++)
		{
			for(int j = 0; j < other.m_cols; j++)
			{
				for(int k = 0; k < this->m_cols; k++)
				{
					result.m_pData[i][j] += this->m_pData[i][k] * other.m_pData[k][j];
				}
			}
		}
		return result;
	}
	
	//------------------------------------
	// Merges two matrices columnwise
	// 'in' Matrix will be placed to the tight of this Matrix
	//--------------------------------------
	CMatrix merge(CMatrix& in){
		CMatrix ret("MER", m_rows, m_cols+ in.m_cols);
		for (int i=0;i< m_rows;i++){
			for (int j=0;j< m_cols+in.m_cols;j++){
				if (j<m_cols)
				   ret.m_pData[i][j] = m_pData[i][j];
			        else
			           ret.m_pData[i][j] = in.m_pData[i][j-m_cols];	
			}
		}
		return ret;
	}
	
	//-----------------------------------------------------------------
	// Gauss elimination routine for solving a linear system Ax=b
	// this Matrix holds the [A | b], where A is nxn Matrix, b is nx1 vector
	// Gauss elimination returns x which satisfies Ax=b
	//---------------------------------------------------------------------
	CMatrix  gauss() {
		double **A= m_pData;
		int n = m_rows;
		CMatrix x ("soln", n,1);
		for (int i=0; i<n; i++) {
			// Search for the maximum value in the ith column
			double maxEl = abs(A[i][i]);
			int maxRow = i;
			for (int k=i+1; k<n; k++) {
				if (abs(A[k][i]) > maxEl) {
					maxEl = abs(A[k][i]);
					maxRow = k;
				}
			}
			// Swap the row containing the maximum value with the current row 
			for (int k=i; k<n+1;k++) {
				double tmp = A[maxRow][k];
				A[maxRow][k] = A[i][k];
				A[i][k] = tmp;
			}
			// Make all rows below this one has zero value in current column
			for (int k=i+1; k<n; k++) {
				if (A[i][i] == 0){
					cout << "This system of equations does not have a unique solution" << endl;
					exit(0);
				}
				double c = -A[k][i]/A[i][i];
				for (int j=i; j<n+1; j++) {
					if (i==j) {
						A[k][j] = 0;
					} else {
						A[k][j] += c * A[i][j];
					}
				}
			}
		}
		// Solve equation Ax=b for x using back substitution
		for (int i=n-1; i>=0; i--) {
			if (A[i][i] == 0){
					cout << "This system of equations does not have a unique solution" << endl;
					exit(0);
			}
			x.m_pData[i][0] = A[i][n]/A[i][i];
			for (int k=i-1;k>=0; k--) {
				A[k][n] -= A[k][i] * x.m_pData[i][0];
			}
		}
		
		return x;
	}
	
	
	friend std::istream& operator >> (std::istream &is, CMatrix &m);
	friend std::ostream& operator << (std::ostream &os, const CMatrix &m);   
};

//---------------------------------------------------------
// Overloaded output operator
//---------------------------------------------------------
std::ostream& operator << (std::ostream &os,const CMatrix &m)
{
	for(int i = 0; i < m.m_rows; i++)	
	{
		for(int j = 0; j < m.m_cols; j++)	
		{
			char buf[32];
			double data = m.m_pData[i][j];
			if( m.m_pData[i][j] > -0.00001 &&
				m.m_pData[i][j] < 0.00001)
			data = 0;
			sprintf(buf, "%10.2lf ", data);
			os <<  buf;
		}
		os << "\n";	
	}
	os << "\n\n";
	return os;
};



//-----------------------------------------------------------------------
// This is the class which handles Multiple Linear Regression.
// It basically reads input data, calculates the regression parameters
// based on the training data,computes the target values for the test data
// and finally outputs them.
//----------------------------------------------------------------------
class mlr{
	
private:	

	CMatrix *fmat; //feature matrix, training data 
	CMatrix *tmat; // corresponding output vector for training data  
	CMatrix *testMat; // feature matrix, test data

public:
	
	//We are satisfied with the default constructor, hence no constructor
	
	//------------------------------
	// Destructor
	//------------------------------
	~mlr(){
		delete fmat;
		delete tmat;
		delete testMat;
	}
	
	//------------------
	//Allocate fmat
	//------------------
	void create_fmat(const char * name, int r, int c){
	  fmat = new CMatrix(name,r,c);
	}
	
	//-----------------------
	//Allocate tmat
	//-----------------------
	void create_tmat(const char * name, int r, int c){
	  tmat = new CMatrix(name,r,c);
	}
	
	//------------------------
	// Allocate testMat
	//-------------------------
	void create_testMat(const char *name, int r , int c){
		testMat = new CMatrix(name, r,c);
	}
	
	//-----------------------------------------------------------------
	// We do the regression here and compute the output for test data 
	//----------------------------------------------------------------
	void fit_data(){
		CMatrix X_t = (fmat->Transpose()); //X' = transpose of fmat
		CMatrix X_t_X = X_t*(*fmat);       //X' * X
		
		CMatrix X_t_target = X_t*(*tmat);  //X' * t  (t= training outputs)
		
		//Create [X'*X | X'*t] for gauss elimination routine
		CMatrix merged = X_t_X.merge(X_t_target); 
		
		// run Gauss elimination with the merged Matrix, and get the 
		// regression parameters weights
		CMatrix weights = merged.gauss();
		
		//Finally compute the outputs for the test data, i.e. predicted values
		CMatrix predicted = (*testMat)*weights;
		
		//output the predicted values
	     	cout << predicted;
	}

	
	//-----------------------------------------------------------------
	// Routine for reading the input data from the STDIN (or any istream
	//------------------------------------------------------------------
	void read_data(istream &is)
	{
		int dim_features; // feature dimension 
		int n_samples; // number of samples 
		int n_test_samples; // number of test samples 
		int d1, d2, d3; //temp values
		
		//read the feature dimension and number of training samples
		// and check whether they are in the correct range
		is >> d1; 
		is >> d2;
		dim_features=d1;
		n_samples= d2;
		if ((dim_features > 10) || (dim_features < 1)) {
			cout << "feature dimension =" << dim_features  << "not in the allowed range" << endl;
			exit(0);
		} 
		if ((n_samples > 100) || (n_samples < 5)) {
			cout << "number of samples =" << n_samples <<  "not in the allowed range" << endl;
			exit(0);
		} 
		
		//Now we can allocate for fmat and tmat
		create_fmat("F", n_samples, dim_features+1); 
		create_tmat("T", n_samples, 1);
		
		float value=0;
		
		//account for unread characters in the buffer
		is.clear();
		is.ignore(2056, '\n');
		
		//Now read the training samples and fill fmat and tmat
		for (int k=0;k<n_samples;k++){
			
			fmat->m_pData[k][0]=1.0; //first column is 1 corresponding to the bias coefficient
			for (int l=1;l<=dim_features;l++){
				is >> value;
				if ((value <0) || (value > 1.0)){
					cout << "Feature value = " << value << " not in the allowed range" << endl;
					exit(0);
				}
				fmat->m_pData[k][l] = value;
			}
			is >> value; 
			if ((value <0) || (value > 1.0E6)){
				cout << "Target value = " << value << " not in the allowed range" << endl;
				exit(0);
			}
			tmat->m_pData[k][0]= value;
		} 
		
		// Now read the number of test samples and check whether it is in the correct range
		is >> n_test_samples;
		if ((n_test_samples > 100) || (n_test_samples < 1)) {
			cout << "number of test samples =" << n_test_samples <<  "not in the allowed range" << endl;
			exit(0);
		} 
		
		//Allocate the testMat for testdata
		create_testMat("T", n_test_samples, dim_features+1);
		is.clear();
		is.ignore(2056, '\n');
		
		//And read the test samples
		for (int k=0;k<n_test_samples;k++){
			testMat->m_pData[k][0]=1.0; //first column is 1 corresponding to the bias coefficient
			for (int l=1;l<=dim_features;l++){
				is >> value;
				if ((value <0) || (value > 1.0)){
					cout << "Feature value = " << value << " not in the allowed range" << endl;
					exit(0);
				}
				testMat->m_pData[k][l] = value;
			}
		}
		
	}
};


//-----------------------------------------------------
//  The main program
//---------------------------------------------------
int main() {
	/* Enter your code here. Read input from STDIN. Print output to STDOUT */
	mlr MLR;  //Create the MLR object
/*	
	ifstream filein;
	filein.open("data2.txt");
	MLR.read_data(filein);
*/	
	MLR.read_data(cin);  // Read input data
	MLR.fit_data();      // run regression and generate the output

	return 0;
}