#ifndef MATRIXDEF
#define MATRIXDEF

/////////////////////////////////////////////////////////////////////////////
///////////// This header file defines the 'matrix' class /////////////
///////////// All matrices created have 'double' entries /////////////
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <string>
#include <iostream>
#include <cassert>


 using std::ostream;


 class matrix
 {

/////////////////// Variables of the matrix class /////////////////////////
//private:
public:
// Define dimensions of matrix
int rows, columns;
// Pointer to first entry of the vector of pointers
// (each of these pointers points to the first entry of each row)
double **xx;


/////////////////// Constructors //////////////////////////////////////////

public:
// Overwritten default constructor
matrix();

// Creates matrix of given dimension
matrix(int no_of_rows,int no_of_columns);

// Overwritten copy constructor
matrix(const matrix& A);

////////////////// Destructor /////////////////////////////////////////////
~matrix();


////////////////// Binary Operators ///////////////////////////////////////

friend matrix operator +(const matrix& A, const matrix& B);
friend matrix operator -(const matrix& A, const matrix& B);
friend matrix operator *(const double& p, const matrix& A);
friend matrix operator *(const matrix& A, const double& p);
friend matrix operator *(const matrix& A, const matrix& B);

friend matrix operator /(const matrix& A, const double& p);
friend matrix operator /(const matrix& b, const matrix& A);

////////////////////// Unary operators ////////////////////////////////////
friend matrix operator +(const matrix& A);
friend matrix operator -(const matrix& A);

// Overload ¬ to mean transpose
friend matrix t(const matrix& A);

/////////////////////// Other operators ///////////////////////////////////

// Overloads the assignment operator, '='
matrix& operator =(const matrix &v);

// Overloads (), so A(i,j) returns the i,j entry a la MATLAB
double &operator() (int i, int j);

// Returns the row dimension of the matrix
friend int rows(matrix A);
// Returns the column dimension of the matrix
friend int columns(matrix A);

////////////////////// Functions that are friends /////////////////////////
// Overloads the '<<' operator to allow easy printing of matrices
friend ostream& operator<<(ostream& output, const matrix& A);

// Create nxn Identity matrix
friend matrix eye(int size);

// Locates largest number below the diagonal, for matrix, A
friend int find_pivot(matrix A, int column);

friend matrix resize(matrix A, int m, int n);
};

matrix eye(int size);

// computes an nxn permutation matrix which swaps rows i and j
matrix permute_r(int n, int i, int j);

#endif
