#ifndef VECTORDEF
#define VECTORDEF

#include "matrix.h"


// 'vec' inherits from 'matrix' in a public way
class vec: public matrix
{
public:

//////////////////////// Constructors ////////////////////////////////////
// default constructor
vec();
// Constructor that takes 1 argument (size of the vec)
vec(int no_of_elements);

////////////////// Binary Operators //////////////////////////////////////

friend vec operator +(const vec& A, const vec& B);
friend vec operator -(const vec& A, const vec& B);

friend vec operator *(const double& p, const vec& A);
friend vec operator *(const vec& A, const double& p);
friend vec operator /(const vec& A, const double& p);


////////////////////// Unary operators //////////////////////////////////
friend vec operator +(const vec& A);
friend vec operator -(const vec& A);


/////////////////////// Other operators /////////////////////////////////

// Overloads (), so x(i) returns the ith entry a la MATLAB
double &operator() (int i);

// Overloads the assignment operator, '=' for vec RHS
vec& operator =(const vec &v);

////////////////////// Functions that are friends //////////////////////

// Default call is norm(v) and returns 2-norm
friend double norm(vec v, int p=2);

friend vec GMRES(matrix A, matrix b, matrix x0, double tol=1e-6);

friend vec GMRESout(matrix A, matrix b, matrix x0, double tol=1e-6);

//friend vec resize(vec v, int m);

};

vec mat2vec(matrix A);

#endif
