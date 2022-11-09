//#include <Rcpp.h>
//using namespace Rcpp;

#include "vector.h" // include the header file
//////////////////////// Constructors /////////////////////////////////////
// in the class, 'vec' there is a constructor of the same name
vec::vec()
// runs the default matrix constructor
  :matrix(){
}

vec::vec(int no_of_elements)

  :matrix(no_of_elements,1)
  {
  
  }

/////////////////////// Other operators //////////////////////////////////
// Overloads (), so x(i) returns the ith entry a la MATLAB
double &vec::operator() (int i){ // Can call or assign values.
    if (i < 1)
    {
    std::cout << "Error: Your index may be too small \n\n";
    }else if (i > rows){
        std::cout << "Error: Your index may be too large \n\n";
        }
    return xx[i - 1][0];
    }

// Operator returns a matrix equal to the RHS
vec& vec::operator =(const vec &v){
    // Destruct previous entries
    for (int i=0; i < rows; i++)
    {
    delete[] xx[i];
    }
    delete[] xx;
    // Assign new dimensions to be equal to that of the RHS
    rows = v.rows;
    columns = v.columns;
    // Allocate the memory as in the constructor
    xx = new double *[v.rows];
    
    for (int i=0; i < v.rows ; i++){
        xx[i] = new double[v.columns];
        }
    
    // Copy the values across from the RHS
    for (int i=0; i < v.rows; i++){
        for (int j=0; j < v.columns ; j++)
        {
        // Set entries to be the same as the RHS matrix
        xx[i][j]=v.xx[i][j];
        }
        }
    return *this;
    
    }

////////////////// Binary Operators ///////////////////////////////////////
vec operator +(const vec& A, const vec& B){
    int m,n;
    m = rows(A);
    n = rows(B);
    if ( m != n ){
        std::cout << "Error: Matrices of different dimensions.";
        std::cout << " Returned first argument";
        return A;
        }else{
        vec v(m);
        for (int i=0; i<m; i++){
            v(i+1) = A.xx[i][0]+B.xx[i][0];
            }
        return v;
        }
    }

vec operator -(const vec& A, const vec& B){
  int m,n;
  m = rows(A);
  n = rows(B);
    if ( m != n ){
        std::cout << "Error: Matrices of different dimensions.";
        std::cout << " Returned first argument";
        return A;
        }else{
        vec v(m);
        for (int i=0; i<m; i++){
            v(i+1) = A.xx[i][0] - B.xx[i][0];
            }
        return v;
        }
    }

vec operator *(const double& p, const vec& A){
    int m = rows(A);
    vec v(m);
    for (int i=0; i<m; i++){
        v(i+1) = p*A.xx[i][0];
        }
    return v;
    }
//
vec operator *(const vec& A, const double& p){
    int m = rows(A);
    vec v(m);
    for (int i=0; i<m; i++){
        v(i+1) = p*A.xx[i][0];
        }
    return v;
    }

vec operator /(const vec& A, const double& p){
    int m = rows(A);
    vec v(m);
    for (int i=0; i<m; i++){
        v(i+1) = A.xx[i][0]/p;
        }
    return v;
    }

////////////////////// Unary operators /////////////////////////////////
vec operator +(const vec& A){
    int m = rows(A);
    vec v(m);
    for (int i=0; i<m; i++){
        v(i+1) = A.xx[i][0];
        }
    return v;
    }


vec operator -(const vec& A){
    int m = rows(A);
    vec v(m);
    for (int i=0; i<m; i++){
        v(i+1) = -A.xx[i][0];
      }
    return v;
    }

////////////////////// Functions that are friends ////////////////////////
// Function that returns the first column of a matrix, A as a vec
vec mat2vec(matrix A){
  std::cout << "mat2vec" << std::endl;
// create vec with same no of rows as A and 1 column
vec v(rows(A));

for (int i=1; i <= rows(A); i++){
        v(i)=A(i,1); // copy only first column
        }
    
    return v;
    }

double norm(vec v, int p){
    // define variables and initialise sum
    double temp, value, sum = 0.0;
    for (int i=1; i <= rows(v); i++){
        // floating point absolute value
        temp = fabs(v(i));
        sum += pow(temp, p);
        }
    
    value = pow(sum, 1.0/((double)(p)));
    return value;
    }
