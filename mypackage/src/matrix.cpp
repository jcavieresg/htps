//#include <Rcpp.h>
//using namespace Rcpp;

#include "matrix.h"
#include "vector.h"

/////////////////// Constructors //////////////////////////////////////////

// Constructor that overrides compiler generated default constructor
matrix::matrix()
{
// constructs an empty object with rows = columns = 0
rows = 0;
columns = 0;
xx = NULL;
}

// Constructor for basic matrix with specified dimensions
matrix::matrix(int no_of_rows, int no_of_columns){
// Matrix dimension fields
rows = no_of_rows;
columns = no_of_columns;

// A is an array of m pointers, each pointing to the
// first entry in the vector (of length, 'no of columns')
xx = new double *[no_of_rows];

// Allocate the memory for the entries of the matrix
for (int i=0; i<no_of_rows; i++){

// Creates 'no of rows' rows of length, 'no of columns'
xx[i] = new double[no_of_columns];
}

for (int i=0; i<no_of_rows; i++){
for (int j=0; j<no_of_columns; j++){
// Initialise the entries of the matrix to zero
      xx[i][j]=0.0;
    }
   }
}

// Copy constructor ??? creates matrix with the same entries as input, A
matrix::matrix(const matrix& A){
rows = A.rows; // Matrix dimension fields
columns = A.columns;

// A is an array of m pointers, each pointing to the
// first entry in the vector (of length, 'A.columns')
xx = new double *[A.rows];

// Allocate the memory for the entries of the matrix
for (int i=0; i<A.rows; i++){
// Creates 'A.rows' rows of length, 'A.columns'
xx[i] = new double[A.columns];
}

for (int i=0; i < A.rows; i++){
for (int j=0; j < A.columns; j++){
// Copy across the entries from matrix A
xx[i][j]=A.xx[i][j];
}
}
}

// Destructor
matrix::~matrix(){
 if ( rows > 0 || columns > 0 ){
 for (int i=0; i < rows; i++){
 delete[] xx[i];
    }
delete[] xx;
 }
}

///////////////////// Binary Operators ///////////////////////////////////

// Overload the + operator to evaluate: A + B, where A and B are matrices
matrix operator +(const matrix& A, const matrix& B){
int m = rows(A), n = columns(A), p = rows(B), q = columns(B);
if ( m != p || n !=q ){
        std::cout << "Error: Matrices of different dimensions";
        std::cout << "Returned first argument";
        return A;
}else{
matrix C(m,n);
for (int i=0; i<m; i++){
for (int j=0; j<n; j++){
C.xx[i][j]=A.xx[i][j]+B.xx[i][j];
}
}
return C;
}
}


// Overload the ??? operator to evaluate: A ??? B, where A and B are matrices
matrix operator -(const matrix& A, const matrix& B)
{
int m,n,p,q;
m = rows(A);
n = columns(A);
p = rows(B);
q = columns(B);
if ( m != p || n !=q ){
std::cout << "Error: Matrices of different dimensions";
std::cout << "Returned first argument";
return A;
}else {
matrix C(m,n);
for (int i=0; i<m; i++){
for (int j=0; j<n; j++){
C.xx[i][j]=A.xx[i][j] - B.xx[i][j];
}
}
return C;
}
}

// Definition of multiplication between a scalar, p and a matrix, A
matrix operator *(const double& p, const matrix& A){
    
// Create a matrix with the same dimensions as A
matrix B(A.rows,A.columns);
for (int i=0; i<A.rows; i++){
for (int j=0; j<A.columns; j++){
B.xx[i][j]= p * A.xx[i][j]; // Multiply each entry by p
}
}
return B;
}

// Definition of multiplication between a matrix, A and a scalar, p
matrix operator *(const matrix& A, const double& p){

// Create a matrix with the same dimensions as A
matrix B(A.rows,A.columns);
for (int i=0; i<A.rows; i++){
for (int j=0; j<A.columns; j++){
B.xx[i][j]= p * A.xx[i][j]; // Multiply each entry by p
}
}
return B;
}


// Definition of division of a matrix, A by a scalar, p i.e. A/p
matrix operator /(const matrix& A, const double& p){

// Create a matrix with the same dimensions as A
matrix B(A.rows,A.columns);

  for (int i=0; i<A.rows; i++){
    for (int j=0; j<A.columns; j++){
            B.xx[i][j]= A.xx[i][j]/p; // Divide each entry by p
            }
        }
    return B;
    }

// Define multiplication for matrices:
matrix operator *(const matrix& A, const matrix& B){

// Use assertion to check matrix dimensions are consistent
assert(A.columns==B.rows);

    // Create a result matrix, C with the correct dimensions
    matrix C(A.rows,B.columns);
    double temp = 0;
    // rows (m)
    for (int i=0; i < A.rows ; i++){
        // columns (q)
        for (int j=0; j < B.columns ; j++){
            // dot product step (n sums)
            for (int k=0; k < A.columns ; k++){
                temp = temp + A.xx[i][k]*B.xx[k][j];
                }
            // Set the C matrix values
            C.xx[i][j]=temp;
            // reset temp
            temp = 0;
            }
        }
    return C;

}

////////////////////// Unary operators ////////////////////////////////////
matrix operator +(const matrix& A){ // Define the unary operator, '+'{
// Create a temporary matrix with the same dimensions as A
matrix B(A.rows,A.columns);
// Set the entires of B to be the same as those in A
for (int i=0; i < A.rows; i++){
        for (int j=0; j < A.columns; j++){
        B.xx[i][j] = A.xx[i][j];
        }
        }
    return B;
}

// Define the unary operator, '???'
matrix operator - (const matrix& A)  // define the unary operator '-'
  {
// Create a temporary matrix with the same dimensions as A
matrix B(A.rows,A.columns);
// Set the entires of B to be the same as those in A
for (int i=0; i < A.rows; i++){
        for (int j=0; j < A.columns; j++){
            B.xx[i][j] = -A.xx[i][j];
      }
        }
    return B;
    }

// Overload the trans operator to mean transpose
matrix t(const matrix& A){
    // Create a temporary matrix with reversed dimensions
    matrix B(A.columns,A.rows);
    // Set the entires of B to be the same as those in A
    for (int i=0; i < A.columns; i++){
        for (int j=0; j < A.rows; j++)
        {
        B.xx[i][j] = A.xx[j][i];
          }
        }
    return B;
    }

/////////////////////////// Other operators //////////////////////////////
// Definition of matrix operator '='
// Operator returns a matrix equal to the RHS
matrix& matrix::operator =(const matrix &A){
// Destruct previous entries
for (int i=0; i < rows; i++)
      {
        delete[] xx[i];
        }
    delete[] xx;
    // Assign new dimensions to be equal to that of the RHS
     rows = A.rows;
     columns = A.columns;
    
    // Allocate the memory as in the constructor
    xx = new double *[A.rows];
    
    for (int i=0; i < A.rows ; i++){
        xx[i] = new double[A.columns];
        }
    // Copy the values across from the RHS
     for (int i=0; i < A.rows; i++){
      for (int j=0; j < A.columns ; j++){
          // Set entries to be the same as the RHS matrix
          xx[i][j]=A.xx[i][j];
          }
      }
    return *this;
    }

// Allows reference to the entries of a matrix in the same way as MATLAB
// Can call or assign values.
double &matrix::operator() (int i, int j){
  if (i < 1 || j < 1){
        std::cout << "Error: One of your indices may have been too small \n\n";
        } else if (i > rows || j > columns){
        std::cout << "Error: One of your indices may have been too large \n\n";
        }
    return xx[i - 1][j - 1];
    }

//////////////////////// Function prototypes /////////////////////////////
int rows(matrix A);
int columns(matrix A);

//////////////////////// Function definitons /////////////////////////////
// Returns the private field, 'rows'
int rows(matrix A){
    return A.rows;
  }

// Returns the private field, 'columns'
int columns(matrix A){
    return A.columns;
}

///////////////////////// Functions that are friends ////////////////////
// Overloads the '<<' operator to allow easy printing of matrices
ostream& operator << (ostream& output, const matrix &A){
    for (int i=0; i<A.rows; i++){
        for (int j=0; j < A.columns; j++)
        {
        output << " " << A.xx[i][j];
        }
        output << "\n";
        }
    output << "\n";
    return output;
    }


matrix eye(int size){
// Create a temporary matrix with the same dimensions as A
matrix temp_eye(size, size);
// Set the entries of B to be the same as those in A
for (int i=0; i < size; i++){
temp_eye.xx[i][i] = 1;
}
return temp_eye;
}


// Function that returns an nxn permutation matrix which swaps rows i and j
matrix permute_r(int n, int i, int j){
// Create nxn identity matrix
matrix I = eye(n);
// Zero the diagonal entries in the given rows
I(i,i)=0;
I(j,j)=0;

// Set the appropriate values to be 1
I(i,j)=1;
I(j,i)=1;
return I;
}

// Function that returns the row number of the largest
// sub???diagonal value of a given column
int find_pivot(matrix A, int column){
// Initialise maxval to be diagonal entry in column, 'column'
double maxval = fabs(A(column,column));
// Initialise rowval to be column
int rowval=column;
for (int i=column+1; i <= A.rows; i++)
      {
      if ( fabs(A(i,column)) > maxval)
      {
      // Update maxval and rowval if bigger than previous maxval
      maxval = fabs(A(i,column));
      rowval = i;
      }
      }
    return rowval;
    }

// Function that returns an mxn matrix with entries that
// are the same as matrix A, where possible
matrix resize(matrix A, int m, int n){
    int p,q;
    matrix Mout(m,n);
    // select lowest of each matrix dimension
    if (m <= A.rows){
        p=m;
    }else{
        p=A.rows;
        }
    if (n <= A.columns)
      {
      q=n;
      }else{
        q=A.columns;
        }
    // copy across relevant values
    for (int i=1; i <= p; i++){
        for (int j=1; j <= q; j++){
            Mout(i,j) = A(i,j);
            }
        }
    return Mout;
}
