//#include <Rcpp.h>
//using namespace Rcpp;

#include "matrix.h"
#include "vector.h"

// Definition of division of a vector, b by a matrix, A i.e. y=b/A
matrix operator /(const matrix& b, const matrix& A)
{
  int n = A.rows;
  
  // Create empty matrices, P & L
  matrix P, L;
  
  // Create and intialise U & Atemp
  matrix Utemp = eye(n);
  matrix Atemp=A;
  
  //std::cout << U << "\n\n";
  
  for (int j=1; j < n ; j++)
  {
    //std::cout << "Need to permute row " << j << " with row ";
    //std::cout << find pivot(Atemp,j) << "\n\n";
    
    // Create appropriate permutation matrix, P
    P = permute_r(n,find_pivot(Atemp,j),j);
    
    Utemp = P*Utemp; // Update U & Atemp
    Atemp = Utemp*A;
    //std::cout << "Permute rows \n\n" << Atemp;
    L = eye(n);
    for (int i=j+1; i <= n ; i++)
    {
      // Check for division by zero
      assert(fabs(Atemp(j,j))>1.0e-15);
      
      // Compute multiplier and store in sub???diagonal entry of L
      L(i,j)= -Atemp(i,j)/Atemp(j,j);
    }
    
    Utemp = L*Utemp;
    
    Atemp = Utemp*A;
    
    //std::cout << "Eliminate sub???diagonal entries \n\n" << Atemp;
    
  }
  
  
  // Now loop through and set to zero any values which are almost zero
  
  for (int j=1; j < n ; j++)
  {
    for (int i=j+1; i <= n ; i++)
    {
      if (fabs(Atemp(i,j)) < 5.0e-16)
      {
        Atemp(i,j)=0;
      }
    }
  }
  
  
  //std::cout << "The matrix U = Utemp*A is then: \n\n" << Atemp;
  
  // So, to solve Ax=b, we do: (Utemp*A)x=Utemp*b i.e.
  // Set U=Utemp*A=Atemp, compute y=Utemp*b and
  // solve Ux=y (upper triangular system ???> back subs)
  
  matrix U = Utemp*A;   //Atemp; gives the same result
  matrix y = Utemp*b;
  
  //std::cout << "The RHS is then Utemp*b: \n\n" << y ;
  
  matrix x(n,1); // Create result vector
  
  // Solve Ux=y by back substitution:
  
  // Compute last entry of vector x (first step in back subs)
  x(n,1)=y(n,1)/U(n,n);
  
  double temp = 0; // Initialise temp
  
  for (int i=n - 1; i >= 1; i--)
  {
    temp = y(i,1);
    for (int j = n; j > i; j--)
    {
      temp = temp - U(i,j)*x(j,1);
    }
    x(i,1)=temp/U(i,i);
  }
  return x;
}
