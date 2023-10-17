//#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "htool/htool.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace htool;
using Rcpp::as;



//==============================================================================
//                     Funciones creadas en RcppArmadillo
//==============================================================================



////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f1(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.75 * exp(- (pow((9*x - 2), 2) + pow((9*y - 2), 2)) / 4);
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f2(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.75 * exp(- (pow((9*x + 1), 2) / 49 + pow((9*y + 1), 2) / 10));
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f3(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.5 * exp(- (pow((9*x - 7), 2) + pow((9*y - 3), 2)) / 4);
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f4(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.2 * exp(- (pow((9*x - 4), 2) + pow((9*y - 7), 2)));
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat testfunction(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = f1(x, y) + f2(x, y) + f3(x, y) - f4(x, y);
  
  return value;
}
//______________________________________________________________________________________



////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat testfunction_random(arma::mat x, arma::mat y){
  int n  = x.n_rows;
  arma::mat value(size(x));

  value = f1(x, y) + f2(x, y) + f3(x, y) - f4(x, y) + arma::randn(n);

  return value;
}
//______________________________________________________________________________________



////////////////////////////////////////////////////////////
// [[Rcpp::export]]
void meshgrid(arma::mat & x, arma::mat & y, arma::vec & xv, arma::vec &yv){
  
  /* Copia los valores de xv en cada fila de x, también copia los valores de 
   yv en cada columna de y. 
   x es una matriz de tamaño (y.n_elem, x.n_elem)
   y es una matriz de tamaño (y_n_elem, x.n_elem)
   */
  
  int i;
  
  const int xv_l = xv.n_elem;
  const int yv_l = yv.n_elem;
  
  /* Copia de los valores de x */
  
  for(i=0; i<xv_l; i++){
    x.col(i).fill(as_scalar(xv(i)));
  }
  
  /* Copia de los valores de y */
  for(i=0; i<yv_l; i++){
    y.row(i).fill(as_scalar(yv(i)));
  }
  
  return;
}
//_____________________________________________________________________________











////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat DistanceMatrix(const arma::mat dsites, const arma::mat ctrs){
  
  int M  = dsites.n_rows;
  int dim = dsites.n_cols;
  int N   = ctrs.n_rows;
  
  arma::mat DM_data(M, N, arma::fill::zeros);
  
  for (int i = 0; i < M; i++){
    for (int j = 0; j < N; j++){
      DM_data(i, j) = norm(dsites.row(i) - ctrs.row(j));
    }
  }
  return(DM_data);
}
//_______________________________________________________________________________



////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat euclidean_dist(const arma::mat dsites){
  
  int n  = dsites.n_rows;
  // int dim = dsites.n_cols;
  // int N   = ctrs.n_rows;
  arma::mat DM_data(n, n, fill::zeros);
  //arma::mat DM_data(M, N, arma::fill::zeros);
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < i; j++){
      DM_data(i, j) = sqrt(sum(pow(dsites.row(i) - dsites.row(j), 2)));
      DM_data(j, i) = DM_data(i, j);
    }
  }
  return(DM_data);
}
//_______________________________________________________________________________



////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat radialFunction(arma::mat & r, const int RBFtype, const double R, const int shape){
  
  r = r / R;
  
  arma::mat phi(r.n_rows, r.n_cols);
  phi.zeros();
  
  if(RBFtype == 1)               // linear spline (R1)
    phi = r;
  
  if(RBFtype == 2){               // Thin Plate Spline (TPS)
    arma::uvec rGt0 = arma::find(r > 0);
    //phi.elem(rGt0) = (1/(8*M_PI))*(arma::pow(r.elem(rGt0),2) % log(r.elem(rGt0)));
    phi.elem(rGt0) = arma::pow(shape*r.elem(rGt0),2) % log(shape*r.elem(rGt0));
  }
  
  if(RBFtype == 3)               // Gaussian (GS)
    phi = exp(-arma::pow(shape*r,2));
  
  if(RBFtype == 4)               // Multiquadric (MQ)
    phi = sqrt(1 + arma::pow(shape*r,2));
  
  if(RBFtype == 5)               // Inverse multiquadric (IMQ)
    phi = 1 / sqrt(1 + arma::pow(shape*r,2));
  
  // if(RBFtype == 6 & r < 1)      // Compact support
  //   phi = arma::pow((1 - r), 2);
  //
  
  return(phi);
}
//_____________________________________________________________________________







/////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List PLS(arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval,
               const int alpha, const int shape){
  
  int N   = dsites.n_rows;
  int dim = dsites.n_cols;
  int M   = ctrs.n_rows;
  int ncol = 1 + dim;
  
  arma::mat DM = DistanceMatrix(dsites, ctrs);
  
  // CM de linea 14 de programa 19.2
  arma::mat IM0 = radialFunction(DM, RBFtype, R, shape);          // Matriz E del libro Green and Silverman
  
  // Matriz de TPS penalizada por alpha*I 
  arma::mat IM = IM0 + eye(IM0.n_rows, IM0.n_cols)*alpha; 
  
  
  // Polinomial precision (1, s1, s2)
  arma::mat PM  = join_rows(ones(N,1), dsites);                  // termino polinomial (Matriz T) 
  arma::mat PtM = trans(PM);          // matriz transpuesta termino polinomial (T transpuesta)
  
  // CM de linea 16 programa 19.2
  arma::mat aux  = join_rows(IM,PM);
  arma::mat aux2 = join_rows(PtM, zeros(ncol, ncol));
  
  // Matriz grande A
  arma::mat A = join_cols(aux, aux2);  
  
  
  // Simulacion de una variable respuesta en base a una funcion declarada previamente
  // linea 17 de programa 19.2
  arma::colvec rhs = testfunction(dsites.col(0), dsites.col(1));
  
  // rhs2 + un ruido Normal (0,1)
  arma::colvec rhs2 = rhs + 0.03*randn(rhs.n_rows, rhs.n_cols);
  
  // linea 19 del programa 19.2 
  // Vector de variable respuesta con fila 1 = rhs3 y fila 2 matriz 3x3 de ceros
  arma::colvec rhs3 = join_cols(rhs2,zeros(3, rhs.n_cols));
  
  
  arma::vec grid = linspace(0, 1, neval);
  
  arma::mat xs(neval, neval);
  arma::mat ys(neval, neval);
  arma::mat epoints(neval*neval, neval*neval);
  
  meshgrid(xs, ys, grid, grid);
  
  epoints = join_rows(vectorise(xs), vectorise(ys));
  
  
  // Calcular matriz de distancia entre los puntos de evaluacion y los centros
  // linea 22 del programa 19.2
  arma::mat DM_eval = DistanceMatrix(epoints, ctrs);
  
  // Matriz de evaluacion para datos de evaluacion
  // EM de linea 23 programa 19.2
  //arma::mat CMe = tps(ep, DM_eval); 
  
  arma::mat IMe0 = radialFunction(DM_eval, RBFtype, R, shape); 
  
  // PM de linea 24 del programa 19.2
  arma::mat PMe = join_rows(ones(neval*neval, 1), epoints);
  
  // EM de linea 24 lado derecho programa 19.2 (aparece como EM = [EM PM])
  arma::mat IMe = join_rows(IMe0, PMe);
  
  arma::colvec Pf  = IMe * solve(A, rhs3);
  arma::colvec Pf0 = solve(A, rhs3);
  
  // arma::colvec yint = IMe * solve(IM, rhs3);
  // arma::colvec yridge = IMe * solve(IM, rhs3);
  // arma::colvec ysmooth = IM * solve(IM, rhs3);
  
  // Calcular solucíon exacta, por ejemplo: evaluar 'testfunction' en los puntos de evaluacion
  // linea 26 programa 19.2
  arma::colvec exact = testfunction(epoints.col(0), epoints.col(1));
  //arma::colvec exact = testfunction_random(epoints.col(0), epoints.col(1));
  
  // Calcular el error maximo en la grilla
  // linea 27 del programa 19.2
  double maxerr = norm(Pf - exact, "inf");
  
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("response") = rhs3,
                            Rcpp::Named("Epa") = IM,
                            Rcpp::Named("Pf")       = Pf,
                            Rcpp::Named("Pf0")       = Pf0,
                            Rcpp::Named("maxerr")  = maxerr,
                            Rcpp::Named("epoints")  = epoints,
                            Rcpp::Named("DM")  = DM,
                            Rcpp::Named("DM_eval")  = DM_eval,
                            Rcpp::Named("PMe")  = PMe,
                            Rcpp::Named("IMe")  = IMe);
  
}






//////////////////////////////////////////////////////////////
//                          'Htool' 
//////////////////////////////////////////////////////////////

// la matriz E+\alpha de [Green/Silverman]
class TPS_Epa:
  public IMatrix<double>{
  const std::vector<R3>& p1;
  const std::vector<R3>& p2;
  
public:
  // Constructor
  TPS_Epa(const std::vector<R3>& p10, const std::vector<R3>& p20 ):
  IMatrix(p10.size(),p20.size()),p1(p10),p2(p20), lambda(0.0){}
  
  TPS_Epa(const std::vector<R3>& p10, const std::vector<R3>& p20, double a ):
  IMatrix(p10.size(),p20.size()),p1(p10),p2(p20), lambda(a){}
  
  // Virtual function to overload
  double get_coef(const int& k, const int& j)const {
    if ( norm2(p1[j]-p2[k]) < 1e-12)
      return lambda;
    else
      return norm2(p1[j]-p2[k])*norm2(p1[j]-p2[k])*log(norm2(p1[j]-p2[k]));
  }
  
  // Matrix vector product
  std::vector<double> operator*(std::vector<double> a){
    std::vector<double> result(p1.size(),0);
    for (int j=0;j<p1.size();j++){
      for (int k=0;k<p2.size();k++){
        result[j]+=this->get_coef(j,k)*a[k];
      }
    }
    return result;
  }
  
  // Frobenius norm
  double norm(){
    double norm = 0;
    for (int j=0;j<p1.size();j++){
      for (int k=0;k<p2.size();k++){
        norm+=this->get_coef(j,k);
      }
    }
    return norm;
  }
  
private: 
  double lambda;
};




class TPS_Mat:
  public IMatrix<double>{
  
  const std::vector<R3>& p1;
  TPS_Epa Epa;
  HMatrix<double,partialACA,GeometricClustering> H_Epa;
  
public: 
  TPS_Mat(const std::vector<R3>& p1_, double a):
  IMatrix(p1_.size(),p1_.size()), p1(p1_), H_Epa(Epa, p1_), Epa(p1_, p1_, a){}
  
  // Virtual function to overload
  // NO USAR, NO ES CORRECTO!!! Solo para no obtener error al compilar
  double get_coef(const int& k, const int& j)const {
    if ( ( j < p1.size() ) && ( k < p1.size() ) )
    {
      if ( j == k)
        return 0.0;
      else
        return norm2(p1[j]-p1[k])*norm2(p1[j]-p1[k])*log(norm2(p1[j]-p1[k]));
    }
    {
      return 0;
      // implementar la evaluaci?n de T y su transpuesta
    }
  }
  
  //Matrix vector product
  std::vector<double> operator*(std::vector<double> x){
    
    unsigned int n = x.size();
    
    std::vector<double> y;
    y.reserve(n);
    
    y.resize(n-3,0);
    
    std::vector<double>::const_iterator first = x.begin();
    std::vector<double>::const_iterator last = x.end()-3;
    std::vector<double> x1(first, last);
    
    y = H_Epa*x1;
    
    y.resize(n,0.0);
    
    // el producto con transpuesta(T)
    for ( unsigned int j=0; j < n-3; j++ )
    {
      y[j] = y[j] + x[n-3] + x[n-2]*p1[j][0] + x[n-1]*p1[j][1];
    }
    
    // el producto con T
    for ( unsigned int j=0; j < n-3; j++ )
    {
      y[n-3] = y[n-3] + x[j];
      y[n-2] = y[n-2] + x[j]*p1[j][0];
      y[n-1] = y[n-1] + x[j]*p1[j][1];
    }
    
    return y;
  }
};






/////////////////////////////////////////////////////////
// [[Rcpp::export]]
void  testTps(double Epsilon,
              double Eta,
              int MinClusterSize,
              double lambda) {
  
  clock_t start, lap, end;
  std::ofstream size; size.open("size.dat"); size.close();
  std::ofstream timeh; timeh.open("timeh.dat"); timeh.close();
  std::ofstream timef; timef.open("timef.dat"); timef.close();
  
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  
  // Htool parameters
  SetEpsilon(Epsilon);
  SetEta(Eta);
  SetMinClusterSize(MinClusterSize);
  

  //Funciona para rcpp y rcpparmadillo
  int K=5;
  for ( unsigned int k=5; k<K; ++k)
  {
    unsigned int N = pow(2,k);
    std::cout << " N = " << N << std::endl;
    size.open("size.dat", std::ios::out | std::ios::app);
    size << N << "\n";
    size.close();
    std::vector<R3> p(N);     // p es del tama?o de los indices previamente declarados con I
    for(int j=0; j < N; j++){
      p[j][0] = j;
      p[j][1] = j;
      p[j][2] = 1.;
    }
    
    start = clock();
    TPS_Mat A(p,lambda);

    std::vector<double> x(N+3,1), result(N+3,0);
    //start = clock();
    result = A*x;
    end = clock();
    std::cout << "Time for assembling and matrix-vector product (hmatrix): "
              << (end-start)*1000.0 / CLOCKS_PER_SEC << std::endl;
    timeh.open("timeh.dat", std::ios::out | std::ios::app);
    timeh << (end-start)*1000.0 / CLOCKS_PER_SEC << "\n";
    timeh.close();
    
    TPS_Epa Af(p,p, lambda);
    std::vector<double> xf(N,1),resultf(N,0);
    start = clock();
    resultf = Af*xf;
    end = clock();
    std::cout << "Time for matrix-vector product (full matrix, only E): "
              << (end-start)*1000.0 / CLOCKS_PER_SEC << std::endl;
    timef.open("timef.dat", std::ios::out | std::ios::app);
    timef << (end-start)*1000.0 / CLOCKS_PER_SEC << "\n";
    timef.close();
    }
}

int main(void)
{
  testTps(0.002,10,10,1);
  return 1.;
}




// ////////////////////////////////////////////////////////////
// vec GMRES(TPS_Mat& A, std::vector<double>& b, std::vector<double>& x0, double tol){
//   
//   //double tol=1e???6
//   unsigned int n=x0.size(); // size of linear system
//   
//   // determine initial residual, r0 in vec form
//   std::vector<double> aux(n,0), aux2(n,0);
//   aux = A*x0;
//   
//   std::cout << aux.size() << std::endl;
//   
//   
//   vec r0(n);
//   for ( unsigned int j=0; j<n; ++j )
//     r0(j+1) = b[j] - aux[j];
//   
//   //std::cout << "initial residual vec, r0 = b???A*x0 : \n\n" << r0;
//   
//   // need this in least square part later
//   double normr0 = norm(r0);
//   
//   // initialise to enter while loop
//   double residual = 1.0;
//   
//   // intialise vec v
//   vec v = r0 / normr0;
//   
//   //std::cout << "initial vec v = r0 / | | ro | | : \n\n" << v;
//   
//   // Arnoldi / GMRES step index
//   int k = 1;
//   
//   // Declare Givens rotation matrix, initialise Jtotal;
//   matrix J, Jtotal;
//   Jtotal = eye(2);
//   
//   // intialise H, declare tempMat, V, w
//   matrix H(1,1), Htemp, HH, bb(1,1), c, cc;
//   matrix tempMat, V, Vold, hNewCol;
//   vec w(n), vj(v.rows);
//   
//   bb(1,1) = normr0;
//   
//   // initialise matrix V (matrix of orthogonal basis vecs)
//   V = v;
// 
//   while ( (residual > tol) && ( k < 100) )
//   {
//     //std::cout<< " \n\n";
//     
//     // update Vold (used for checking Arnoldi later)
//     Vold = V;
//     
//     H = resize(H, k+1,k);
//     
//     
//     // Arnoldi steps (using Gram???Schmidt process)
//     // w = mat2vec(A*v);
//     // we need to multiply with A using std::vector<double>
//     for ( unsigned int j=0; j<n; ++j )
//       aux[j] = v(j+1);
//     
//     aux2 = A*aux;
//     
//     for ( unsigned int j=0; j<n; ++j )
//       w(j+1) = aux2[j];
//     
//     
//     //std::cout<< "(k = " << k <<") : vec w=Av : \n\n" << w;
//     
//     for (int j=1; j <= k; j++)
//     {
//       for (int i=1; i <= rows(V); i++)
//       {// set the vec vj to be jth column of V
//         vj(i)=V(i,j);
//       }
//       
//       tempMat = t(vj)*w;
//       
//       // these two lines calculate the inner product
//       H(j,k)= tempMat(1,1);
//       //std::cout<< "H("<<j<<","<<k<<")= "<<H(j,k)<<"\n\n";
//       
//       w = w - H(j,k)*vj;
//       //std::cout<< "Gramm-Schmidt update of vec w: \n\n" << w;
//     }
//     
//     H(k +1, k) = norm(w);
//     //std::cout<< "H(" << k+1 << "," << k << ")= " << H(k+1,k) << "\n\n";
//     
//     v = w / H(k +1, k);
//     //std::cout<< "(k = " << k <<") :new vec v: \n\n" << v;
//     
//     // add one more column to matrix V
//     V = resize(V, rows(V),k+1);
//     
//     for (int i=1; i <= rows(V); i++)
//     {
//       // copy entries of v to new column of V
//       V(i, k + 1) = v(i);
//     }
//     
//     //std::cout<< "(k = " << k << ") :latest matrix, V: \n\n" << V;
//     
//     //std::cout << "(k = " << k <<") :latest matrix, H: \n\n" << H;
//     
//     //std::cout << "check: AV[k] = V[k+1]H: \n\n" << A*Vold << V*H;
//     
//     /////////////////////////////// Least squares step //////////////////////
//     
//     if (k==1)
//     {
//       // First pass through, Htemp=H
//       Htemp = H;
//     }
//     else
//     {
//       // for subsequent passes, Htemp=Jtotal*H
//       Jtotal = resize(Jtotal, k +1, k + 1);
//       Jtotal(k + 1, k + 1) = 1;
//       Htemp = Jtotal*H;
//     }
//     
//     // Form next Givens rotation matrix
//     J = eye( k - 1);
//     J = resize(J,k+1,k+1);
//     
//     J(k,k)=Htemp(k,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
//     J(k,k+1)=Htemp(k+1,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
//     J(k+1,k)= -Htemp(k+1,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
//     J(k+1,k+1)=Htemp(k,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
//     
//     //std::cout<< "J: \n\n" << J;
//     
//     // combine together with previous Givens rotations
//     Jtotal = J*Jtotal;
//     
//     //std::cout<< "Check orthogonality of Jtotal \n\n" << Â¬Jtotal*Jtotal;
//     
//     HH = Jtotal*H;
//     
//     for (int i=1; i <= k+1; i++)
//     {
//       for (int j=1; j <= k; j++)
//       {
//         // set all 'small' values to zero
//         if (abs(HH(i,j)) < 1e-15)
//         {
//           HH(i,j)=0;
//         }
//       }
//     }
//     
//     //std::cout<< "Check Jtotal*H is upper triangular: \n\n" << HH;
//     
//     bb = resize(bb,k+1,1);
//     
//     //std::cout<< "bb: \n\n" << bb;
//     
//     c=Jtotal*bb;
//     
//     //std::cout<< "c=J*bb: \n\n" << c;
//     
//     residual = fabs(c(k+1,1));
//     
//     std::cout<< k << "th residual: \n\n" << residual << "\n\n";
//     k++;
//   }
//    
//   
//   std::cout<< "GMRES iteration converged in " << k - 1 << " steps\n\n";
//   // Extract upper triangular square matrix
//   HH = resize(HH, rows(HH) - 1, columns(HH));
// 
//   //std::cout<< "HH: \n\n" << HH;
//   
//   cc = resize(c, rows(HH),1);
//   
//   //std::cout<< "cc: \n\n" << cc;
//   
//   
//   //yy = cc/HH; // solve linear system
//   matrix yy = cc / HH; // solve linear system
//   
//   vec y = mat2vec(yy);
//   
//   //std::cout<< "y: \n\n" << y;
//   
//   // chop the newest column off of matrix V
//   V = resize(V, rows(V), columns(V) - 1);
//   
//   // reuse r0
//   for ( unsigned int j=0; j<n; ++j )
//     r0(j+1) = x0[j];
//   
//   vec x = mat2vec(r0 + V*y);
//   
//   return x;
//   //return Rcpp::List::create(Rcpp::Named("x_coefficients") = x);
// }


//////////////////////////////////////////////////
//              Conjugate gradient
//////////////////////////////////////////////////
arma::vec cg_arma_full(arma::mat A, arma::vec b, float tol = 1e-14, int maxIter = 1000) {
  /* 
   Input:
   A: matrix (Aqu_ A deber_a ser la matriz H)
   b: vector
   Output
   x: vector
   */
  // get number of columns of A
  int C = A.n_cols ;
  int R = A.n_rows ;
  // initiate solution x as zeros
  arma::vec x(C) ;
  x.zeros() ; 
  // arma::vec oneVec(C);
  // oneVec.ones() ;
  
  arma::vec r = b - A * x;
  arma::vec p = r;
  double rs_old = sum( r % r );
  double rs_new=1;
  
  
  arma::vec Ap(R);
  double alpha, beta;
  // vector version of alpha
  // arma::vec alphaVec(1);
  
  unsigned int iter=0;
  for(iter = 0; (iter < maxIter) && (rs_new > tol); iter++){
    Ap = A * p;
    
    alpha = rs_old / sum(p % Ap);
    // alphaVec.fill(alpha); 
    x += alpha * p;
    r -= alpha * Ap;
    rs_new = sum(r % r);
    beta = rs_new / rs_old; 
    
    p = r + beta * p;
    rs_old = rs_new;
    if (iter >= maxIter){
      Rcout << "cg did not converge." << endl;
    }
  }
  
  std::cout << " CG iterations:" << iter << std::endl;
  return x;
  
} 

//=======================================================================================
//                          Conjugated Gradient Hmat
//=======================================================================================
// solve (E - S*M1)x = b with CG
arma::vec cg_arma_hmat(HMatrix<double,partialACA,GeometricClustering> &E, arma::mat &M1, arma::mat &S, arma::vec &b, float tol = 1e-14, int maxIter = 1000) {
  
  // get size of problem
  int N = b.n_rows;
  
  // initialize solution x as zeros
  arma::vec x(N);
  x.zeros();
  arma::vec xaux1(N), xaux2(N);
  xaux1.zeros(N);
  xaux2.zeros(N);
  
  std::vector<double> xh(N,0.0);
  std::vector<double> xhaux(N,0.0);
  
  //------- r = b - A*x
  xaux1 = M1*x;
  xaux2 = S*xaux1;
  arma::vec r = b + xaux2;
  
  xhaux = E*xh;
  
  for (unsigned int j=0; j<N; ++j )
    r[j] = r[j] - xhaux[j];
  //------- r = b - A*x
  
  arma::vec p = r;
  double rs_old = sum( r % r );
  double rs_new=1;
  
  
  arma::vec Ap(N);
  double alpha, beta;
  // vector version of alpha
  // arma::vec alphaVec(1);
  unsigned int iter=0;
  for(iter = 0; (iter < maxIter) && (rs_new > tol); iter++){
    
    //------- Ap = A*p
    xaux1 = M1*p;
    xaux2 = S*xaux1;
    
    for (unsigned int j=0; j<N; ++j )
      xh[j] = p[j];
    xhaux = E*xh;
    
    for (unsigned int j=0; j<N; ++j )
      Ap[j] = xhaux[j] - xaux2[j];
    //------- Ap = A*p
    
    alpha = rs_old / sum(p % Ap);
    // alphaVec.fill(alpha); 
    x += alpha * p;
    r -= alpha * Ap;
    rs_new = sum(r % r);
    beta = rs_new / rs_old; 
    
    p = r + beta * p;
    rs_old = rs_new;
    if (iter >= maxIter){
      Rcout << "cg did not converge." << endl;
    }
  }
  
  std::cout << " CG iterations:" << iter << std::endl;
  return x;
}

// [[Rcpp::export]]
void mpiinit()
{
  MPI_Init(NULL,NULL);
}

// [[Rcpp::export]]
void mpifinalize()
{
  MPI_Finalize();
}


// solve the Schur complement version with H-matrices
// [[Rcpp::export]]
arma::vec calculateTPS_hmat(arma::mat &sites, arma::vec &values, arma::mat &M1, arma::mat &S,
                            double Epsilon, double Eta, int MinClusterSize, double lambda){
  
  unsigned int N = sites.n_rows;        
  
  // Htool parameters
  SetEpsilon(Epsilon);
  SetEta(Eta);
  SetMinClusterSize(MinClusterSize);
  
  std::vector<R3> sitios1(N);
  for(unsigned int j=0; j < N; j++){
    sitios1[j][0] = sites(j,0);
    sitios1[j][1] = sites(j,1);
    sitios1[j][2] = 1.;
  }
  
  // levantar la matriz H
  TPS_Epa Epa(sitios1,sitios1,lambda);
  HMatrix<double,partialACA,GeometricClustering> H_Epa(Epa,sitios1);
  
  arma::vec result = cg_arma_hmat(H_Epa,M1,S,values);
  //arma::vec result(values.n_rows);
  //result.zeros();
  
  //MPI_Finalize();
  return result; 
}

// solve the full Schur complement version with CG
// [[Rcpp::export]]
arma::vec calculateTPS_full( arma::vec &values, arma::mat &full ){
  
  //arma::vec result = cg_arma(H_Epa,M1,S,values);
  arma::vec result = cg_arma_full(full,values);
  
  return result; 
}