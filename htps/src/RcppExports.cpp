// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// f1
arma::mat f1(arma::mat x, arma::mat y);
RcppExport SEXP _testTPS_f1(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f1(x, y));
    return rcpp_result_gen;
END_RCPP
}
// f2
arma::mat f2(arma::mat x, arma::mat y);
RcppExport SEXP _testTPS_f2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// f3
arma::mat f3(arma::mat x, arma::mat y);
RcppExport SEXP _testTPS_f3(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f3(x, y));
    return rcpp_result_gen;
END_RCPP
}
// f4
arma::mat f4(arma::mat x, arma::mat y);
RcppExport SEXP _testTPS_f4(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f4(x, y));
    return rcpp_result_gen;
END_RCPP
}
// testfunction
arma::mat testfunction(arma::mat x, arma::mat y);
RcppExport SEXP _testTPS_testfunction(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(testfunction(x, y));
    return rcpp_result_gen;
END_RCPP
}
// testfunction_random
arma::mat testfunction_random(arma::mat x, arma::mat y);
RcppExport SEXP _testTPS_testfunction_random(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(testfunction_random(x, y));
    return rcpp_result_gen;
END_RCPP
}
// meshgrid
void meshgrid(arma::mat& x, arma::mat& y, arma::vec& xv, arma::vec& yv);
RcppExport SEXP _testTPS_meshgrid(SEXP xSEXP, SEXP ySEXP, SEXP xvSEXP, SEXP yvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type xv(xvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yv(yvSEXP);
    meshgrid(x, y, xv, yv);
    return R_NilValue;
END_RCPP
}
// DistanceMatrix
arma::mat DistanceMatrix(const arma::mat dsites, const arma::mat ctrs);
RcppExport SEXP _testTPS_DistanceMatrix(SEXP dsitesSEXP, SEXP ctrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ctrs(ctrsSEXP);
    rcpp_result_gen = Rcpp::wrap(DistanceMatrix(dsites, ctrs));
    return rcpp_result_gen;
END_RCPP
}
// euclidean_dist
arma::mat euclidean_dist(const arma::mat dsites);
RcppExport SEXP _testTPS_euclidean_dist(SEXP dsitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dsites(dsitesSEXP);
    rcpp_result_gen = Rcpp::wrap(euclidean_dist(dsites));
    return rcpp_result_gen;
END_RCPP
}
// radialFunction
arma::mat radialFunction(arma::mat& r, const int RBFtype, const double R, const int shape);
RcppExport SEXP _testTPS_radialFunction(SEXP rSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(radialFunction(r, RBFtype, R, shape));
    return rcpp_result_gen;
END_RCPP
}
// PLS
Rcpp::List PLS(arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval, const int alpha, const int shape);
RcppExport SEXP _testTPS_PLS(SEXP dsitesSEXP, SEXP ctrsSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP nevalSEXP, SEXP alphaSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ctrs(ctrsSEXP);
    Rcpp::traits::input_parameter< int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type neval(nevalSEXP);
    Rcpp::traits::input_parameter< const int >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(PLS(dsites, ctrs, RBFtype, R, neval, alpha, shape));
    return rcpp_result_gen;
END_RCPP
}
// testTps
void testTps(double Epsilon, double Eta, int MinClusterSize, double lambda);
RcppExport SEXP _testTPS_testTps(SEXP EpsilonSEXP, SEXP EtaSEXP, SEXP MinClusterSizeSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Epsilon(EpsilonSEXP);
    Rcpp::traits::input_parameter< double >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< int >::type MinClusterSize(MinClusterSizeSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    testTps(Epsilon, Eta, MinClusterSize, lambda);
    return R_NilValue;
END_RCPP
}
// mpiinit
void mpiinit();
RcppExport SEXP _testTPS_mpiinit() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    mpiinit();
    return R_NilValue;
END_RCPP
}
// mpifinalize
void mpifinalize();
RcppExport SEXP _testTPS_mpifinalize() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    mpifinalize();
    return R_NilValue;
END_RCPP
}
// calculateTPS_hmat
arma::vec calculateTPS_hmat(arma::mat& sites, arma::vec& values, arma::mat& M1, arma::mat& S, double Epsilon, double Eta, int MinClusterSize, double lambda);
RcppExport SEXP _testTPS_calculateTPS_hmat(SEXP sitesSEXP, SEXP valuesSEXP, SEXP M1SEXP, SEXP SSEXP, SEXP EpsilonSEXP, SEXP EtaSEXP, SEXP MinClusterSizeSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type sites(sitesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type Epsilon(EpsilonSEXP);
    Rcpp::traits::input_parameter< double >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< int >::type MinClusterSize(MinClusterSizeSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateTPS_hmat(sites, values, M1, S, Epsilon, Eta, MinClusterSize, lambda));
    return rcpp_result_gen;
END_RCPP
}
// calculateTPS_full
arma::vec calculateTPS_full(arma::vec& values, arma::mat& full);
RcppExport SEXP _testTPS_calculateTPS_full(SEXP valuesSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type full(fullSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateTPS_full(values, full));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_testTPS_f1", (DL_FUNC) &_testTPS_f1, 2},
    {"_testTPS_f2", (DL_FUNC) &_testTPS_f2, 2},
    {"_testTPS_f3", (DL_FUNC) &_testTPS_f3, 2},
    {"_testTPS_f4", (DL_FUNC) &_testTPS_f4, 2},
    {"_testTPS_testfunction", (DL_FUNC) &_testTPS_testfunction, 2},
    {"_testTPS_testfunction_random", (DL_FUNC) &_testTPS_testfunction_random, 2},
    {"_testTPS_meshgrid", (DL_FUNC) &_testTPS_meshgrid, 4},
    {"_testTPS_DistanceMatrix", (DL_FUNC) &_testTPS_DistanceMatrix, 2},
    {"_testTPS_euclidean_dist", (DL_FUNC) &_testTPS_euclidean_dist, 1},
    {"_testTPS_radialFunction", (DL_FUNC) &_testTPS_radialFunction, 4},
    {"_testTPS_PLS", (DL_FUNC) &_testTPS_PLS, 7},
    {"_testTPS_testTps", (DL_FUNC) &_testTPS_testTps, 4},
    {"_testTPS_mpiinit", (DL_FUNC) &_testTPS_mpiinit, 0},
    {"_testTPS_mpifinalize", (DL_FUNC) &_testTPS_mpifinalize, 0},
    {"_testTPS_calculateTPS_hmat", (DL_FUNC) &_testTPS_calculateTPS_hmat, 8},
    {"_testTPS_calculateTPS_full", (DL_FUNC) &_testTPS_calculateTPS_full, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_testTPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
