// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// Function to find inverse Cholesky factor
// Idea of the code is adopted from Golub & van Loan (1996)

// [[Rcpp::export]]
arma::mat inv_chol(arma::mat A){
  int n = A.n_rows;
  double nr = A(0,0);
  // Normalizing the matrix (if required):
  if(nr!=1){
    A = A/nr;
  }
  arma::rowvec av = A.row(0);
  arma::colvec avec = av.t();
  arma::colvec r = avec.tail(n-1);
  arma::colvec y = zeros(n);
  arma::mat R = eye(n,n);
  // Initial values:
  y(0) = - r(0);
  double a = - r(0);
  double b = 1;
  // Updating b:
  b = (1 - (a*a))*b;
  R(0,1) = y(0);
  R.col(1) = R.col(1)/pow(b,0.5);
  // Updating other columns:
  //int k = 0;
  for(int k = 0; k <= (n-3); k++){
    arma::colvec rsub = r.head(k+1);
    arma::colvec rrev = reverse(rsub);
    arma::colvec ysub = y.head(k+1);
    arma::colvec yrev = reverse(ysub);
    // Updating a:
    a = -(r(k+1) + sum(rrev.t()*ysub))/b;
    // Updating y:
    y.head(k+1) = ysub + a*yrev; y(k+1) = a;
    // Updating b:
    b = (1 - (a*a))*b;
    // Updating the columns:
    arma::colvec ysub2 = y.head(k+2);
    arma::colvec yrev2 = reverse(ysub2);
    R.submat(0,k+2,k+1,k+2) = yrev2;
    R.col(k+2) = R.col(k+2)/pow(b,0.5);
  }
  R = R/pow(nr,0.5);
  return R;
}

