// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include "coda.h"

using namespace Rcpp;

//' @export
 // [[Rcpp::export]]
arma::cube c_conditional_obasis(arma::mat& C){

  int n = C.n_cols;
  int D = C.n_rows;
  int d = D - 1;

  //Rcpp::Rcout << "D: " << D << " n: " << n << std::endl;
  arma::vec d0 = arma::vec(n);
  arma::umat I0 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    d0(k) = arma::accu(C.col(k) == 0);
    I0.col(k).head(d0(k)) = arma::find(C.col(k) == 0);
  }

  arma::vec d1 = d-d0;
  arma::vec D1 = 1+d1;
  arma::umat I1 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    I1.col(k).head(D1(k)) = arma::find(C.col(k) > 0);
  }

  arma::cube B = arma::zeros(d,D,n);
  for(int k=0; k<n; k++){
    if(D1(k) == 0){
      B.slice(k) = ilr_basis_default(D).t();
    }else{
      if(d0(k) > 0){
        if(d0(k) > 1){
          arma::uvec i0 = I0.col(k).head(d0(k));
          for(unsigned int i = 0; i < d0(k) - 1; i++){
            unsigned int I1 = i + 1;
            unsigned int I2 = i + 2;
            double l = 1/std::sqrt((double)(I1*I2));
            double r = - std::sqrt((double)I1/I2);
            for(unsigned int j = 0; j < I1; j++){
              B(i,i0(j),k) = l;
            }
            B(i,i0(I1),k) = r;
          }
        }
        for(unsigned int j = 0; j < D; j++){
          if(C(j,k) == 0){
            B(d0(k) - 1, j, k) = +1/d0(k) * std::sqrt(d0(k) * (D1(k))/D);
          }else{
            B(d0(k) - 1, j, k) = -1/(D1(k)) * std::sqrt(d0(k) * (D1(k))/D);
          }
        }
        if(D1(k) > 1){
          arma::uvec i1 = I1.col(k).head(D1(k));
          for(unsigned int i = 0; i < D1(k) - 1; i++){
            unsigned int I1 = i + 1;
            unsigned int I2 = i + 2;
            double l = 1/std::sqrt((double)(I1*I2));
            double r = - std::sqrt((double)I1/I2);
            for(unsigned int j = 0; j < I1; j++){
              B(d0(k)+i,i1(j),k) = l;
            }
            B(d0(k)+i,i1(I1),k) = r;
          }
        }
      }else{
        B.slice(k) = ilr_basis_default(D1(k)).t();
      }
    }

  }
  return(B);
}

