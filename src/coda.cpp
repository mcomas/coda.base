// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "coda.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


// [[Rcpp::export]]
arma::mat ilr_basis_default(unsigned int dim){
  arma::mat B = arma::mat(dim, dim-1);

  for(unsigned int i = 0; i < dim - 1; i++){
    unsigned int I1 = i + 1;
    unsigned int I2 = i + 2;
    double l = 1/sqrt(I1*I2);
    double r = - sqrt((double)I1/I2);
    for(unsigned int j = 0; j < I1; j++){
      B(j,i) = l;
    }
    B(I1,i) = r;
    for(unsigned int j = I2; j < dim; j++){
      B(j,i) = 0;
    }
  }
  return(B);
}

// [[Rcpp::export]]
arma::mat ilr_basis_simplex(unsigned int dim){
  arma::mat B = arma::mat(dim, dim-1);

  for(unsigned int i = 0; i < dim - 1; i++){
    unsigned int I1 = i + 1;
    unsigned int I2 = i + 2;
    double l = exp(1/sqrt(I1*I2));
    double r = 1/exp(sqrt((double)I1/I2));
    for(unsigned int j = 0; j < I1; j++){
      B(j,i) = l;
    }
    B(I1,i) = r;
    for(unsigned int j = I2; j < dim; j++){
      B(j,i) = 1;
    }
    B.col(i) = B.col(i) / (I1 * l + r + dim - I2);
  }
  return(B);
}


// [[Rcpp::export]]
arma::mat ilr_to_alr(unsigned int dim){
  unsigned int k = dim-1;
  arma::mat B = ilr_basis_default(dim);
  arma::mat ILR_TO_ALR = B(arma::span(0,k-1),arma::span(0,k-1));
  for(unsigned int i = 0;i < k;i++){
    ILR_TO_ALR.col(i) = ILR_TO_ALR.col(i) - B(k,i);
  }
  return(ILR_TO_ALR);
}

// [[Rcpp::export]]
arma::mat clr_coordinates(arma::mat X){
  arma::mat LOGX = log(X);

  for(unsigned int i = 0; i < X.n_rows; i++){
    LOGX.row(i) = LOGX.row(i) - mean(LOGX.row(i));
  }

  return(LOGX);
}

// [[Rcpp::export]]
arma::mat inv_clr_coordinates(arma::mat clrX){
  arma::mat X = arma::exp(clrX);
  arma::mat sumX = arma::sum(X, 1);
  for(unsigned int i = 0; i < X.n_cols; i++){
    X.col(i) = X.col(i) / sumX;
  }
  return(X);
}

// [[Rcpp::export]]
arma::mat coordinates(arma::mat X, arma::mat B){
  arma::mat logX = log(X);
    return(logX * B);
}

// [[Rcpp::export]]
arma::mat coordinates_default(arma::mat X){
  arma::mat logX = log(X);
  arma::mat B = ilr_basis_default(X.n_cols);
  return(logX * B);
}

// [[Rcpp::export]]
arma::mat inv_ilr_coordinates(arma::mat ilrX){
  arma::mat B = ilr_basis_default(ilrX.n_cols + 1);

  return inv_clr_coordinates(ilrX * B.t());
}

