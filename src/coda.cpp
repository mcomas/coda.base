#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "coda.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat pinv(arma::mat X){
  return arma::pinv(X);
}

// [[Rcpp::export]]
arma::mat c_variation_array(arma::mat X, bool include_means = false){
  unsigned int K = X.n_cols;
  arma::mat lX = log(X);
  arma::mat varray = arma::mat(K,K);
  varray.diag().zeros();
  arma::mat Xcov = cov(lX);
  if(include_means){
    arma::mat Xmeans = arma::mean(lX, 0);
    for(unsigned i = 0; i < K; i++){
      for(unsigned j = 0; j < i; j++){
        varray(i,j) = Xmeans(j) - Xmeans(i);
        varray(j,i) = Xcov(i,i) + Xcov(j,j) - 2*Xcov(i,j);
      }
    }
  }else{
    for(unsigned i = 0; i < K; i++){
      for(unsigned j = 0; j < i; j++){
        varray(i,j) = varray(j,i) = Xcov(i,i) + Xcov(j,j) - 2*Xcov(i,j);
      }
    }
  }
  return varray;
}

// [[Rcpp::export]]
arma::mat alr_basis_default(unsigned int dim){
  arma::mat B = arma::mat(dim, dim-1);
  for(unsigned int i = 0; i < dim - 1; i++){
    for(unsigned int j = 0; j < dim - 1; j++){
      B(i,j) = i == j ? 1 : 0;
    }
  }
  for(unsigned int j = 0; j < dim - 1; j++){
    B(dim-1,j) = -1;
  }
  return(B);
}

// [[Rcpp::export]]
arma::mat clr_basis_default(unsigned int dim){
  arma::mat B = arma::mat(dim, dim);
  double dim_inv = -1/(double)dim;
  for(unsigned int i = 0; i < dim; i++){
    for(unsigned int j = 0; j < dim; j++){
      B(i,j) = i == j ? 1 + dim_inv : dim_inv;
    }
  }
  return(B);
}

// [[Rcpp::export]]
arma::mat ilr_basis_default(unsigned int dim){
  arma::mat B = arma::mat(dim, dim-1);

  for(unsigned int i = 0; i < dim - 1; i++){
    unsigned int I1 = i + 1;
    unsigned int I2 = i + 2;
    double l = 1/std::sqrt((double)(I1*I2));
    double r = - std::sqrt((double)I1/I2);
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
    double l = exp(1/std::sqrt((double)(I1*I2)));
    double r = 1/exp(std::sqrt((double)I1/I2));
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
arma::mat clr_coordinates(arma::mat &X){
  arma::mat LOGX = log(X);
  arma::mat m = mean(LOGX, 1);

  LOGX.each_col() -= m;
   // for(unsigned int i = 0; i < X.n_rows; i++){
   //   LOGX.row(i) = LOGX.row(i) - mean(LOGX.row(i));
   // }

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
arma::mat alr_coordinates(arma::mat &X, int denominator){
  int idenom = X.n_cols-1;
  arma::mat res(X.n_rows, idenom);
  arma::mat logX = log(X);
  for(unsigned int j = 0; j < res.n_cols; j++){
    for(unsigned int i = 0; i< res.n_rows; i++){
      res(i,j) = logX(i,j) - logX(i,idenom);
    }
  }
  //arma::mat res = logX(arma::span::all, arma::span(0,X.n_cols-2));
  //arma::vec vdenom = logX.col(X.n_cols-1);
  //res.each_col() -= vdenom;
  return(res);
}

// [[Rcpp::export]]
arma::mat matrix_coordinates(arma::mat X, arma::mat B){
    return(log(X) * B);
}

// [[Rcpp::export]]
arma::mat sparse_coordinates(arma::mat X, arma::sp_mat B){
  return(log(X) * B);
}

// [[Rcpp::export]]
arma::mat coordinates_basis(arma::mat X, arma::mat B, bool sparse = false){
  if(sparse){
    arma::sp_mat sB(B);
    arma::mat logX = log(X);
    return(logX * sB);
  }else{
    arma::mat logX = log(X);
    return(logX * B);
  }
}

// [[Rcpp::export]]
arma::mat ilr_coordinates(arma::mat &X){
  arma::mat logX = log(X);
  arma::mat B = ilr_basis_default(X.n_cols);
  return(logX * B);
}

// [[Rcpp::export]]
arma::mat inv_ilr_coordinates(arma::mat ilrX){
  arma::mat B = ilr_basis_default(ilrX.n_cols + 1);

  return inv_clr_coordinates(ilrX * B.t());
}

