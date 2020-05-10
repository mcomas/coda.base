#include <RcppArmadillo.h>

void f(int mu, int nu, int sigma, arma::uvec& I, arma::uvec& A);
void b(int mu, int nu, int sigma, arma::uvec& I, arma::uvec& A);

void f(int mu, int nu, int sigma, arma::uvec& I, arma::uvec& A){
  if(mu == 2){
    Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
  }else{
    f(mu-1, nu-1, (mu+sigma) % 2, I, A);
  }
  if(nu == mu + 1){
    // Rcpp::Rcout << "a" << std::endl;
    A[mu-1] = mu - 1;
    Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
    while(A[nu-1] > 0){
      A[nu-1] = A[nu-1] - 1;
      Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
    }
    // Rcpp::Rcout << "b" << std::endl;
  }else if(nu > mu + 1){ // nu > mu + 1
    if( (mu + sigma) % 2 ){
      A[nu-2] = mu - 1;
    }else{
      A[mu-1] = mu - 1;
    }
    if( (A[nu-1] + sigma) % 2 ){
      b(mu, nu - 1, 0, I, A);
    }else{
      f(mu, nu - 1, 0, I, A);
    }
    while( A[nu-1] > 0){
      A[nu-1] = A[nu-1] - 1;
      if( (A[nu-1] + sigma) % 2 ){
        b(mu, nu - 1, 0, I, A);
      }else{
        f(mu, nu - 1, 0, I, A);
      }
    }
  }
}

void b(int mu, int nu, int sigma, arma::uvec& I, arma::uvec& A){
  if( nu == mu + 1){
    while(A[nu-1] < mu - 1){
      Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
      A[nu-1] = A[nu-1] + 1;
    }
    Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
    A[mu-1] = 0;
    // Rcpp::Rcout << "d" << std::endl;
  }else if(nu > mu + 1){ // nu > mu + 1:
    if( (A[nu-1] + sigma) % 2 ){
      f(mu, nu - 1, 0, I, A);
    }else{
      b(mu, nu - 1, 0, I, A);
    }
    while(A[nu-1] < mu - 1){
      A[nu-1] = A[nu-1] + 1;
      if( (A[nu-1] + sigma) % 2 ){
        f(mu, nu - 1, 0, I, A);
      }else{
        b(mu, nu - 1, 0, I, A);
      }
    }
    if( (mu + sigma) % 2 ){
      A[nu-2] = 0;
    }else{
      A[mu-1] = 0;
    }
  }
  // Rcpp::Rcout << "a" << std::endl;
  if(mu == 2){
    Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
  }else{
    b(mu - 1, nu - 1, (mu + sigma) % 2, I, A);
  }
}

// [[Rcpp::export]]
void combinations2(unsigned n) {
  arma::uvec I = arma::uvec(n);
  arma::uvec A = arma::uvec(n);

  A.fill(0);
  A[n-2] = 1;
  A[n-1] = 2;
  f(3, n, 0, I, A);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
combinations(5)
*/
