#include <RcppArmadillo.h>


void f(int mu, int nu, int sigma, arma::uvec A){
  if(mu == 2){
    Rcpp::Rcout << A.t() << std::endl;
  }else{
    f(mu-1, nu-1, mu+sigma % 2, A);
  }
  if(nu == mu + 1){

  }
}

void b(int mu, int nu, int sigma){

}

// [[Rcpp::export]]
arma::uvec combinations(unsigned n) {
  unsigned m = 3;
  arma::uvec A = arma::uvec(n);
  A.fill(0);
  for(int i = 0; i < m; i++){
    A(n-m+i) = i;
  }
  f(m, n, 0, A);
  return(A);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
a = t(combinations(5))
*/
