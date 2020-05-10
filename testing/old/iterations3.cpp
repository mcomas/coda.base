#include <RcppArmadillo.h>

void f(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int *p);
void b(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int *p);

inline void visit(arma::uvec& I, arma::uvec& A,
                  std::vector<arma::uvec>& P, int *p){
  // Rcpp::Rcout << "A: " << A.tail(A.n_elem-1).t();
  Rcpp::Rcout << "L:";
  for(int i = 0; i<p[1];i++) Rcpp::Rcout << " " << P[1][i];
  Rcpp::Rcout << "   R:";
  for(int i = 0; i<p[2];i++) Rcpp::Rcout << " " << P[2][i];
  Rcpp::Rcout << std::endl;
}
inline void update(unsigned e, unsigned s1, unsigned s2,
                   arma::uvec& I, arma::uvec& A,
                   std::vector<arma::uvec>& P, int *p){
  // Reduce origin set size
  p[s1] = p[s1] - 1;
  // Remove e from s1 add position I[e]
  P[s1][I[e]] = P[s1][p[s1]];
  // Index of element P[s1][I[e]] to I[e]
  I[P[s1][I[e]]] = I[e];

  // Add element to s2
  P[s2][p[s2]] = e;
  // Index of element e to p[s2]
  I[e] = p[s2];
  // Increment destination set size
  p[s2] = p[s2] + 1;

}
void f(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int *p){
  unsigned e, s1, s2;
  if(mu == 2){
    visit(I, A, P, p);
  }else{
    f(mu-1, nu-1, (mu+sigma) % 2, I, A, P, p);
  }
  if(nu == mu){
    e = mu-1;
    s1 = 0;
    s2 = mu-1;
    A[e] = s2;
    update(e - 1, s1, s2, I, A, P, p);
    visit(I, A, P, p);
    while(A[nu] > 0){
      e = nu;
      s1 = A[nu];
      s2 = A[nu]-1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
      visit(I, A, P, p);
    }
  }else if(nu > mu){ // nu > mu + 1
    if( (mu + sigma) % 2 ){
      e = nu-1;
      s1 = 0;
      s2 = mu-1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    }else{
      e = mu-1;
      s1 = 0;
      s2 = mu-1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    }
    if( (A[nu] + sigma) % 2 ){
      b(mu, nu - 1, 0, I, A, P, p);
    }else{
      f(mu, nu - 1, 0, I, A, P, p);
    }
    while( A[nu] > 0){
      e = nu;
      s1 = A[nu];
      s2 = A[nu] - 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);

      if( (A[nu] + sigma) % 2 ){
        b(mu, nu - 1, 0, I, A, P, p);
      }else{
        f(mu, nu - 1, 0, I, A, P, p);
      }
    }
  }
}

void b(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int *p){
  unsigned e, s1, s2;
  if( nu == mu){
    while(A[nu] < mu - 1){
      visit(I, A, P, p);
      e = nu;
      s1 = A[nu];
      s2 = A[nu] + 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    }
    visit(I, A, P, p);
    e = mu-1;
    s1 = A[mu-1];
    s2 = 0;
    A[e] = s2;
    update(e - 1, s1, s2, I, A, P, p);
  }else if(nu > mu){ // nu > mu + 1:
    if( (A[nu] + sigma) % 2 ){
      f(mu, nu - 1, 0, I, A, P, p);
    }else{
      b(mu, nu - 1, 0, I, A, P, p);
    }
    while(A[nu] < mu - 1){
      e = nu;
      s1 = A[nu];
      s2 = A[nu] + 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
      if( (A[nu] + sigma) % 2 ){
        f(mu, nu - 1, 0, I, A, P, p);
      }else{
        b(mu, nu - 1, 0, I, A, P, p);
      }
    }
    if( (mu + sigma) % 2 ){
      e = nu-1;
      s1 = A[e];
      s2 = 0;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    }else{
      e = mu-1;
      s1 = A[e];
      s2 = 0;
      A[e] = 0;
      update(e - 1, s1, s2, I, A, P, p);
    }
  }
  if(mu == 2){
    visit(I, A, P, p);
  }else{
    b(mu - 1, nu - 1, (mu + sigma) % 2, I, A, P, p);
  }
}

// [[Rcpp::export]]
void combinations3(unsigned n) {
  unsigned N = n - 1;
  std::vector<arma::uvec> P = std::vector<arma::uvec>(3);
  P[0] = arma::uvec(N);
  P[1] = arma::uvec(N);
  P[2] = arma::uvec(N);
  int p[3];
  p[0] = N-2;
  p[1] = 1;
  p[2] = 1;
  for(int i = 0;i < N-2; i++) P[0][i] = i;
  P[1][0] = N-2;
  P[2][0] = N-1;
  arma::uvec I = arma::uvec(N);
  for(int i = 0;i < N-2; i++) I(i) = i;
  I(N-2) = 0;
  I(N-1) = 0;

  arma::uvec A = arma::uvec(n);
  A.fill(0);
  A[n-2] = 1;
  A[n-1] = 2;
  f(3, n-1, 0, I, A, P, p);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
combinations3(5)
*/
