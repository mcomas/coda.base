#ifndef coda_base_H
#define coda_base_H

#include <RcppArmadillo.h>

arma::mat clr_coordinates(arma::mat &X);

arma::mat inv_clr_coordinates(arma::mat clrX);

arma::mat ilr_basis_default(unsigned int dim);

arma::mat ilr_basis_simplex(unsigned int dim);

arma::mat ilr_coordinates(arma::mat &X);

arma::mat ilr_coordinates_with_basis(arma::mat X, arma::mat B);

arma::mat inv_ilr_coordinates(arma::mat ilrX);

arma::mat ilr_to_alr(unsigned int dim);

#endif
