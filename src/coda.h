#ifndef CODA_BASE_H
#define CODA_BASE_H

#include <RcppArmadillo.h>

arma::mat pinv(const arma::mat& X);

arma::mat c_variation_array(const arma::mat& X,
                            bool include_means = false,
                            bool ml_covariance = false);

arma::mat alr_basis_default(unsigned int dim);
arma::mat clr_basis_default(unsigned int dim);
arma::mat ilr_basis_default(unsigned int dim);
arma::mat ilr_basis_simplex(unsigned int dim);

arma::mat clr_coordinates(const arma::mat& X);
arma::mat inv_clr_coordinates(const arma::mat& clrX);
arma::mat alr_coordinates(const arma::mat& X, unsigned int denominator);
arma::mat ilr_coordinates(const arma::mat& X);
arma::mat inv_ilr_coordinates(const arma::mat& ilrX);

arma::mat matrix_coordinates(const arma::mat& X, const arma::mat& B);
arma::mat sparse_coordinates(const arma::mat& X, const arma::sp_mat& B);

arma::mat ilr_to_alr(unsigned int dim);

#endif
