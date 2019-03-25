#ifndef __SIM_MEASURES_H_INCLUDED__
#define __SIM_MEASURES_H_INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

void as_pct(int, const Eigen::SparseMatrix<double>&, std::vector<double>& res);

void sim_product(int, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&,
                 std::vector<double>&, std::vector<bool>&);

void sim_min(int, const Eigen::SparseMatrix<double>& m1, const Eigen::SparseMatrix<double>& m2,
                     std::vector<double>&, std::vector<bool>&);

void sim_softprod(int, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, 
                 std::vector<double>&, std::vector<bool>&, const Eigen::SparseMatrix<double>&);


#endif