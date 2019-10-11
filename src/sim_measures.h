#ifndef __SIM_MEASURES_H_INCLUDED__
#define __SIM_MEASURES_H_INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

void as_pct(int, const Eigen::SparseMatrix<double>&, std::vector<double>& res);

void pnorm_filter(std::vector<double>& res, bool log_trans, bool nz, double max_p);
void pbeta_filter(std::vector<double>& res, bool nz, double max_p);
void pdisparity_filter(std::vector<double>& res, double k, double max_p);

void sim_product(int, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&,
                 std::vector<double>&, std::vector<bool>&);

void sim_maxproduct(int, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&,
                 std::vector<double>&, std::vector<bool>&);

void sim_min(int, const Eigen::SparseMatrix<double>& m1, const Eigen::SparseMatrix<double>& m2,
                     std::vector<double>&, std::vector<bool>&);

void sim_softprod(int, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, 
                 std::vector<double>&, std::vector<bool>&, const Eigen::SparseMatrix<double>&);

void sim_lookup(int, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&,
                 std::vector<double>&, std::vector<bool>&);

#endif