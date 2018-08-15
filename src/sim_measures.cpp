#include <Rcpp.h>
#include <RcppEigen.h>
#include "sim_measures.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

void as_pct(double sum_value, std::vector<double>& res){
  if (sum_value > 0) {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      res[res_i] = res[res_i] / sum_value;
    }  
  }
}

void sim_product(int i, const Eigen::SparseMatrix<double>& m1, const Eigen::SparseMatrix<double>& m2, 
             std::vector<double>& res, std::vector<bool>& use_pair){
  for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
    for (Eigen::SparseMatrix<double>::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += it1.value() * it2.value();
    }
  }
} 

void sim_min(int i, const Eigen::SparseMatrix<double>& m1, const Eigen::SparseMatrix<double>& m2,
                 std::vector<double>& res, std::vector<bool>& use_pair){
  for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
    for (Eigen::SparseMatrix<double>::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += std::min(it1.value(), it2.value());
    }
  }
} 

void sim_product_pct(int i, const Eigen::SparseMatrix<double>& m1, const Eigen::SparseMatrix<double>& m2,
                 std::vector<double>& res, std::vector<bool>& use_pair){
  double sum_value=0;
  for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
    sum_value += it1.value();
    for (Eigen::SparseMatrix<double>::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += it1.value() * it2.value();
    }
  }
  as_pct(sum_value, res);
} 

void sim_min_pct(int i, const Eigen::SparseMatrix<double>& m1, const Eigen::SparseMatrix<double>& m2,
                     std::vector<double>& res, std::vector<bool>& use_pair){
  double sum_value=0;
  for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
    sum_value += it1.value();
    for (Eigen::SparseMatrix<double>::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += std::min(it1.value(), it2.value());
    }
  }
  as_pct(sum_value, res);
} 
