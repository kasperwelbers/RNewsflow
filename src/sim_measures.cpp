#include <Rcpp.h>
#include <RcppEigen.h>
#include "sim_measures.h"
#include "RNewsflow_types.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

void as_pct(int i, const SpMat& m, std::vector<double>& res){
  double sum_value=0;
  for (SpMat::InnerIterator it(m,i); it; ++it) {
    sum_value += it.value();
  }
  if (sum_value > 0) {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      res[res_i] = res[res_i] / sum_value;
    }  
  }
}


void pnorm_filter(std::vector<double>& res, bool log_trans=false, bool nz = false, double max_p=1) {
  NumericVector v = Rcpp::wrap(res);
  if (nz) v = v[v>0];
  if (log_trans) v = log(v+1);
  double M = mean(v);
  double SD = sd(v);
  double p;
  
  for (int i=0; i < res.size(); i++) {
    if (SD == 0 || res[i] < M) {
      res[i] = 0;
    } else {
      p = R::pnorm(res[i], M, SD, 0, 0);
      if (p > max_p) res[i] = 0;
    }
  }
}

void pbeta_filter(std::vector<double>& res, bool nz = false, double max_p=1) {
  NumericVector v = Rcpp::wrap(res);
  if (nz) v = v[v>0];
  double M = mean(v);
  double V = var(v);
  double p;

  double alpha = ((1 - M) / V - 1 / M) * pow(M, 2);
  double beta  = alpha * (1 / M - 1);
  
  for (int i=0; i < res.size(); i++) {
    p = R::pbeta(res[i], alpha, beta, 0, 0);
    if (p > max_p) res[i] = 0;
  }
}

void pdisparity_filter(std::vector<double>& res, double k, double max_p=1) {
  NumericVector v = Rcpp::wrap(res);
  double vsum = sum(v);
  double p;

  for (int i=0; i < res.size(); i++) {
    if (k == 0) {
      res[i] = 0;
    } else {
      p = res[i] / vsum;
      p = pow(1 - p, k - 1);
      if (p > max_p) res[i] = 0;
    }
  }
}


void sim_product(int i, const SpMat& m1, const SpMat& m2, 
             std::vector<double>& res, std::vector<bool>& use_pair){
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    for (SpMat::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += it1.value() * it2.value();
    }
  }
} 

void sim_maxproduct(int i, const SpMat& m1, const SpMat& m2,
             std::vector<double>& res, std::vector<bool>& use_pair){
  // special case: return highest product (for use in query_lookup)
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    for (SpMat::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += std::max(it1.value() * it2.value(), res[it2.row()]);
    }
  }
} 

void sim_min(int i, const SpMat& m1, const SpMat& m2,
                 std::vector<double>& res, std::vector<bool>& use_pair){
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    for (SpMat::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += std::min(it1.value(), it2.value());
    }
  }
} 

void sim_softprod(int i, const SpMat& m1, const SpMat& m2, 
                 std::vector<double>& res, std::vector<bool>& use_pair, const SpMat& simmat){
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    for (SpMat::InnerIterator simmat_it(simmat,it1.row()); simmat_it; ++simmat_it) {
      for (SpMat::InnerIterator it2(m2,simmat_it.row()); it2; ++it2) {
        if (!use_pair[it2.row()]) continue;
        res[it2.row()] += it1.value() * it2.value() * simmat_it.value();
      }
    }
  }
} 

void sim_lookup(int i, const SpMat& m1, const SpMat& m2, 
                 std::vector<double>& res, std::vector<bool>& use_pair){
  // for query lookup in which the value in it2 does not matter
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    for (SpMat::InnerIterator it2(m2,it1.row()); it2; ++it2) {
      if (!use_pair[it2.row()]) continue;
      res[it2.row()] += it1.value();
    }
  }
} 
