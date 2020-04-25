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


void sim_cp_lookup(int i, const SpMat& m1, const SpMat& m2, 
                  std::vector<double>& res, std::vector<bool>& use_pair, const SpMat& simmat, bool normalize=false){
  // this measure looks at all unique pairs of terms in m1 and m2, and sums the scores of these pairs in simmat
  // the simmat is thus more accurately a matrix of term and term cooccurrence weights, but we use the term 
  // simmat simply because that's how it was called for the softcosine measure. (in time we might tidy this up) 
  
  // In this function the normalization is built-in, because that's just hell of a lot easier and also quite efficient
  double doc_sum = 0;
  
  // first find all nz terms in m1
  std::vector<bool> m1_terms(m1.rows());
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    m1_terms[it1.row()] = true;
  }
  
  // then iterate over all terms in m1
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    
    // remember if the first term of the cooc (x) occurs in m2 doc
    std::vector<bool> has_x(m2.rows());
    for (SpMat::InnerIterator it2(m2, it1.row()); it2; ++it2) {
      has_x[it2.row()] = true;
    }
    
    // then iterate over the cooc scores y for x

    // test: if a term has more than 2 cooccurence (including with itself)
    // we need to count it again to 
    double xvalue;
    int xterm_count = 0;
    std::vector<int> xterm_count_it2(m2.rows());
    
    for (SpMat::InnerIterator simmat_it(simmat,it1.row()); simmat_it; ++simmat_it) {
      if (it1.row() > simmat_it.row()) continue;
      if (it1.row() == simmat_it.row())
        xvalue = simmat_it.value();
      
      // skip if y did not occur in m1
      if (!m1_terms[simmat_it.row()]) continue;
      
      
      if (normalize) 
        doc_sum += simmat_it.value();   
      
      xterm_count += 1;
      if (xterm_count > 2) {
        doc_sum += xvalue;
      }
      
      // iterate over all y terms in m2.
      for (SpMat::InnerIterator it2(m2,simmat_it.row()); it2; ++it2) {
        
        // only write to documents that are valid pairs AND in which x occurred 
        if (!use_pair[it2.row()]) continue;
        if (!has_x[it2.row()]) continue;
        
        
        // use the exact value in simmat (so preprocessing of simmat should determine weights)
        res[it2.row()] += simmat_it.value();
        
        xterm_count_it2[it2.row()] += 1;
        if (xterm_count_it2[it2.row()] > 2)
          res[it2.row()] += xvalue;
      }
    }
  }
  if (normalize && doc_sum > 0) {
    for (int i=0; i < res.size(); i++) 
      res[i] = res[i] / doc_sum;
  }
} 

void sim_cp_lookup2(int i, const SpMat& m1, const SpMat& m2, 
                   std::vector<double>& res, std::vector<bool>& use_pair, const SpMat& simmat, bool normalize=false){
  // alternative idea.
  // make it so that the observed frequencies in m1 and m2 are used. For each combination of terms, take the min of the values
  // this way, for the same term (diagonal of simmat) the outcome is simply the observed value, making this a generalization of vector space models
  // for combinations of terms, the min observed frequency is the number of times the two terms co-occurr in the document.
  
  double doc_sum = 0;
  
  // first find all nz terms in m1
  std::vector<double> m1_terms(m1.rows());
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    m1_terms[it1.row()] = it1.value();
  }
  
  // then iterate over all terms in m1
  for (SpMat::InnerIterator it1(m1,i); it1; ++it1) {
    
    // remember if the first term of the cooc (x) occurs in m2 doc
    std::vector<double> yt1_vec(m2.rows());
    for (SpMat::InnerIterator it2(m2, it1.row()); it2; ++it2) {
      yt1_vec[it2.row()] = it2.value();
    }
    
    // then iterate over the cooc scores y for x
    
    for (SpMat::InnerIterator simmat_it(simmat,it1.row()); simmat_it; ++simmat_it) {
      
      //if (xi < yi) break;
      if (it1.row() < simmat_it.row()) break;
      // skip if y did not occur in m1
      if (m1_terms[simmat_it.row()] == 0) continue;
      
      double xt1 = it1.value();
      double xt2 = m1_terms[simmat_it.row()];
      double x = std::min(xt1, xt2);
    
      if (normalize) doc_sum += x * simmat_it.value();
      
      // iterate over all y terms in m2.
      for (SpMat::InnerIterator it2(m2,simmat_it.row()); it2; ++it2) {
        
        // only write to documents that are valid pairs AND in which x occurred 
        if (!use_pair[it2.row()]) continue;
        if (yt1_vec[it2.row()] == 0) continue;
        
        double yt1 = yt1_vec[it2.row()];
        double yt2 = it2.value();
        double y = std::min(yt1, yt2);
        
        // use the exact value in simmat (so preprocessing of simmat should determine weights)
        res[it2.row()] += x * y * simmat_it.value();
      }
    }
  }
  if (normalize && doc_sum > 0) {
    for (int i=0; i < res.size(); i++) 
      res[i] = res[i] / doc_sum;
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
