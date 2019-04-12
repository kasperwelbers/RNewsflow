#include <RcppEigen.h>
#include <Rcpp.h>
#include "RNewsflow_types.h"

#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

std::string prep_term(std::string term, bool parentheses) {
  if (parentheses) term = "(" + term + ")";
  return term;
}

// [[Rcpp::export]]
List term_union_cpp(SpMat& m, SpMat& simmat, std::vector<std::string> terms, std::vector<bool> parentheses, bool verbose = true, std::string sep = " | ") {
  // simmat is assumed to be symmetrical
  if (simmat.cols() != m.cols()) stop("number of columns in m and simmat not identical");
  if (simmat.rows() != simmat.cols()) stop("simmat not symmetrical");
  if (terms.size() != m.cols()) stop("number of terms not identical to number of columns in m");

  //Rcout << "1" << std::endl;
  Triplet tl;
  tl.reserve(m.nonZeros());
  
  int cluster = 0;
  std::vector<int> clusters(m.cols());
  std::vector<bool> in_cluster(m.cols());
  
  int nclust = 0;
  std::vector<std::string> clustnames;
  clustnames.reserve(m.cols());
  std::string term;
  
  Progress p(m.cols(), verbose);   
  
  bool new_cluster;
  for (int r = 0; r < m.cols(); r++) {
    new_cluster = true;
    for (SpMat::InnerIterator check_clusters(simmat, r); check_clusters; ++check_clusters) {
      if (in_cluster[check_clusters.row()]) {
        cluster = clusters[check_clusters.row()];
        new_cluster = false;
      }
    }
    if (new_cluster) {
      cluster = nclust;
      nclust++;
      clustnames.push_back("");
    }
    
    if (clustnames[cluster] == "")
      clustnames[cluster] = prep_term(terms[r], parentheses[r]);
    else
      clustnames[cluster] = clustnames[cluster] + sep + prep_term(terms[r], parentheses[r]);
    
    for (SpMat::InnerIterator termit(simmat, r); termit; ++termit) {
      if (termit.row() < r) continue;
      in_cluster[termit.row()] = true;
      clusters[termit.row()] = cluster;
      
      for (SpMat::InnerIterator it(m, termit.row()); it; ++it) {
        tl.push_back(Eigen::Triplet<double>(it.row(), cluster, 1));      
      }
    }  
    p.increment(1);
    if (tl.capacity() < 0.5*tl.size()) tl.reserve(tl.capacity() * 2);
  }

  SpMat newm(m.rows(), nclust);
  newm.setFromTriplets(tl.begin(), tl.end());
  
  List out;
  out["m"] = newm;
  out["colnames"] = clustnames;
  return out;
}

// [[Rcpp::export]]
List term_intersect_cpp(SpMat& m, SpMat& simmat, std::vector<std::string> terms, std::vector<bool> parentheses, bool verbose = true, std::string sep = " & ") {
  // simmat is assumed to be symmetrical
  if (simmat.cols() != m.cols()) stop("number of columns in m and simmat not identical");
  if (simmat.rows() != simmat.cols()) stop("simmat not symmetrical");
  if (terms.size() != m.cols()) stop("number of terms not identical to number of columns in m");
  
  
  Triplet tl;
  tl.reserve(m.nonZeros());
  
  int col = 0;
  std::vector<std::string> colnames;
  colnames.reserve(simmat.nonZeros());
  std::string term;
  
  Progress p(m.cols(), verbose);   
  int j;
  for (int i = 0; i < m.cols(); i++) {
    for (SpMat::InnerIterator it(simmat, i); it; ++it) {
      j = it.row();
      if (j < i) continue;
      
      if (i == j) 
        colnames.push_back(prep_term(terms[i], parentheses[i]));
      else
        colnames.push_back(prep_term(terms[i], parentheses[i]) + sep + prep_term(terms[j], parentheses[j]));
        
      SpMat::InnerIterator itx(m, i);
      SpMat::InnerIterator ity(m, j);
      while (itx && ity) {
        if (itx.row() < ity.row()) {
          ++itx;
          continue;
        }
        if (itx.row() > ity.row()) {
          ++ity;
          continue;
        }
        tl.push_back(Eigen::Triplet<double>(itx.row(), col, 1));      
        ++itx;
        ++ity;
      }
      col++;
    }
    p.increment(1);
    if (tl.capacity() < 0.5*tl.size()) tl.reserve(tl.capacity() * 2);
  }

  SpMat newm(m.rows(), col);
  newm.setFromTriplets(tl.begin(), tl.end());
  
  List out;
  out["m"] = newm;
  out["colnames"] = colnames;
  return out;
}


/*** R
m = rsparsematrix(10,5,0)
m[1,1] = 1; m[1,2] = 1; m[3,3] = 1; m[3,5] = 1
m[7,1] = 1; m[7,2] = 1; m[6,4] = 1; m[6,5] = 1
m

colnames(m) = c('a','b','c','d','e')

simmat = tcrossprod_sparse(t(m), min_value=0.5, normalize='l2')
term_union_cpp(m,simmat,colnames(m),parentheses=T)
term_intersect_cpp(m,simmat,colnames(m), parentheses=T)
*/


    

