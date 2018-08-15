#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"
#include "sim_measures.h"

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


void fill_triples(std::vector<Eigen::Triplet<double> >& tl, std::vector<double>& res, 
                  const std::vector<std::tuple<double,double,int> >& index1, const std::vector<std::tuple<double,double,int> >& index2,
                  int offset, int i, bool use_threshold, double min_value, int top_n){
  // fill the triples based on the results per row (i), using one of two options: top_n filtering or regular
  // the use_pair makes it somewhat complicated, but is necessary because res is a vector with only the results for pairs used for the current row. 
  // the positions in use_pair are the actual positions, and the true values match the positions in res.
  if (top_n > 0 && top_n < res.size()) {
    std::vector<std::pair<double,int> > res_index = index_and_sort_top_n<double>(res, top_n, offset);
    for (int res_i = 0; res_i < top_n; res_i++) {
      if (use_threshold && res_index[res_i].first < min_value) continue;
      if (res_index[res_i].first == 0) continue;
      tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[i]), std::get<2>(index2[res_index[res_i].second]), res_index[res_i].first));
    }
  } else {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      if (use_threshold && res[res_i] < min_value) continue;
      if (res[res_i] == 0) continue;
      tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[i]), std::get<2>(index2[res_i+offset]), res[res_i]));
    }
  }
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> batched_tcrossprod_cpp(Eigen::SparseMatrix<double>& m1, Eigen::SparseMatrix<double>& m2,
                                                   IntegerVector group1, IntegerVector group2, 
                                                   NumericVector order1, NumericVector order2,
                                                   bool use_threshold=true, double min_value=0, int top_n=0, bool diag=true, bool only_upper=false, bool rowsum_div=false, bool l2norm=false, std::string crossfun="prod",
                                                   int lwindow=0, int rwindow=0,
                                                   bool verbose=false, int batchsize = 10000) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  std::vector<std::tuple<double,double,int> > index1, index2;
  index1 = create_index(group1,order1);
  index2 = create_index(group2,order2);
  
  // temp hack to prevent useless batches. Once everything is tested we should make a separate tcrossprod_cpp function without group and order
  bool not_batched = unique(group1).size() == 1 && unique(group2).size() == 1 && unique(order1).size() == 1 && unique(order2).size() == 1;
  if (not_batched) batchsize = m1.rows();
    
  m1 = sm_prepare(m1, index1, true, l2norm);   // sort and transpose m1
  m2 = sm_prepare(m2, index2, false, l2norm);  // only sort m2
  
  int rows = m1.cols();
  int cols = m2.rows();
  if (rows != cols && only_upper) stop("using 'only_upper = true' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  if (rows != cols && !diag) stop("using 'diag = false' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  
  if (crossfun != "prod" && crossfun != "min") stop("Not a valid crossfun (currently supports prod and min)");
  
  
  std::vector<Eigen::Triplet<double> > tl;
  if (top_n > 0) {         // for top_n we know max required capacity, otherwise guess (and update in loop)
    tl.reserve(rows*top_n);
  } else {
    tl.reserve(rows*10);  
  }
  
  Progress p(rows, verbose);   
  
  Eigen::SparseMatrix<double> m2_batch;
  
  std::vector<bool> use_pair;
  std::vector<int> positions;
  int i_end;
  int offset = 0;
  
  std::pair<int,int> i_indices;
  std::pair<int,int> batch_indices;
  
  for (int i = 0; i < rows; i++) {
    // CREATE BATCH
    if (i % batchsize == 0) {
      i_end = i + (batchsize - 1);
      if (i_end > rows) i_end = rows - 1;
      batch_indices = find_positions(index2, std::get<0>(index1[i]),
                                     std::get<0>(index1[i_end]),
                                     std::get<1>(index1[i]) + lwindow,
                                     std::get<1>(index1[i_end]) + rwindow);
      m2_batch = m2.middleRows(batch_indices.first, batch_indices.second - batch_indices.first);
      offset = batch_indices.first;
    } 
    std::vector<double> res(m2_batch.rows());
    
    
    // FIND ALL MATCHES FOR i IN BATCH
    std::vector<bool> use_pair(m2_batch.rows());
    double group1_val = std::get<0>(index1[i]);
    double order1_val = std::get<1>(index1[i]);
    double group2_val, order2_val;
    for (int j = 0; j < use_pair.size(); j++) {
      group2_val = std::get<0>(index2[j + offset]);
      if (group2_val != group1_val) continue;
      order2_val = std::get<1>(index2[j + offset]);
      if (order2_val < order1_val+lwindow) continue;
      if (order2_val > order1_val+rwindow) continue;
      if (!diag) if (i == j + offset) continue;
      if (only_upper) if (i > j + offset) continue;
      use_pair[j] = true;
    }
    
    // MANAGE TRIPLET VECTOR CAPACITY
    if (tl.capacity() < tl.size() + res.size()) {
      double pct_to_go = 1 - i/double(rows);
      tl.reserve(tl.capacity() * (1.2 + 0.8*pct_to_go));
    }
    
    // CROSSPROD (updates res by reference)
    if (crossfun == "prod" && !rowsum_div) sim_product(i, m1, m2_batch, res, use_pair);
    if (crossfun == "prod" && rowsum_div) sim_product_pct(i, m1, m2_batch, res, use_pair);
    if (crossfun == "min" && !rowsum_div) sim_min(i, m1, m2_batch, res, use_pair);
    if (crossfun == "min" && rowsum_div) sim_min_pct(i, m1, m2_batch, res, use_pair);

    // SAVE VALUES WITH CORRECT POSITIONS
    fill_triples(tl, res, index1, index2, offset, i, use_threshold, min_value, top_n);
    if (Progress::check_abort())
      stop("Aborted");
    p.increment(1);
  }
  
  Eigen::SparseMatrix<double> out(rows,cols);
  out.setFromTriplets(tl.begin(), tl.end());
  return out;
}








