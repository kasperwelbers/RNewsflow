#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"
#include "sim_measures.h"
#include "RNewsflow_types.h"

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// fill the triples based on the results per row (i), using one of two options: top_n filtering or regular
// the use_pair makes it somewhat complicated, but is necessary because res is a vector with only the results for pairs used for the current row. 
// the positions in use_pair are the actual positions, and the true values match the positions in res.
void fill_triples(Triplet& tl, std::vector<double>& res, 
                  const Index& index1, const Index& index2,
                  int offset, int i, bool use_min, NumericVector min_value, bool use_max, NumericVector max_value, int top_n){
  bool min_value_vec = min_value.size() > 1;
  bool max_value_vec = max_value.size() > 1;
  int min_value_i = 0;
  int max_value_i = 0;
  
  if (top_n > 0 && top_n < res.size()) {
    std::vector<std::pair<double,int> > res_index = index_and_sort_top_n<double>(res, top_n, offset);
    for (int res_i = 0; res_i < top_n; res_i++) {
      if (min_value_vec) min_value_i = std::get<2>(index1[i]);
      if (max_value_vec) max_value_i = std::get<2>(index1[i]);
      if (use_min && res_index[res_i].first < min_value[min_value_i]) continue;
      if (use_max && res_index[res_i].first > max_value[max_value_i]) continue;
      if (res_index[res_i].first == 0) continue;
      tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[i]), std::get<2>(index2[res_index[res_i].second]), res_index[res_i].first));
    }
  } else {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      if (min_value_vec) min_value_i = std::get<2>(index1[i]);
      if (max_value_vec) max_value_i = std::get<2>(index1[i]);
      if (use_min && res[res_i] < min_value[min_value_i]) continue;
      if (use_max && res[res_i] > max_value[max_value_i]) continue;
      if (res[res_i] == 0) continue;
      tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[i]), std::get<2>(index2[res_i+offset]), res[res_i]));
    }
  }
}

// used in batched_tcrossprod_cpp
void fill_pair_information(std::vector<bool>& use_pair, std::vector<bool>& is_lag, 
                                        int i, int offset, Index& index1, Index& index2, SpMat& m2_batch, bool diag, bool only_upper, int lwindow, int rwindow) {
  use_pair = std::vector<bool>(m2_batch.rows());
  is_lag = std::vector<bool>(m2_batch.rows());
  double group1_val = std::get<0>(index1[i]);
  double order1_val = std::get<1>(index1[i]);
  double group2_val, order2_val;
  for (int j = 0; j < use_pair.size(); j++) {
    group2_val = std::get<0>(index2[j + offset]);
    order2_val = std::get<1>(index2[j + offset]);
    if (order2_val < order1_val) is_lag[j] = true;
    if (group2_val != group1_val) continue;
    if (order2_val < order1_val+lwindow) continue;
    if (order2_val > order1_val+rwindow) continue;
    if (!diag) if (i == j + offset) continue;
    if (only_upper) if (i > j + offset) continue;
    use_pair[j] = true;
  }
}

void fill_row_attributes(int i, bool row_attr, bool col_attr, bool lag_attr,
                         Index& index1, Index& index2, int offset,
                         std::vector<double> res, std::vector<bool> use_pair, std::vector<bool> is_lag, 
                         NumericVector& row_n, NumericVector& row_sum, NumericVector& row_nz,
                         NumericVector& col_n, NumericVector& col_sum, NumericVector& col_nz,
                         NumericVector& lag_n, NumericVector& lag_sum, NumericVector& lag_nz) {
  if (row_attr) {
    row_n[std::get<2>(index1[i])] = count(use_pair.begin(), use_pair.end(), true);
    row_sum[std::get<2>(index1[i])] = sum_std_vec(res);
    row_nz[std::get<2>(index1[i])] = nz_std_vec(res);
  }
  if (col_attr) {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      col_n[std::get<2>(index2[res_i+offset])] += use_pair[res_i];
      col_sum[std::get<2>(index2[res_i+offset])] += res[res_i];
      if (!res[res_i] == 0) col_nz[std::get<2>(index2[res_i+offset])] += 1;
    }
  }
  if (lag_attr) {
    for (int pair_i = 0; pair_i < use_pair.size(); pair_i++) {
      if (use_pair[pair_i]) {
        if (is_lag[pair_i]) {
          lag_n[std::get<2>(index1[i])] += 1;
          lag_sum[std::get<2>(index1[i])] += res[pair_i];
          if (!res[pair_i] == 0) lag_nz[std::get<2>(index1[i])] += 1;
        } 
      }
    }
  }
}


// [[Rcpp::export]]
List batched_tcrossprod_cpp(SpMat& m1, SpMat& m2, 
                             IntegerVector group1, IntegerVector group2, 
                             NumericVector order1, NumericVector order2,
                             const SpMat& simmat,
                             bool use_min=true, NumericVector min_value=0, bool use_max=false, NumericVector max_value=1, int top_n=0, bool diag=true, bool only_upper=false, 
                             bool rowsum_div=false, std::string pvalue="NA", double max_p=1, 
                             std::string normalize="none", std::string crossfun="prod",
                             int lwindow=0, int rwindow=0,
                             bool row_attr=false, bool col_attr=false, bool lag_attr=false,
                             bool verbose=false, int batchsize = 1000) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
  if (min_value.size() != 1) 
    if (min_value.size() != m1.rows()) stop("min_value needs to be a scalar or a vector in which each value matches a row in m1");
    
  Index index1, index2;
  index1 = create_index(group1,order1);
  index2 = create_index(group2,order2);
  
  // temp hack to prevent useless batches. Once everything is tested we should make a separate tcrossprod_cpp function without group and order
  bool not_batched = unique(group1).size() == 1 && unique(group2).size() == 1 && unique(order1).size() == 1 && unique(order2).size() == 1;
  if (not_batched) batchsize = m1.rows();
    
  m1 = sm_prepare(m1, index1, simmat, true, normalize);   // sort and transpose m1
  m2 = sm_prepare(m2, index2, simmat, false, normalize);  // only sort m2
  
  int rows = m1.cols();
  int cols = m2.rows();
  if (rows != cols && only_upper) stop("using 'only_upper = true' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  if (rows != cols && !diag) stop("using 'diag = false' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  
  if (crossfun != "prod" && crossfun != "min" && crossfun != "softprod" && crossfun != "maxproduct" && crossfun != "lookup") stop("Not a valid crossfun (currently supports prod, min, softcos and maxproduct)");
  if (normalize != "none" && normalize != "l2" && normalize != "softl2") stop("Not a valid normalize function (currently supports none, l2 and softl2)");
  
  Triplet tl;
  if (top_n > 0) {         // for top_n we know max required capacity, otherwise guess (and update in loop)
    tl.reserve(rows*top_n);
  } else {
    tl.reserve(rows*5);  
  }
  
  Progress p(rows, verbose);   
  
  SpMat m2_batch;
  SpMat batch_simmat;
  std::vector<bool> use_pair;
  std::vector<int> positions;
  int i_end;
  int offset = 0;
  std::pair<int,int> batch_indices;
  
  
  // optional row attributes
  NumericVector row_n, row_sum, row_nz, col_n, col_sum, col_nz, lag_n, lag_sum, lag_nz;
  if (row_attr) {
    row_n =   NumericVector(rows);
    row_sum = NumericVector(rows);
    row_nz =  NumericVector(rows);
  }
  if (col_attr) {
    col_n =   NumericVector(cols);
    col_sum = NumericVector(cols);
    col_nz = NumericVector(cols);
  }
  if (lag_attr) {
    lag_n =   NumericVector(rows);
    lag_sum = NumericVector(rows);
    lag_nz =  NumericVector(rows);
  }

  
  for (int i = 0; i < rows; i++) {
    
    // CREATE BATCH
    if (i % batchsize == 0) {
      i_end = i + (batchsize - 1);
      if (i_end > rows) i_end = rows - 1;
      
      batch_indices = find_positions(index2, std::get<0>(index1[i]),
                                     std::get<0>(index1[i_end]),
                                     std::get<1>(index1[i]) + lwindow,
                                     std::get<1>(index1[i_end]) + rwindow);
      
      // note that batch_indices second is the position after the last item that matches the search (like x.end())
      // therefore batch_indices.second - batch_indices.first correctly gives the number of rows after the row at batch_indices.first.
      if (batch_indices.first >= batch_indices.second) continue;
      m2_batch = m2.middleRows(batch_indices.first, batch_indices.second - batch_indices.first);
      
      if (crossfun == "softprod") {
        batch_simmat = batch_simmat_prepare(m2_batch, simmat);
        //batch_softcos_mag = batch_softcos_mag_prepare(m2_batch, simmat);
      }
      
      offset = batch_indices.first;
    } 
    if (batch_indices.first >= batch_indices.second) continue;
    
    std::vector<double> res(m2_batch.rows());

    // FIND ALL MATCHES FOR i IN BATCH (look up once, instead of for each iteration in the crossprod)
    std::vector<bool> use_pair, is_lag; 
    fill_pair_information(use_pair, is_lag, i, offset, index1, index2, m2_batch, diag, only_upper, lwindow, rwindow);
    
    // MANAGE TRIPLET VECTOR CAPACITY
    if (tl.capacity() < tl.size() + res.size()) {
      //Rcout << tl.capacity() << std::endl;
      double pct_to_go = 1 - i/double(rows);
      tl.reserve((tl.size() + res.size()) * (1.2 + 0.8*pct_to_go));
    }
    
    // CROSSPROD (fills res by reference)
    if (crossfun == "prod") sim_product(i, m1, m2_batch, res, use_pair);
    if (crossfun == "maxproduct") sim_maxproduct(i, m1, m2_batch, res, use_pair);      
    if (crossfun == "min") sim_min(i, m1, m2_batch, res, use_pair);
    if (crossfun == "softprod") sim_softprod(i, m1, m2_batch, res, use_pair, batch_simmat);      
    if (crossfun == "lookup") sim_lookup(i, m1, m2_batch, res, use_pair);
    
    if (rowsum_div) as_pct(i, m1, res); 
    
    // IF PVALUE IS USED
    if (max_p < 1) {
      if (pvalue == "normal") pnorm_filter(res, false, false, max_p);
      if (pvalue == "beta") pbeta_filter(res, false, max_p);
      if (pvalue == "lognormal") pnorm_filter(res, true, false, max_p);
      if (pvalue == "nz_normal") pnorm_filter(res, false, true, max_p);
      if (pvalue == "nz_lognormal") pnorm_filter(res, true, true, max_p);
      
      double k = count(use_pair.begin(), use_pair.end(), true);
      if (pvalue == "disparity") pdisparity_filter(res, k, max_p);
    }
    
    // SAVE COL/ROW ATTRIBUTES WITH CORRECT POSITIONS (using index)
    // (fill by reference)
    fill_row_attributes(i, row_attr, col_attr, lag_attr,
                        index1, index2, offset,
                        res, use_pair, is_lag,
                        row_n, row_sum, row_nz,
                        col_n, col_sum, col_nz,
                        lag_n, lag_sum, lag_nz);
    
      
    // SAVE VALUES WITH CORRECT POSITIONS (using offset and index)
    fill_triples(tl, res, index1, index2, offset, i, use_min, min_value, use_max, max_value, top_n);
    if (Progress::check_abort())
      stop("Aborted");
    
    p.increment(1);
  }
  
  SpMat out(rows,cols);
  out.setFromTriplets(tl.begin(), tl.end());
  
  List row_attributes;
  if (row_attr) {
    row_attributes["row_n"] = row_n;
    row_attributes["row_sum"] = row_sum;
    row_attributes["row_nz"] = row_nz;
  }
  if (col_attr) {
    row_attributes["col_n"] = col_n;
    row_attributes["col_sum"] = col_sum;
    row_attributes["col_nz"] = col_nz;
  }
  if (lag_attr) {
    row_attributes["lag_n"] = lag_n;
    row_attributes["lag_sum"] = lag_sum;
    row_attributes["lag_nz"] = lag_nz;
  }
    
  return List::create(Named("cp") = out, Named("margin") = row_attributes);
}


