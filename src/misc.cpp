#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"
#include "RNewsflow_types.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

std::vector<double> get_row_l2(SpMat& m) {
  // calculate column sums for an Eigen matrix
  std::vector<double> out(m.rows());
  
  for (int k=0; k < m.outerSize(); ++k) {
    for (SpMat::InnerIterator it(m, k); it; ++it) {
      out[it.row()] += pow(it.value(),2);
    }
  }
  for (int i=0; i < out.size(); i++) {
    out[i] = pow(out[i], 0.5);
  }
  return(out);
}



std::vector<double> softcos_row_mag(const SpMat& m, const SpMat& simmat, bool verbose) {
  SpMat m1(m.transpose());
  std::vector<double> out(m1.cols());
  std::vector<double> tmp_sim(simmat.rows());
  
  Progress p(simmat.cols(), verbose); 
  for (int i = 0; i < simmat.cols(); i++) {
    tmp_sim = std::vector<double>(simmat.rows());
    for (SpMat::InnerIterator simmat_it(simmat, i); simmat_it; ++simmat_it) {
      tmp_sim[simmat_it.row()] = simmat_it.value();
    }
    for (SpMat::InnerIterator it1(m,i); it1; ++it1) {
      for (SpMat::InnerIterator it2(m1,it1.row()); it2; ++it2) {
        out[it1.row()] += it1.value() * it2.value() * tmp_sim[it2.row()];
      }
    }
    if (Progress::check_abort())
      stop("Aborted");
    p.increment(1);
  }
  
  for (int i=0; i < out.size(); i++) {
    out[i] = pow(out[i], 0.5);
  }
  return(out);
}

SpMat batch_simmat_prepare(SpMat& m2_batch, const SpMat& simmat) {
  std::vector<bool> non_zero_term(m2_batch.cols());
  for (int nz_i = 0; nz_i < m2_batch.cols(); nz_i++) {
    for (SpMat::InnerIterator nz_it(m2_batch, nz_i); nz_it; ++nz_it) {
      non_zero_term[nz_i] = true;
      continue;
    }
  }
  
  std::vector<Eigen::Triplet<double>> tl(simmat.nonZeros());
  for (int i = 0; i < simmat.cols(); i++) {
    for (SpMat::InnerIterator simmat_it(simmat, i); simmat_it; ++simmat_it) {
      if (non_zero_term[simmat_it.row()]) {
        tl.push_back(Eigen::Triplet<double>(simmat_it.row(), simmat_it.col(), simmat_it.value()));
      }
    }
  }
  
  SpMat out(simmat.cols(), simmat.rows());
  out.setFromTriplets(tl.begin(), tl.end());
  return(out);
}


SpMat sm_prepare(SpMat& m, Index index, const SpMat& simmat, bool transpose, std::string normalize) {
  if (m.rows() != index.size()) stop("number of rows is not equal to length of index");
  
  std::vector<int> tvec(index.size());
  for (int i = 0; i < index.size(); i++) tvec[std::get<2>(index[i])] = i;
  
  std::vector<double> l2;
  bool do_normalize = normalize != "none";
  if (normalize == "l2") l2 = get_row_l2(m);
  if (normalize == "softl2") l2 = softcos_row_mag(m, simmat, false);
    
  int row, col;
  std::vector<Eigen::Triplet<double>> tl(m.nonZeros());
  for (int k=0; k < m.outerSize(); ++k) {
    for (SpMat::InnerIterator it(m, k); it; ++it) {
      if (transpose) {
        row = it.col();
        col = tvec[it.row()];
      } else {
        row = tvec[it.row()];
        col = it.col();
      }
      if (do_normalize) 
        tl.push_back(Eigen::Triplet<double>(row, col, it.value() / l2[it.row()]));
      else 
        tl.push_back(Eigen::Triplet<double>(row, col, it.value()));
    }
  }

  
  SpMat out;
  if (transpose) 
    out = SpMat(m.cols(), m.rows());
  else 
    out = SpMat(m.rows(), m.cols());
    
  out.setFromTriplets(tl.begin(), tl.end());
  return(out);
}


Index create_index(Rcpp::IntegerVector group, 
                                                         Rcpp::NumericVector order) {
  std::vector<double> g;
  std::vector<double> o;
  g = as<std::vector<double> >(group);
  o = as<std::vector<double> >(order);
  Index out = index_and_sort<double,double>(g,o); 
  return(out);
}


bool search_group_l ( const std::tuple<double,double,int>& vec, const double& val) {return (val > std::get<0>(vec));}
bool search_group_u ( const double& val, const std::tuple<double,double,int>& vec) {return (val < std::get<0>(vec));}

bool search_order_l ( const std::tuple<double,double,int>& vec, const double& val) {return (val > std::get<1>(vec));}
bool search_order_u ( const double& val, const std::tuple<double,double,int>& vec) {return (val < std::get<1>(vec));}

std::pair<int,int> find_positions(Index& xi, const double min_g, const double max_g, const double min_o, const double max_o) {
  Index::iterator first_group_low, first_group_up, first, last_group_low, last_group_up, last;

  first_group_low = std::lower_bound(xi.begin(), xi.end(), min_g, search_group_l);
  first_group_up = std::upper_bound(xi.begin(), xi.end(), min_g, search_group_u);
  first = std::lower_bound(first_group_low, first_group_up, min_o, search_order_l);
  
  last_group_low = std::lower_bound(xi.begin(), xi.end(), max_g, search_group_l);
  last_group_up = std::upper_bound(xi.begin(), xi.end(), max_g, search_group_u);
  last = std::upper_bound(last_group_low, last_group_up, max_o, search_order_u);

  Index::iterator last_old = std::lower_bound(last_group_low, last_group_up, max_o, search_order_l);
  
  int firstv = (first - xi.begin());
  int lastv = (last - xi.begin());
  // note that last isn't strictly the position, it's the position after the last position.
  
  std::pair<int,int> out(firstv, lastv);
  return(out);
}

double sum_std_vec(std::vector<double>& x) {
  double out = 0;
  for (auto& v : x)
    out += v;
  return(out);
}

double nz_std_vec(std::vector<double>& x) {
  double out = 0;
  for (auto& v : x)
    if (!(v == 0)) out += 1;
  return(out);
}

/*** R
*/
