#include <Rcpp.h>
#include <RcppEigen.h>
#include "misc.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

std::vector<double> get_row_l2(Eigen::SparseMatrix<double>& m) {
  // calculate column sums for an Eigen matrix
  std::vector<double> out(m.rows());
  
  for (int k=0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      out[it.row()] += pow(it.value(),2);
    }
  }
  for (int i=0; i < out.size(); i++) {
    out[i] = pow(out[i], 0.5);
  }
  return(out);
}

std::vector<double> softcos_row_mag(const Eigen::SparseMatrix<double>& m1, NumericMatrix& simmat) {
  std::vector<double> out(m1.cols());
  
  for (int k = 0; k < m1.cols(); k++){
    for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,k); it1; ++it1) {
      for (Eigen::SparseMatrix<double>::InnerIterator it2(m1,k); it2; ++it2) {
        if (simmat(it1.row(), it2.row()) == 0) continue; 
        out[k] += it1.value() * it2.value() * simmat(it1.row(), it2.row());
      }
    }
  }

  for (int i=0; i < out.size(); i++) {
    out[i] = pow(out[i], 0.5);
  }
  return(out);
}


Eigen::SparseMatrix<double> sm_prepare(Eigen::SparseMatrix<double>& m, std::vector<std::tuple<double,double,int> > index, bool transpose, bool l2norm) {
  if (m.rows() != index.size()) stop("number of rows is not equal to length of index");
  
  std::vector<int> tvec(index.size());
  for (int i = 0; i < index.size(); i++) tvec[std::get<2>(index[i])] = i;
  
  std::vector<double> l2;
  if (l2norm) l2 = get_row_l2(m);
  
  int row, col;
  std::vector<Eigen::Triplet<double>> tl(m.nonZeros());
  for (int k=0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      if (transpose) {
        row = it.col();
        col = tvec[it.row()];
      } else {
        row = tvec[it.row()];
        col = it.col();
      }
      if (l2norm) 
        tl.push_back(Eigen::Triplet<double>(row, col, it.value() / l2[it.row()]));
      else 
        tl.push_back(Eigen::Triplet<double>(row, col, it.value()));
    }
  }

  
  Eigen::SparseMatrix<double> out;
  if (transpose) 
    out = Eigen::SparseMatrix<double>(m.cols(), m.rows());
  else 
    out = Eigen::SparseMatrix<double>(m.rows(), m.cols());
    
  out.setFromTriplets(tl.begin(), tl.end());
  return(out);
}


std::vector<std::tuple<double,double,int> > create_index(Rcpp::IntegerVector group, 
                                                         Rcpp::NumericVector order) {
  std::vector<double> g;
  std::vector<double> o;
  g = as<std::vector<double> >(group);
  o = as<std::vector<double> >(order);
  return(index_and_sort<double,double>(g,o));
}


bool search_group_l ( const std::tuple<double,double,int>& vec, const double& val) {return (val > std::get<0>(vec));}
bool search_group_u ( const double& val, const std::tuple<double,double,int>& vec) {return (val < std::get<0>(vec));}

bool search_order_l ( const std::tuple<double,double,int>& vec, const double& val) {return (val > std::get<1>(vec));}
bool search_order_u ( const double& val, const std::tuple<double,double,int>& vec) {return (val < std::get<1>(vec));}

std::pair<int,int> find_positions(std::vector<std::tuple<double,double,int> >& xi, const double min_g, const double max_g, const double min_o, const double max_o) {
  std::vector<std::tuple<double,double,int> >::iterator first_group_low, first_group_up, first, last_group_low, last_group_up, last;

  first_group_low = std::lower_bound(xi.begin(), xi.end(), min_g, search_group_l);
  first_group_up = std::upper_bound(xi.begin(), xi.end(), min_g, search_group_u);
  first = std::lower_bound(first_group_low, first_group_up, min_o, search_order_l);
  
  last_group_low = std::lower_bound(xi.begin(), xi.end(), max_g, search_group_l);
  last_group_up = std::upper_bound(xi.begin(), xi.end(), max_g, search_group_u);
  last = std::lower_bound(last_group_low, last_group_up, max_o, search_order_l);
  
  int firstv = first - xi.begin();
  int lastv = (last - xi.begin());
  // note that last isn't strictly the position, it's the position after the last position.
  
  std::pair<int,int> out(firstv, lastv);
  return(out);
}


/*** R
*/
