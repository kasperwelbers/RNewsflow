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
  for (int i=0; i < out.size(); i++) out[i] = pow(out[i], 0.5);
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


/*** R
*/
