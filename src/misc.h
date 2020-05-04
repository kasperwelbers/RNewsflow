#ifndef __MISC_H_INCLUDED__
#define __MISC_H_INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

std::vector<double> get_row_l2(Eigen::SparseMatrix<double>& m);
std::vector<double> softcos_row_mag(const Eigen::SparseMatrix<double>& m, const Eigen::SparseMatrix<double>& simmat, bool verbose);
  
  
bool search_group_l ( const std::tuple<double,double,int>& vec, const double& val) ;
bool search_group_u ( const double& val, const std::tuple<double,double,int>& vec) ;
bool search_order_l ( const std::tuple<double,double,int>& vec, const double& val) ;
bool search_order_u ( const double& val, const std::tuple<double,double,int>& vec) ;
std::pair<int,int> find_positions(std::vector<std::tuple<double,double,int> >& xi, const double min_g, const double max_g, const double min_o, const double max_o) ;

template <typename T1, typename T2>
bool sort_tuple0( const std::tuple<T1,T2,int>& a, const std::tuple<T1,T2,int>& b) {return std::get<0>(a) < std::get<0>(b);}

template <typename T1, typename T2>
bool sort_tuple1( const std::tuple<T1,T2,int>& a, const std::tuple<T1,T2,int>& b) {return std::get<1>(a) < std::get<1>(b);}

template <typename T>
bool sort_pair_low ( const std::pair<T,int>& a, const std::pair<T,int>& b) {return a.first < b.first;}

template <typename T>
bool sort_pair_high ( const std::pair<T,int>& a, const std::pair<T,int>& b) {return a.first > b.first;}

template <typename T1, typename T2>
std::vector<std::tuple<T1,T2,int> > index_and_sort(std::vector<T1>& x, std::vector<T2>& y) {
  if (!(x.size() == y.size())) stop("x and y are not of the same length");
  std::vector<std::tuple<T1,T2,int> > xi(x.size());
  
  for (int i = 0; i < x.size(); i++) {
    xi[i] = std::tuple<T1,T2,int>(x[i], y[i], i);
  }
  std::sort(xi.begin(), xi.end(), sort_tuple1<T1,T2>);
  std::stable_sort(xi.begin(), xi.end(), sort_tuple0<T1,T2>);
  return(xi);
}

template <typename T>
std::vector<std::pair<T,int> > index_and_sort_top_n(std::vector<T>& x, int top_n, int offset) {
  // partial sort to get top_n highest scores, but after creating pair to remember the original indices.
  std::vector<std::pair<T,int> > xi(x.size());
  for (int i = 0; i < x.size(); i++) {
    xi[i] = std::pair<T,int>(x[i], i+offset);
  }
  partial_sort(xi.begin(), xi.begin()+top_n, xi.end(), sort_pair_high<T>);
  return(xi);
}

Eigen::SparseMatrix<double> batch_simmat_prepare(Eigen::SparseMatrix<double>& m2_batch, const Eigen::SparseMatrix<double>& simmat);

Eigen::SparseMatrix<double> sm_prepare(Eigen::SparseMatrix<double>& m, std::vector<std::tuple<double,double,int> > index, const Eigen::SparseMatrix<double>& simmat, bool transpose, std::string normalize);

std::vector<std::tuple<double,double,int> > create_index(Rcpp::IntegerVector group, 
                                                         Rcpp::NumericVector order);

double sum_std_vec(std::vector<double>& x);
double nz_std_vec(std::vector<double>& x);
  
#endif