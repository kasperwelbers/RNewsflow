#ifndef __MISC_H_INCLUDED__
#define __MISC_H_INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

std::vector<double> get_row_l2(Eigen::SparseMatrix<double>& m);
  
template <typename T1, typename T2>
bool search_l ( const std::tuple<double,double,int>& a, const std::pair<T1,T2>& b) {return (std::get<0>(a) <= b.first && std::get<1>(a) < b.second) ;}

template <typename T1, typename T2>
bool search_r ( const std::tuple<double,double,int>& a, const std::pair<T1,T2>& b) {return (std::get<0>(a) <= b.first && std::get<1>(a) <= b.second) ;}
//template <typename T1, typename T2>
//bool search_r (const std::pair<T1,T2>&a , const std::tuple<double,double,int>&  b) {return (a.first < std::get<0>(b) && a.second < std::get<1>(b) );}


template <typename T>
std::pair<int,int> find_positions(std::vector<std::tuple<double,double,int> >& xi, const T min_g, const T max_g, const T min_o, const T max_o) {
  typename std::vector<std::tuple<double,double,int> >::iterator low,up;
  
  std::pair<T,T> lv(min_g, min_o);
  std::pair<T,T> uv(max_g, max_o);
  low = std::lower_bound(xi.begin(), xi.end(), lv, search_l<T,T>);
  up = std::lower_bound(xi.begin(), xi.end(), uv, search_r<T,T>);
  
  std::pair<T,T> out(low - xi.begin(), up - xi.begin());
  return(out);
}

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
  if (!x.size() == y.size()) stop("x and y are not of the same length");
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


Eigen::SparseMatrix<double> sm_prepare(Eigen::SparseMatrix<double>& m, std::vector<std::tuple<double,double,int> > index, bool transpose, bool l2norm);

std::vector<std::tuple<double,double,int> > create_index(Rcpp::IntegerVector group, 
                                                         Rcpp::NumericVector order);
  
#endif