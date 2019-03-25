#ifndef __WINDOW_ITERATOR_H_INCLUDED__
#define __WINDOW_ITERATOR_H_INCLUDED__

#include <RcppEigen.h>
#include <Rcpp.h>
#include "RNewsflow_types.h"

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

class windowIterator {
  Index index1, index2;
  int lwindow, rwindow, lwindow_border, rwindow_border;
  int ls, le, rs, re; // left start, left end, right start, right end
public:
  std::vector<double> leftsum, rightsum;
  double lefttotal, righttotal;
  int pos = 0;
  bool done = false;
  windowIterator(Index i1, Index i2, int l, int r, int lwb = 0, int rwb = 0): index1(i1), index2(i2), lwindow(l), rwindow(r), lwindow_border(lwb), rwindow_border(rwb) {
    if (l > 0) stop("left window has to be zero or lower");
    if (lwb > 0) stop("left window border has to be zero or lower");
    if (r < 0) stop("right window has to be zero or higher");
    if (rwb < 0) stop("right window border has to be zero or higher");
  }
  void start(const SpMat& m);
  void increment(const SpMat& m);
};



#endif