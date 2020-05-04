#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"
#include "sim_measures.h"
#include "window_iterator.h"

std::vector<double> lr_chi(const SpMat& m, windowIterator& wi, double min_chi, double min_ratio, double smooth = 0.05, bool yates=true) {
  int nterms = wi.leftsum.size();
  std::vector<double> res(nterms);
   
   
  double lefttotal = wi.lefttotal + (nterms * smooth);
  double righttotal = wi.righttotal + (nterms * smooth);
  
  double a, b, c, d;
  double x, n, chi, ratio;
  for (SpMat::InnerIterator it(m, wi.pos); it; ++it) {
    a = wi.leftsum[it.row()] + smooth;
    b = wi.rightsum[it.row()] + smooth;
    c = lefttotal - a;
    d = righttotal - b;
      
    n = a+b+c+d;
    x = a*d - b*c;
    if (yates) x = abs(x) - n/2;
    chi = (n*pow(x,2)) / ((a+c) * (b+d) * (a+b) * (c+d));
    if (chi < min_chi) continue;
    ratio = (b/d) / (a/c);
    if (ratio < min_ratio) continue;
    res[it.row()] = ratio;
  }
  return(res);
}


// [[Rcpp::export]]
SpMat window_corp_comp(SpMat& m1, SpMat& m2, 
                             NumericVector order1, NumericVector order2,
                             int lwindow=0, int rwindow=0, 
                             double min_chi=0, double min_ratio=0, double smooth=1) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
  
  Index index1, index2;

  IntegerVector group1(order1.size()); // dummy 
  IntegerVector group2(order2.size()); // dummy 
  index1 = create_index(group1,order1);
  index2 = create_index(group2,order2);
  
  
  SpMat simmat(1,1); // dummy, to use sm_prepare from the crossprod function (tidy up eventually)
  m1 = sm_prepare(m1, index1, simmat, true, "none");   // sort and transpose m1
  m2 = sm_prepare(m2, index2, simmat, true, "none");  // only sort m2
  
  windowIterator wi(index1, index2, lwindow, rwindow);
  wi.start(m2);
    
  Triplet tl;
  tl.reserve(m1.nonZeros());  
  
  while (!wi.done) {
    //if (Progress::check_abort()) stop("Aborted");

    std::vector<double> res = lr_chi(m1, wi, min_chi, min_ratio, smooth);
    for (int res_i = 0; res_i < res.size(); res_i++) {
      if (res[res_i] == 0) continue;
      tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[wi.pos]), res_i, res[res_i]));
    }
    wi.increment(m2);
  }
  
  SpMat out(m1.cols(),m1.rows());
  out.setFromTriplets(tl.begin(), tl.end());
  return out;
}


/*** R
library(quanteda)
m1 = dfm(c('hond kat paard', 'aap noot mies'))
m2 = dfm(c('hond kat aap noot', 'hond kat aap noot', 'kat paard noot mies', 'aap paard mies', 'kat paard noot', 'noot mies mies', 'hond kat paard','paard paard paard paard paard mies'))

m2 = m2[,c(colnames(m1))]

dates = seq.Date(as.Date('2010-01-01'), as.Date('2010-01-10'), by=1)
dates = as.POSIXct(dates)
date1 = dates[c(5,6)]   
date2 = dates[c(1,2,3,4,7,8,9,10)]  

window_corp_comp(m1,m2,date1,date2,-3*60*60*24, 3*60*60*24)
*/


