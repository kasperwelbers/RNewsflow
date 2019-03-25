#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"
#include "sim_measures.h"
#include "window_iterator.h"

void add_query(NumericVector& doc, NumericVector& term, NumericVector& term2, NumericVector& score,
               double docv, double termv, double term2v, double scorev) {
  doc.push_back(docv);
  term.push_back(termv);
  term2.push_back(term2v);
  score.push_back(scorev);
} 

//std::vector<double> create_queries(List& l, const SpMat& m, windowIterator& wi, double max_prob) {
//  int nterms = wi.leftsum.size();
//  
//  NumericVector doc, term, term2, score;
//  
//  double prob;
//  std::vector<bool> used(m.cols());
//  //for (SpMat::InnerIterator it(m, wi.pos); it; ++it) {
//  //  if (used[it.row()]) continue;
//  //  a = wi.leftsum[it.row()];
//  return(res);
//}


// [[Rcpp::export]]
List autoquery(SpMat& m1, SpMat& m2, 
                       NumericVector order1, NumericVector order2,
                       int lwindow_start=0, int lwindow_end=0, 
                       double max_prob=0.001) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
  
  Index index1, index2;
  
  IntegerVector group1(order1.size()); // dummy 
  IntegerVector group2(order2.size()); // dummy 
  index1 = create_index(group1,order1);
  index2 = create_index(group2,order2);

  SpMat simmat(1,1); // dummy, to use sm_prepare from the crossprod function (tidy up eventually)
  m1 = sm_prepare(m1, index1, simmat, true, "none");   // sort and transpose m1
  m2 = sm_prepare(m2, index2, simmat, true, "none");  // only sort m2
  
  windowIterator wi(index1, index2, lwindow_start, 0, lwindow_end);
  wi.start(m2);
  
  List l(10);
  if(Rf_isNull(l[0])) l[0] = NumericVector(10);
  NumericVector v = l[0];
  v[4] = 5;
  //Rcout << l[0] << std::endl;
  return(l);
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

x = autoquery(m1,m2,date1,date2,-3*60*60*24, -1*60*60*24)
print(x)
print(object.size(x))
*/


