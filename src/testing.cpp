#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

bool search_l ( const double& vec, const double& val) {return (val > vec);}
bool search_r ( const double& val, const double& vec) {return (val < vec);}

// [[Rcpp::export]]
void findit(std::vector<double> x, double start, double end, double v1, double v2) {
  const std::vector<double>::iterator lowg = std::lower_bound(x.begin(), x.end(), start, search_l);
  const std::vector<double>::iterator highg = std::upper_bound(x.begin(), x.end(), end, search_r);
  
  Rcout << (lowg == x.begin()) << " " << (highg == x.end()) << std::endl;
  
  const std::vector<double>::iterator low = std::lower_bound(lowg, highg, v1, search_l);
  const std::vector<double>::iterator high = std::upper_bound(lowg, highg, v2, search_r);
  
  double lowv = low - x.begin();
  double highv = high - x.begin();
  if (high == x.end()) highv--;
  Rcout << lowv << " " << highv << std::endl;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x = c(1,2,3,3,3,3,3,5,5,5,6,6,7,7)
length(x)
findit(x, 3, 3, 2, 5)
*/
