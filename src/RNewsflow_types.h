#include <RcppEigen.h>

typedef Eigen::SparseMatrix<double> SpMat;
typedef std::vector<std::tuple<double,double,int> > Index;
typedef std::vector<Eigen::Triplet<double> > Triplet;
