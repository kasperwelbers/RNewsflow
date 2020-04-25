testthat::context('Lookup with cooccurrence')


test_that("Lookup with cooccurrence", {
  library(RNewsflow)
  set.seed(1)
  
  m = Matrix::rsparsematrix(10,4,0.4)
  m = (m != 0) * 1
  
  simmat = RNewsflow:::prepare_cp_lookup_matrix(m, m, id_from = 'm2')
  test = RNewsflow::tcrossprod_sparse(m, min_value=0, crossfun = 'cp_lookup', verbose=F)
  
  m[4,] & m[9,]
  should_be = sum(diag(simmat)[c(1,4)]) + simmat[4,1]
  testthat::expect_true(should_be == test[4,9])
  
  dtm = rnewsflow_dfm
  test = RNewsflow::tcrossprod_sparse(m, min_value=0, crossfun = 'cp_lookup_norm', verbose=F)
  max(test)
  
})
  