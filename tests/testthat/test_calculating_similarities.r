testthat::context('Calculating similarities')

is_same <- function(m1,m2){
  testthat::expect_true(identical(as(as(m1, 'dgTMatrix'), 'dgCMatrix'),
                                  as(as(m2, 'dgTMatrix'), 'dgCMatrix')))
}

is_same_m1nonzero <- function(m1,m2) {
  ## only compares the nonzero values in m1 to m2 (for testing date/group and batching)
  m2[which(m1 == 0, arr.ind = T)] = 0
  m2 = Matrix::drop0(m2)
  is_same(m1,m2)
}


test_that("Matrix multiplication", {
  library(RNewsflow)
  set.seed(1)
  m = Matrix::rsparsematrix(5,10,0.5)
  
  cp = tcrossprod_sparse(m, min_value = NULL)
  cp_correct = Matrix::tcrossprod(m)
  is_same(cp,cp_correct)
  
  tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = F)
  tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F)
  tcrossprod_sparse(m, min_value = 0.2, only_upper = T, diag = F)
  tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F, top_n = 1)
  
  ## heavy lifting
  #m = abs(Matrix::rsparsematrix(1000000,20000,0.001))
  #date = seq.Date(as.Date('1000-01-01'), as.Date('4010-01-10'), by=1)[1:nrow(m)] ## just to make 1,000,000 days
  #cp = tcrossprod_sparse(m, date=date, lwindow = 15, rwindow=15, verbose=T)
  
  ## using two matrices
  
  m1 = Matrix::rsparsematrix(5,10,0.5)
  m2 = Matrix::rsparsematrix(8,10,0.5)
  cp = tcrossprod_sparse(m1, m2, min_value = NULL, only_upper = F, diag = T)
  cp_correct = Matrix::tcrossprod(m1,m2)
  is_same(cp,cp_correct)
  
  

  ## filtering by group/date. Documents have to be in the same group, or within the given date range
  
  m = Matrix::rsparsematrix(10,10,0.5)
  tcrossprod_sparse(m, group = c(1,1,1,2,2,2,3,3,3,3), batchsize = 1)
  
  date = seq.Date(as.Date('2010-01-01'), as.Date('2010-01-10'), by=1)
  tcrossprod_sparse(m, date = date, lwindow = -1, rwindow = 1)
  tcrossprod_sparse(m, date = date, lwindow = -2, rwindow = 2)
  
  cp = tcrossprod_sparse(m, date = date, lwindow = -1, rwindow = 1)
  cp_correct = Matrix::tcrossprod(m)
  is_same_m1nonzero(cp,cp_correct) 
  
  cp = tcrossprod_sparse(m, date = date, lwindow = -2, rwindow = 3)
  cp_correct = Matrix::tcrossprod(m)
  is_same_m1nonzero(cp,cp_correct) 
  
  
  ## with batches
  cp = tcrossprod_sparse(m, date = date, lwindow = -1, rwindow = 1, batchsize = 1, verbose = F)
  cp_correct = Matrix::tcrossprod(m)
  is_same_m1nonzero(cp,cp_correct) 

  ## using jacard
  m = abs(Matrix::rsparsematrix(10,10,0.5))
  m[m>0] = 1
  
  is_same(t(tcrossprod_sparse(m)) / rowSums(m),
          tcrossprod_sparse(m, rowsum_div = T))
  
  ## cosine (using l2 Norm parameter)
  m = abs(Matrix::rsparsematrix(10,10,0.5))
  cp = tcrossprod_sparse(m, l2norm=T)
  
  mnorm = sqrt(Matrix::rowSums(m^2))
  m2 = as(m, 'dgTMatrix')
  m2@x = m2@x / mnorm[m2@i+1]  
  cp_correct = Matrix::tcrossprod(m2)
  
  is_same(cp,cp_correct)
  
  ## using the min function
  m = abs(Matrix::rsparsematrix(10,10,0.5))
  
  tcrossprod_sparse(m, crossfun = 'min')
  tcrossprod_sparse(m, crossfun = 'min', rowsum_div = T) ## what percentage of words in a document also occurs in the other document?
  
})