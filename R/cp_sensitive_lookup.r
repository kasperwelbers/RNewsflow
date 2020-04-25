cp_idf_matrix <- function(dtm, min_docfreq=1, max_docprob=1, min_obs_exp=0) {
  #dtm = dtm > 0
  #simmat = Matrix::crossprod(dtm)
  #simmat@x[simmat@x / nrow(dtm) > 0.01] = 0
  
  simmat = term_cooccurence_docprob(dtm, max_docfreq = max_docprob * nrow(dtm), min_obs_exp=min_obs_exp,
                                    min_docfreq=min_docfreq)
  
  #diag(simmat) = 0   ## !! disables single terms
  simmat = Matrix::drop0(simmat)
  simmat = Matrix::tril(simmat) ## only the lower triangle is used in the cp_lookup similarity measure
  
  simmat = as(simmat, 'dgTMatrix')
  simmat@x = log(nrow(dtm) / simmat@x)
  
  #n = nrow(dtm)
  #weightfun <- function(x) log((n - x + 0.5) / (x+0.5), base=2)^2
  #weightfun <- function(x) log((n - x + 0.5) / (x+0.5))^2
  #simmat@x = weightfun(simmat@x) / weightfun(1)

  #weightfun <- function(x, base=1.02) base/(base^x)
  #simmat@x = weightfun(simmat@x, base=1.02)
  #simmat@x[simmat@x < 0.0001] = 0
  #simmat = Matrix::drop0(simmat)
  
  term_idf = diag(simmat)
  cooc_indices = summary(simmat)[,c('i','j')]
  simmat@x = simmat@x - (term_idf[cooc_indices$i] + term_idf[cooc_indices$j])
  ###simmat@x[simmat@x < 0] = 0
  diag(simmat) = term_idf
  
  
  simmat = Matrix::drop0(simmat)
  #simmat = simmat^2
  methods::as(simmat, 'dgCMatrix')
}

prepare_cp_lookup_matrix <- function(m, m2, id_from = c('m','m2','both')) {
  if (identical(m,m2)) {
    simmat = cp_idf_matrix(m)
  } else {
    if (id_from == 'm') {
      simmat = cp_idf_matrix(m)  
    } else {
      terms = colnames(m)
      if (!identical(terms, colnames(m2))) 
        m2 = methods::as(reindexTerms(m2, terms), 'dgCMatrix')
      
      if (id_from == 'both')
        simmat = cp_idf_matrix(rbind(m, m2))
      if (id_from == 'm2')
        simmat = cp_idf_matrix(m2)
      
      ## also filter out all term combinations that do not occur in m
      nz_in_m = crossprod(m) > 0
      simmat = simmat * nz_in_m
    }
    
  }
  simmat
}


function() {
  dtm = rnewsflow_dfm
  #RNewsflow:::cp_idf_matrix(dtm[,1:10])
  test = tcrossprod_sparse(dtm, min_value=0, crossfun = 'cp_lookup_norm', verbose=T)
  test
  max(test)
  x = Matrix::rsparsematrix(10,10,0.1)
  y = Matrix::rsparsematrix(10,10,0.1)
  max(test)
  hist(test@x)
  length(test@x)
  sum(cp_mat)
  
  n = 100
  x = 5
  y = 7
  z = 2
  xy = 5
  xz = 1
  
  
  
  log(n/x)
  log(n/y)
  q = log(n/xy) - (log(n/x) + log(n/y))
  
  q + log(n/x) + log(n/y)
  log(n/xy)
  
  
  q1 = log(n/xy) - (log(n/x) + log(n/y))
  q2 = log(n/xz) - (log(n/x) + log(n/z))
  
  log(n/x) + log(n/y) + log(n/z) + q1 + q2
  
  log(n/x) + log(n/y) + q1 + log(n/x) + log(n/z) + q2
  
  log(n/xy) + log(n/xz)
  
}