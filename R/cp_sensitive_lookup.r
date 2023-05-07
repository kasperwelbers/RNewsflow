cp_idf_matrix <- function(dtm, min_docfreq=1, max_docprob=1, min_obs_exp=0) {
  simmat = term_cooccurence_docprob(dtm, max_docfreq = max_docprob * nrow(dtm), min_obs_exp=min_obs_exp,
                                    min_docfreq=min_docfreq)
  
  simmat = Matrix::drop0(simmat)
  simmat = Matrix::tril(simmat) ## only the lower triangle is used in the cp_lookup similarity measure
  
  simmat = methods::as(methods::as(simmat, 'generalMatrix'), 'TsparseMatrix')
  simmat@x = log(nrow(dtm) / simmat@x)
  
  term_idf = diag(simmat)
  cooc_indices = summary(simmat)[,c('i','j')]
  simmat@x = simmat@x - (term_idf[cooc_indices$i] + term_idf[cooc_indices$j])
  diag(simmat) = term_idf
  
  simmat = Matrix::drop0(simmat)
  methods::as(simmat, 'CsparseMatrix')
}

prepare_cp_lookup_matrix <- function(m, m2, idf_from = c('m','m2','both')) {
  if (identical(m,m2)) {
    simmat = cp_idf_matrix(m)
  } else {
    if (idf_from == 'm') {
      simmat = cp_idf_matrix(m)  
    } else {
      terms = colnames(m)
      if (!identical(terms, colnames(m2))) 
        m2 = as(reindexTerms(m2, terms), 'CsparseMatrix')
      
      if (idf_from == 'both')
        simmat = cp_idf_matrix(rbind(m, m2))
      if (idf_from == 'm2')
        simmat = cp_idf_matrix(m2)
      
      ## also filter out all term combinations that do not occur in m
      nz_in_m = crossprod(m) > 0
      simmat = simmat * nz_in_m
    }
    
  }
  simmat
}
