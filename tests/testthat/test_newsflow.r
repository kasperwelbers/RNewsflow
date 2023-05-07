testthat::context('Using newsflow')


test_that("rnewsflow", {
  library(RNewsflow)
  library(quanteda)
  options(Matrix.warnDeprecatedCoerce = 2)
  
  
  dtm = rnewsflow_dfm
  docvars(dtm)$ym = docvars(dtm)$date
  docvars(dtm)$date = NULL
  
  el1 = compare_documents(dtm, date_var = 'ym', min_similarity = 0.1)
  testthat::expect_equal(nrow(el1$d), 50518)
  el2 = compare_documents(dtm, date_var = 'ym', min_similarity = 0.5)
  testthat::expect_equal(nrow(el2$d), 4238)

  g = as_document_network(el1)
  
  ## meta with margin attributes
  test = compare_documents(dtm, date_var = 'ym', min_similarity = 0.5)
  testthat::expect_true(all(c('from_n','from_sum','lag_n','lag_sum') %in% colnames(test$from_meta)))
  testthat::expect_true(all(c('to_n','to_sum') %in% colnames(test$to_meta)))
  
  ## test deduplication
  x = delete_duplicates(dtm, date_var='ym', similarity = 0.9, tf_idf = T)
  
  
  ## test if lag_attr is correct
  tokens = quanteda::tokens(c('a b c d e', 'b c d f', 'a b d e', 'b c f', 'f g h'))
  m = quanteda::dfm(tokens)
  
  docvars(m, 'date') = seq.POSIXt(as.POSIXct('2010-01-01'), as.POSIXct('2010-01-05'), by='day')
  docvars(m, 'date')[2] = docvars(m, 'date')[1]
  g = newsflow_compare(m, hour_window=c(-1000,1000), only_complete_window = F)
  testthat::expect_equal(sum(V(g)$lag_sum), 
                         sum(E(g)$weight[E(g)$hourdiff < 0]))
  
})

