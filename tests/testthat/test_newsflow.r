testthat::context('Using newsflow')


test_that("rnewsflow", {
  library(RNewsflow)
  library(quanteda)
  
  dtm = rnewsflow_dfm
  
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.1)
  testthat::expect_equal(igraph::ecount(test), 50518)
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.5)
  testthat::expect_equal(igraph::ecount(test), 4238)
  
  ## test meta
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.5)
  va = vertex.attributes(test)
  names(va)
  
  ## meta with margin attributes
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.5)
  va = vertex.attributes(test)
  testthat::expect_true(all(c('from_n','from_sum','lag_n','lag_sum') %in% names(va)))
  testthat::expect_true(all(c('to_n','to_sum') %in% names(va)))
  

  ## test matrix output
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.5, return_as = 'matrix')
  testthat::expect_s4_class(test, 'dgCMatrix')
  testthat::expect_true(all(c('row_n','row_sum','lag_n','lag_sum') %in% names(attr(test, 'row_meta'))))
  testthat::expect_true(all(c('col_n','col_sum') %in% names(attr(test,'col_meta'))))
  
  ## test edgelist
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.5, return_as = 'edgelist')
  testthat::expect_s3_class(test, 'data.table')
  testthat::expect_true(all(c('from_n','from_sum','lag_n','lag_sum') %in% names(attr(test, 'from_meta'))))
  testthat::expect_true(all(c('to_n','to_sum') %in% names(attr(test,'to_meta'))))
  
  ## test deduplication
  x = delete.duplicates(dtm, date.var='date', similarity = 0.9, tf.idf = T)
  
  
  ## test if lag_attr is correct
  m = quanteda::dfm(c('a b c d e', 'b c d f', 'a b d e', 'b c f', 'f g h'))
  docvars(m, 'date') = seq.POSIXt(as.POSIXct('2010-01-01'), as.POSIXct('2010-01-05'), by='day')
  docvars(m, 'date')[2] = docvars(m, 'date')[1]
  g = newsflow.compare(m, return_as='edgelist', hour.window=c(-1000,1000), only.complete.window = F)
  fm = attr(g, 'from_meta')
  testthat::expect_equal(sum(fm$lag_sum), sum(g$weight[g$hourdiff < 0]))
  
  
  
})

