testthat::context('Using newsflow')


test_that("rnewsflow", {
  library(RNewsflow)
  library(quanteda)
  
  dtm = rnewsflow_dfm
  
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.1)
  testthat::expect_equal(igraph::ecount(test), 50518)
  test = newsflow.compare(dtm, date.var = 'date', min.similarity = 0.5)
  testthat::expect_equal(igraph::ecount(test), 4238)
  
  x = delete.duplicates(dtm, date.var='date', similarity = 0.9, tf.idf = T)
})

