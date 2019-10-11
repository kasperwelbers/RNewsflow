#' Create a network of document similarities over time
#' 
#' This is a wrapper for the \code{\link{compare_documents}} function, specialised for the case of analyzing documents over time.
#' The difference is that using date_var is mandatory, and the output is returned as an igraph network (using \code{\link{as_document_network}}).
#' 
#' @param dtm         A quanteda \link[quanteda]{dfm}. Note that it is common to first weight the dtm(s) before calculating document similarity,
#'                    For this you can use quanteda's \link[quanteda]{dfm_tfidf} and \link[quanteda]{dfm_weight} 
#' @param dtm_y       Optionally, another dtm. If given, the documents in dtm will be compared to the documents in dtm_y. 
#' @param date_var    The name of the column in meta that specifies the document date. default is "date". The values should be of type POSIXct, or coercable with as.POSIXct.
#'                    If given, the hour_window argument is used to only compare documents within a time window. 
#' @param hour_window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. 
#'                    For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#'                    It is possible to specify time windows down to the level of seconds by using fractions (hours / 60 / 60).
#' @param group_var   Optionally,  The name of the column in meta that specifies a group (e.g., source, sourcetype). If given, only documents within the same group will be compared.
#' @param measure     The measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), the assymetrical measures "overlap_pct" (percentage of term scores in the document 
#'                    that also occur in the other document), "overlap" (like overlap_pct, but as the sum of overlap instead of the percentage) and the symmetrical soft cosine measure (experimental).
#'                    The regular crossprod (inner product) is also supported.
#'                    If the dtm's are prepared with the create_queries function, the special "query_lookup" and "query_lookup_pct" can be used.
#' @param tf_idf      If TRUE, weigh the dtm (and dtm_y) by term frequency - inverse document frequency. For more control over weighting,
#'                    we recommend using quanteda's \link[quanteda]{dfm_tfidf} or \link[quanteda]{dfm_weight} on dtm and dtm_y. 
#' @param min_similarity A threshold for similarity. lower values are deleted. For all available similarity measures zero means no similarity.
#' @param n_topsim    An alternative or additional sort of threshold for similarity. Only keep the [n_topsim] highest similarity scores for x. Can return more than [n_topsim] similarity scores in the case of duplicate similarities.
#' @param only_complete_window If True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' @param ...         Other arguments passed to \code{\link{compare_documents}}.
#' 
#' @return An igraph network. 
#' @export
#' 
#' @examples 
#' rnewsflow_dfm 
#' 
#' dtm = quanteda::dfm_tfidf(rnewsflow_dfm)
#' el = newsflow_compare(dtm, date_var='date', hour_window = c(0.1, 36))
newsflow_compare <- function(dtm, dtm_y=NULL, date_var='date', hour_window=c(-24,24), group_var=NULL, 
                              measure=c('cosine','overlap_pct','overlap','crossprod','softcosine','query_lookup','query_lookup_pct'), tf_idf=F,
                              min_similarity=0, n_topsim=NULL, only_complete_window=T, ...){
  
  measure = match.arg(measure)
  
  if (is.null(date_var)) stop('date_var has to be given')
  if (!'date' %in% colnames(quanteda::docvars(dtm))) stop(sprintf('date_var "%s" is not a valid column in quanteda::docvars(dtm)', date_var))
  if (!is.null(dtm_y))
    if (!'date' %in% colnames(quanteda::docvars(dtm_y))) stop(sprintf('date_var "%s" is not a valid column in quanteda::docvars(dtm_y)', date_var))
  
  el = compare_documents(dtm, dtm_y, date_var=date_var, hour_window=hour_window, group_var=group_var, 
                         measure=measure, tf_idf=tf_idf, min_similarity=min_similarity, n_topsim=n_topsim, only_complete_window=only_complete_window, ...) 
  as_document_network(el)
}

