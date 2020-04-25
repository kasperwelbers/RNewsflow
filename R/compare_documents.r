#' Compare the documents in a dtm
#' 
#' This function calculates document similarity scores using a vector space approach. The most important
#' benefit is that it includes options for limiting the number of comparisons that need to be made and filtering
#' the results, that are efficiently implemented in a custom inner product calculation. This makes it possible 
#' to compare a huge number of documents, especially for cases where only documents witihin a given time window need
#' to be compared.
#' 
#' By default, the function performs a regular tcrossprod of the dtm (with itself or with dtm_y). The following parameters can be
#' set to limit comparisons and filter output:
#' \itemize{
#'    \item{If the 'date_var' is specified. The given hour_window is used to only compare documents within the specified time distance.}
#'    \item{If the 'group_var' is specified, only documents for which the group is identical will be compared.}
#'    \item{With the 'min_similarity' argument, the output can be filtered with a minimum similarity threshold. For the inner product of two
#'          DTMs the size of the output matrix is often the main bottleneck for comparing many documents, because it generally increases exponentially with
#'          the number of documents in the DTMs. Even a low similarity threshold can greatly reduce the size of the output}
#'    \item{As an alternative or additional filter, you can limit the results for each row in dtm to the highest top_n similarity scores}
#' }
#'  
#' Margin attributes are also included in the output in the from_meta and to_meta data.tables (see details).
#' If copy_meta = TRUE, The dtm docvars are also included in from_meta and to_meta.
#' 
#' @param dtm         A quanteda \link[quanteda]{dfm}. Note that it is common to first weight the dtm(s) before calculating document similarity,
#'                    For this you can use quanteda's \link[quanteda]{dfm_tfidf} and \link[quanteda]{dfm_weight} 
#' @param dtm_y       Optionally, another dtm. If given, the documents in dtm will be compared to the documents in dtm_y. 
#' @param date_var    Optionally, the name of the column in docvars that specifies the document date. The values should be of type POSIXct, or coercable with as.POSIXct.
#'                    If given, the hour_window argument is used to only compare documents within a time window. 
#' @param hour_window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. 
#'                    For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#'                    It is possible to specify time windows down to the level of seconds by using fractions (hours / 60 / 60).
#' @param group_var   Optionally,  The name of the column in docvars that specifies a group (e.g., source, sourcetype). If given, only documents within the same group will be compared.
#' @param measure     The measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), the assymetrical measures "overlap_pct" (percentage of term scores in the document 
#'                    that also occur in the other document), "overlap" (like overlap_pct, but as the sum of overlap instead of the percentage) and the symmetrical soft cosine measure (experimental).
#'                    The regular crossprod (inner product) is also supported.
#'                    If the dtm's are prepared with the create_queries function, the special "query_lookup" and "query_lookup_pct" can be used.
#' @param tf_idf      If TRUE, weigh the dtm (and dtm_y) by term frequency - inverse document frequency. For more control over weighting,
#'                    we recommend using quanteda's \link[quanteda]{dfm_tfidf} or \link[quanteda]{dfm_weight} on dtm and dtm_y. 
#' @param min_similarity A threshold for similarity. lower values are deleted. For all available similarity measures zero means no similarity.
#' @param n_topsim    An alternative or additional sort of threshold for similarity. Only keep the [n_topsim] highest similarity scores for x. Can return more than [n_topsim] similarity scores in the case of duplicate similarities.
#' @param only_complete_window If True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' @param copy_meta   If TRUE, copy the dtm docvars to the from_meta and to_meta data.tables
#' @param backbone_p  Apply backbone filtering with a "disparity" filter (see \href{https://www.pnas.org/content/106/16/6483.full}{Serrano et al.}).
#'                    It is different from the original disparity filter algorithm in that it only looks at outward edges. Also, the outward degree k is
#'                    measured as all possible edges (within a window), not just the non-zero edges.
#' @param simmat      If softcosine is used, a symmetrical matrix with the similarity scores of terms. If NULL, the cosine similarity of terms in dtm will be used
#' @param simmat_thres A large, dense simmat can lead to memory problems and slows down computation. A pragmatig (though not mathematically pure) solution is to use a threshold to prune small similarities. 
#' @param verbose     If TRUE, report progress
#' 
#' @details 
#' Margin attributes are added to the meta data.
#' The reason for including this is that some values that are normally available in a similarity matrix are missing if certain filter options are used.
#' If group or date is used, we don't know how many columns a rows has been compared to (normally this is all columns).
#' If a min/max or top_n filter is used, we don't know the true row sums (and thus row means).
#' The meta data therefore includes the "row_n", "row_sum", "col_n", and "col_sum".
#' In addition, there are "lag_n" and "lag_sum". this is a special case where row_n and row_sum are calculated for only matches where the column date < row date (lag).
#' This can be used for more refined calculations of edge probabilities before and after a row document.
#' 
#' @return A S3 class for RNewsflow_edgelist, which is a list with the edgelist, from_meta and to_meta data.tables.
#' @export
#' 
#' @examples 
#' dtm = quanteda::dfm_tfidf(rnewsflow_dfm)
#' el = compare_documents(dtm, date_var='date', hour_window = c(0.1, 36))
#' 
#' 
#' d = data.frame(text = c('a b c d e', 
#'                         'e f g h i j k',
#'                         'a b c'),
#'                date = as.POSIXct(c('2010-01-01','2010-01-01','2012-01-01')), 
#'                stringsAsFactors=FALSE)
#' corp = quanteda::corpus(d, text_field='text')
#' dtm = quanteda::dfm(corp)
#' 
#' g = compare_documents(dtm)
#' g
#' 
#' g = compare_documents(dtm, measure = 'overlap_pct')
#' g
compare_documents <- function(dtm, dtm_y=NULL, date_var=NULL, hour_window=c(-24,24), group_var=NULL, 
                              measure=c('cosine','overlap_pct','overlap','crossprod','softcosine','query_lookup','query_lookup_pct','cp_lookup','cp_lookup_norm'), tf_idf=F,
                              min_similarity=0, n_topsim=NULL, only_complete_window=T, copy_meta=F,
                              backbone_p=1, simmat=NULL, simmat_thres=NULL, verbose=FALSE){
    
    measure = match.arg(measure)
    batchsize = 1000
  
    ########### prepare dtm
    if (!methods::is(dtm, 'dfm')) stop('dtm has to be a quanteda dfm')
    if (tf_idf) dtm = quanteda::dfm_tfidf(dtm)
    meta = quanteda::docvars(dtm)
    meta$document_id = rownames(dtm)
    date = get_date(meta, date_var)
    group = get_group(meta, group_var)

    ########### prepare dtm_y
    if (!is.null(dtm_y)) {
      if (!methods::is(dtm_y, 'dfm')) stop('dtm_y has to be a quanteda dfm')
      if (tf_idf) dtm_y = quanteda::dfm_tfidf(dtm_y)
      if (!identical(quanteda::featnames(dtm), quanteda::featnames(dtm_y))) 
        dtm_y = quanteda::dfm_match(dtm_y, quanteda::featnames(dtm))
      
      meta_y = quanteda::docvars(dtm_y)
      meta_y$document_id = rownames(dtm_y)
      date_y = get_date(meta_y, date_var)
      if (is.null(date) || is.null(date_y)) date = date_y = NULL
      group_y = get_group(meta_y, group_var)
    } else {
      dtm_y = NULL
      meta_y = NULL
      date_y = NULL
      group_y = NULL
    }
  
    
    ####### inspect and filter date range
    if (!is.null(date)) {
      mindate = if (!is.null(date_y)) min(date_y) else min(date)
      maxdate = if (!is.null(date_y)) max(date_y) else max(date)
      if (only_complete_window) {
        mindate = mindate - as.difftime(hour_window[1], units = 'hours')
        maxdate = maxdate - as.difftime(hour_window[2], units = 'hours')
      }
      left_out = date < mindate
      right_out = date > maxdate
      if (any(left_out) || any(right_out)) {
        keep = !left_out & !right_out
        
        message(sprintf('Dropped %s rows from dtm that do not have a complete window:
    dtm date range:          %s - %s
    comparison window range: %s - %s',
                        nrow(dtm)-sum(keep), min(date), max(date), mindate, maxdate))
        if (!any(keep)) stop('No possible comparisons remain with current settings')
        
        dtm = dtm[keep,]
        meta$.complete_window = keep
        #meta = meta[keep,]
        date = date[keep]
        if (!is.null(group)) group = group[keep]
      }
    }
    
    ############ compare
    
    diag = !is.null(dtm_y)
    lag_attr = !is.null(date) && hour_window[1] < 0
    
    if (measure == 'cosine')       cp = tcrossprod_sparse(dtm, dtm_y, max_p=backbone_p, pvalue="disparity", normalize='l2', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2=group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure == 'overlap_pct')  cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = T, max_p=backbone_p, pvalue="disparity", crossfun = 'min', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure == 'overlap')      cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = F, max_p=backbone_p, pvalue="disparity", crossfun = 'min', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure == 'crossprod')    cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = F, max_p=backbone_p, pvalue="disparity", crossfun = 'prod', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure == 'query_lookup') cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = F, max_p=backbone_p, pvalue="disparity", crossfun = 'lookup', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure == 'query_lookup_pct') cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = T, max_p=backbone_p, pvalue="disparity", crossfun = 'lookup', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure %in% c('cp_lookup','cp_lookup_norm')) cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = F, max_p=backbone_p, pvalue="disparity", crossfun = measure, min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                              date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                              row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    if (measure == 'softcosine')   cp = tcrossprod_sparse(dtm, dtm_y, rowsum_div = T, max_p=backbone_p, pvalue="disparity", normalize='softl2', crossfun = 'softprod', min_value = min_similarity, top_n = n_topsim, diag=diag, group=group, group2 = group_y,
                                                          date = date, date2 = date_y, lwindow = hour_window[1], rwindow = hour_window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                          row_attr=T, col_attr=T, lag_attr=lag_attr, verbose=verbose)
    
    
    ## meta data is returned as data.table 
    meta = data.table::as.data.table(meta)
    if (!copy_meta) meta = subset(meta, select='document_id') 
    if (!is.null(meta_y)) {
      meta_y = data.table::as.data.table(meta_y) 
      if (!copy_meta) meta_y = subset(meta_y, select='document_id') 
    } else {
      meta_y = meta
    }
    
    data.table::setcolorder(meta, 'document_id')
    data.table::setcolorder(meta_y, 'document_id')
    
    ## add margin (col/row) attributes 
    marvars = attr(cp, 'margin')
    rowvars = grep('row\\_|lag\\_', names(marvars), value=T)
    colvars = grep('col\\_', names(marvars), value=T)
    meta_i = match(meta$document_id, rownames(cp))
    for (rowvar in rowvars) meta[[gsub('^row', 'from', rowvar)]] = marvars[[rowvar]][meta_i]
    meta_i = match(meta_y$document_id, colnames(cp))
    for (colvar in colvars) meta_y[[gsub('^col', 'to', colvar)]] = marvars[[colvar]][meta_i]
    
    cp = methods::as(cp, 'dgTMatrix')
    if (length(cp@i) == 0) return(NULL)
    
    if (!is.null(date)) {
      date = as.numeric(date)
      if (!is.null(dtm_y)) {
        date_y = as.numeric(date_y)
        hourdiff = (date_y[cp@j+1] - date[cp@i+1]) / (60*60)
      } else {
        hourdiff = (date[cp@j+1] - date[cp@i+1]) / (60*60)
      }
    } else hourdiff = NULL
    
    cp = data.table::data.table(from=rownames(cp)[cp@i+1], to=colnames(cp)[cp@j+1], weight=cp@x)
    if (!is.null(hourdiff)) cp[, hourdiff := hourdiff]
    
    cp = data.table::setorderv(cp, cols='weight', order = -1)
    cp = cp[!as.character(cp$from) == as.character(cp$to),]
    l = as_rnewsflow_edgelist(list(d = cp, 
                                   from_meta = meta, 
                                   to_meta=if (is.null(meta_y)) meta else meta_y))
    return(l)
}

as_rnewsflow_edgelist <- function(l) {
  if (!all(c('d','from_meta','to_meta') %in% names(l))) stop('input is not a proper RNewsflow_edgelist')
  if (!all(sapply(l, methods::is, 'data.table'))) stop('input is not a proper RNewsflow_edgelist')
  if (!all(c('from','to','weight') %in% colnames(l$d))) stop('input is not a proper RNewsflow_edgelist')
  if (!'document_id' %in% colnames(l$from_meta)) stop('input is not a proper RNewsflow_edgelist')
  if (!'document_id' %in% colnames(l$to_meta)) stop('input is not a proper RNewsflow_edgelist')
  class(l) = c('RNewsflow_edgelist', class(l))
  l
}

get_date <- function(meta, date_var) {
  if (any(sapply(meta, methods::is, 'POSIXlt'))) stop('meta data cannot contain a column of the POSIXlt class (long story, just use a different datetime class)')
  if (is.null(date_var)) return(NULL)
  if (!date_var %in% colnames(meta)) stop('The name specified in date_var is not a valid dfm docvar')
  if (!methods::is(meta[[date_var]], 'POSIXct')) stop("Date has to be of type POSIXct (use as.POSIXct)")
  meta[[date_var]]
}

get_group <- function(meta, group_var) {
  if (is.null(group_var)) return(NULL)
  if (!group_var %in% colnames(meta)) stop('The name specified in group_var is not a valid dfm docvar')
  meta[[group_var]]
} 