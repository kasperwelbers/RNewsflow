#' The tcrossprod function for sparse matrices, for people that like drowning in parameters
#'
#' tcrossprod with benefits. Enables limiting row combinations to within specified groups 
#' and date windows, and filters results that do not pass the threshold on the fly.
#' To achieve this, options for similarity measures are included in the function.
#' For example, to get the cosine similarity, you can normalize with "l2" and use the "prod" (product) function for the   
#'
#' This function is called by the document comparison functions (documents.compare, newsflow.compare, delete.duplicates).
#' We only expose it here for additional flexibility, and because it could be usefull outside of the purpose of this package. 
#'
#' @param m A dgCMatrix
#' @param m2 A dgCMatrix
#' @param min_value Optionally, a numerical value, specifying the threshold for including a score in the output. 
#' @param max_value Optionally, a numerical value for the upper limit for including a score in the output.
#' @param only_upper if true, only the upper triangle of the matrix is returned. Only possible for symmetrical output (m and m2 have same number of columns)
#' @param diag if false, the diagonal of the matrix is not returned. Only possible for symmetrical output (m and m2 have same number of columns)
#' @param top_n an integer, specifying the top number of strongest scores for each column in m
#' @param rowsum_div if true, divide crossproduct by column sums of m. (this has to happen within the loop for min_value and top_n filtering)
#' @param zscore if true, transform the similarity scores for each row in m to z-scores. The min_value and max_value filters will in that case apply to the zscore value.
#' @param normalize normalize rows by a given norm score. Default is 'none' (no normalization). 'l2' is the l2 norm (use in combination with 'prod' crossfun for cosine similarity). 'l2soft' is the adaptation of l2 for soft similarity (use in combination with 'softprod' crossfun for soft cosine) 
#' @param crossfun The function used in the vector operations. 
#'                 Normally this is the "prod", for product (dot product). 
#'                 Here we also allow the "min", for minimum value. 
#'                 We use this in our document overlap_pct score.
#'                 In addition, there is the (experimental) softprod, that can be used in combination with softl2 normalization to get the soft cosine similarity.
#'                 And, the "maxproduct" is a special case used in the query_lookup measure, that uses product but only returns the score of the strongest matching term. 
#' @param group Optionally, a character vector that specifies a group (e.g., source) for each row in m. If given, only pairs of rows with the same group are calculated. 
#' @param group2 If m2 and group are used, group2 has to be used to specify the groups for the rows in m2 (otherwise group will be ignored)
#' @param date Optionally, a character vector that specifies a date for each row in m. If given, only pairs of rows within a given date range (see lwindow, rwindow and date_unit) are calculated. 
#' @param date2 If m2 and date are used, date2 has to be used to specify the date for the rows in m2 (otherwise date will be ignored)
#' @param lwindow If date (and date2) are used, lwindow determines the left side of the date window. e.g. -10 means that rows are only matched with rows for which date is within 10 [date_units] before.
#' @param rwindow Like lwindow, but for the right side. e.g. an lwindow of -1 and rwindow of 1, with date_unit is "days", means that only rows are matched for which the dates are within a 1 day distance
#' @param date_unit The date unit used in lwindow and rwindow. Supports "days", "hours", "minutes" and "seconds". Note that refers to the time distance between two rows ("days" doesn't refer to calendar days, but to a time of 24 hours)
#' @param simmat if softcos is used, a symmetric matrix with terms that indicates the similarity of terms (i.e. adjacency matrix). If NULL, a cosine similarity matrix will be created on the go 
#' @param simmat_thres if softcos is used, a threshold for the term similarity. 
#' @param batchsize If group and/or date are used, size of batches.
#' @param verbose if TRUE, report progress
#'
#' @return A dgCMatrix
#' @export
#'
#' @examples
#' set.seed(1)
#' m = Matrix::rsparsematrix(5,10,0.5)
#' tcrossprod_sparse(m, min_value = 0, only_upper = FALSE, diag = TRUE)
#' tcrossprod_sparse(m, min_value = 0, only_upper = FALSE, diag = FALSE)
#' tcrossprod_sparse(m, min_value = 0, only_upper = TRUE, diag = FALSE)
#' tcrossprod_sparse(m, min_value = 0.2, only_upper = TRUE, diag = FALSE)
#' tcrossprod_sparse(m, min_value = 0, only_upper = TRUE, diag = FALSE, top_n = 1)
tcrossprod_sparse <- function(m, m2=NULL, min_value=NULL, max_value=NULL, only_upper=F, diag=T, top_n=NULL, rowsum_div=F, zscore=F, normalize=c('none','l2','softl2'), crossfun=c('prod','min','softprod','maxproduct'), group=NULL, group2=NULL, date=NULL, date2=NULL, lwindow=-1, rwindow=1, date_unit=c('days','hours','minutes','seconds'), simmat=NULL, simmat_thres=NULL, batchsize=1000, verbose=F) {
  date_unit = match.arg(date_unit)
  crossfun = match.arg(crossfun)
  normalize = match.arg(normalize)
  
  if (crossfun == 'min' && min(m) < 0) stop('The "min" crossfun cannot be used if the dtm contains negative values')
  if (is.null(top_n)) top_n = 0
  if (is.null(m2)) {
    m2 = m
    group2 = group
    date2 = date
  } else {
    if (crossfun == 'min' && min(m2) < 0) stop('The "min" crossfun cannot be used if the dtm contains negative values')
  }
  if (is.null(colnames(m))) colnames(m) = 1:ncol(m) ##     mainly for testing (normally there should really be column names)
  if (is.null(colnames(m2))) colnames(m2) = 1:ncol(m2) ## 
  
  if (!is.null(group) && !is.null(group2)) {
    if (!length(group) == nrow(m)) stop('group has to have the same length as the number of rows in m')
    if (!length(group2) == nrow(m2)) stop('group2 has to have the same length as the number of rows in m2')
    group = as.character(group)
    group2 = as.character(group2)
    unique_group = unique(c(group,group2))
    group = match(group, unique_group)
    group2 = match(group2, unique_group)
  } else {
    group = rep(1, nrow(m))
    group2 = rep(1, nrow(m2))
  }
  
  if (!is.null(date) && !is.null(date2)) {
    if (!length(date) == nrow(m)) stop('date has to have the same length as the number of rows in m')
    if (!length(date2) == nrow(m2)) stop('date2 has to have the same length as the number of rows in m2')
    
    date = as.POSIXct(date)
    date2 = as.POSIXct(date2)
    order1 = as.numeric(date)
    order2 = as.numeric(date2)
    startorder = min(c(order1, order2))
    order1 = order1 - startorder
    order2 = order2 - startorder
    if (date_unit == 'seconds') unit_multip = 1
    if (date_unit == 'minutes') unit_multip = 60
    if (date_unit == 'hours') unit_multip = 60 * 60    
    if (date_unit == 'days') unit_multip = 60 * 60 * 24
    lwindow = lwindow * unit_multip
    rwindow = rwindow * unit_multip
    
  } else {
    order1 = rep(1, nrow(m))
    order2 = rep(1, nrow(m2))
  }
  
  if (is.null(min_value)) {
    use_min=FALSE
    min_value=0 ## not used, but can't be NULL
  } else use_min=TRUE
  
  if (is.null(max_value)) {
    use_max=FALSE
    max_value=1 ## not used, but can't be NULL
  } else use_max=TRUE
  
  
  if (crossfun == 'softprod' || normalize == 'l2soft') {
    if (!is.null(simmat)) {
      if (!nrow(simmat) == ncol(simmat)) stop('simmat has to be symmetrical')
      if (!identical(colnames(simmat), colnames(m))) stop('colnames(m) has to be identical to colnames(simmat)')
      if (!is.null(simmat_thres)) m[m < simmat_thres] = 0
      simmat = methods::as(simmat, 'dgCMatrix')
    } else {
      if (identical(m,m2)) {
        simmat = tcrossprod_sparse(t(m), normalize='l2', min_value = simmat_thres, diag=T)
      } else {
        if (!identical(colnames(m), colnames(m2))) {
          terms = union(colnames(m), colnames(m2))
          m = methods::as(reindexTerms(m, terms), 'dgCMatrix')
          m2 = methods::as(reindexTerms(m2, terms), 'dgCMatrix')
        }
        simmat = tcrossprod_sparse(t(rbind(m,m2)),normalize='l2', min_value = simmat_thres, diag=T)
      }
    }
  } else {
    simmat = methods::as(Matrix::spMatrix(0,0), 'dgCMatrix')
  }
  
  cp = batched_tcrossprod_cpp(m, m2, group1=group, group2=group2, order1=order1, order2=order2, simmat=simmat, 
                              use_min=use_min, min_value=min_value, use_max=use_max, max_value=max_value, 
                              top_n=top_n, diag=diag, only_upper=only_upper, 
                              rowsum_div=rowsum_div, zscore=zscore, normalize=normalize, crossfun=crossfun,
                              lwindow=lwindow, rwindow=rwindow, verbose=verbose, batchsize=batchsize)
  rownames(cp) = rownames(m)
  colnames(cp) = rownames(m2)
  cp
}

