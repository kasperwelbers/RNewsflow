reindexTerms <- function(dtm, terms){
  dtm = methods::as(dtm, 'dgTMatrix')
  documents = rownames(dtm)
  dtm = Matrix::spMatrix(nrow(dtm), length(terms), dtm@i+1, match(colnames(dtm)[dtm@j+1], terms), dtm@x)
  dimnames(dtm) = list(documents, terms)
  dtm
}


#' Compare the documents in two corpora/dtms
#' 
#' Compare the documents in corpus dtm.x with reference corpus dtm.y. 
#' 
#' The calculation of document similarity is performed using a vector space model approach. 
#' Inner-product based similarity measures are used, such as cosine similarity.
#' It is recommended to weight the DTM beforehand, for instance using Term frequency-inverse document frequency (tf.idf)
#' 
#' @param dtm A quanteda \link[quanteda]{dfm}. Alternatively, a DocumentTermMatrix from the tm package can be used.
#' @param dtm.y Optional. If given, documents from dtm will only be compared to the documents in dtm.y
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), the assymetrical measures "overlap_pct" (percentage of term scores in the document 
#'                that also occur in the other document), "overlap" (like overlap_pct, but as the sum of overlap instead of the percentage) and the symmetrical soft cosine measure (experimental).
#'                The regular crossprod (inner product) is also supported.
#'                If the dtm's are prepared with the create_queries function, the special "query_lookup" and "query_lookup_pct" can be used.
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param pvalue If used, transform the similarity score to a p-value. The value is reversed, so that higher means more similar 
#'               (and thus the min.similarity still makes sense). Currently supports "normal" and "lognormal" distribution, and the uniform distribution 
#'               used in the "disparity" filter (see \href{https://www.pnas.org/content/106/16/6483.full}{Serrano et al.}). Also "nz_normal" and "nz_lognormal" can be used
#'               to only consider the nonzero values.
#' @param simmat If softcosine is used, a symmetrical matrix with the similarity scores of terms. If NULL, the cosine similarity of terms in dtm will be used
#' @param simmat_thres If softosine is used, a threshold for the similarity scores of terms
#' 
#' @return A data frame with pairs of documents and their similarities. 
#' @export
#' 
#' @import tm
#' 
#' @examples
#' rnewsflow_dfm 
#' 
#' comp = documents.compare(rnewsflow_dfm, min.similarity=0.4)
#' head(comp)
documents.compare <- function(dtm, dtm.y=NULL, measure=c('cosine','overlap_pct','overlap','crossprod','softcosine','query_lookup','query_lookup_pct'), 
                              min.similarity=0, n.topsim=NULL, pvalue=c("none", "normal", "lognormal", "nz_normal", "nz_lognormal", "disparity"), 
                              simmat=NULL, simmat_thres=NULL) {  
  dtm = quanteda::as.dfm(dtm)
  measure = match.arg(measure)

  out = vector('list')

  if(!is.null(dtm.y)){
    dtm.y = quanteda::as.dfm(dtm.y)
    if(!all(colnames(dtm) == colnames(dtm.y))){
      ## if colnames do not match, reindex them.
      terms = unique(c(colnames(dtm), colnames(dtm.y)))
      dtm = reindexTerms(dtm, terms)
      dtm.y = reindexTerms(dtm.y, terms)
    }
  }
  dtm = methods::as(dtm, 'dgCMatrix')
  if (!is.null(dtm.y)) dtm.y = methods::as(dtm.y, 'dgCMatrix')
  
  diag = !is.null(dtm.y)
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, pvalue=pvalue, normalize='l2', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'overlap') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'crossprod') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'softcosine') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, normalize='softl2', crossfun = 'softprod', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'query_lookup') {
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, 
                           simmat=simmat, simmat_thres=simmat_thres)
  }
  if (measure == 'query_lookup_pct') {
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, 
                           simmat=simmat, simmat_thres=simmat_thres)
  }
  cp = methods::as(cp, 'dgTMatrix')
  cp = data.frame(x=rownames(cp)[cp@i+1], y=colnames(cp)[cp@j+1], similarity=cp@x)
  cp[!as.character(cp$x) == as.character(cp$y),]
}


#' Compare the documents in a dtm with a sliding window over time
#' 
#' Given a document-term matrix (DTM) with dates for each document, calculates the document similarities over time using with a sliding window.
#' 
#' The calculation of document similarity is performed using a vector space model approach. 
#' Inner-product based similarity measures are used, such as cosine similarity.
#' It is recommended to weight the DTM beforehand, for instance using Term frequency-inverse document frequency (tf.idf)
#' 
#' Meta data is included in the output. Margin attributes can also be added to meta with the margin_attr argument. see details.
#' 
#' @param dtm         A quanteda \link[quanteda]{dfm}. Alternatively, a DocumentTermMatrix from the tm package can be used, but then the meta parameter needs to be specified manually
#' @param dtm.y       Optionally, another dtm. If given, the documents in dtm will be compared to the documents in dtm.y. This cannot be combined with only.from and only.to
#' @param meta        If dtm is a quanteda dfm, docvars(meta) is used by default (meta is NULL) to obtain the meta data. Otherwise, the meta data.frame has to be given by the user, with the rows of the meta data.frame matching the rows of the dtm (i.e. each row is a document)
#' @param meta.y      Like meta, but for dtm.y (only necessary if dtm.y is used)
#' @param date.var    The name of the column in meta that specifies the document date. default is "date". The values should be of type POSIXct
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param group.var   Optionally,  The name of the column in meta that specifies a group (e.g., source, sourcetype). If given, only documents within the same group will be compared.
#' @param measure     The measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), the assymetrical measures "overlap_pct" (percentage of term scores in the document 
#'                    that also occur in the other document), "overlap" (like overlap_pct, but as the sum of overlap instead of the percentage) and the symmetrical soft cosine measure (experimental).
#'                    The regular crossprod (inner product) is also supported.
#'                    If the dtm's are prepared with the create_queries function, the special "query_lookup" and "query_lookup_pct" can be used.
#' @param min.similarity A threshold for similarity. lower values are deleted. Set to 0.1 by default.
#' @param n.topsim    An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param only.from   A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare only these documents to other documents. 
#' @param only.to     A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare other documents to only these documents.
#' @param only.complete.window If True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' @param pvalue      If used, transform the similarity score to a p-value. The value is reversed, so that higher means more similar 
#'                    (and thus the min.similarity still makes sense). Currently supports "normal" and "lognormal" distribution, and the uniform distribution 
#'                    used in the "disparity" filter (see \href{https://www.pnas.org/content/106/16/6483.full}{Serrano et al.}). Also "nz_normal" and "nz_lognormal" can be used
#'                    to only consider the nonzero values.
#' @param return_as   Detemine whether output is returned as an "edgelist", "igraph" network or sparse "matrix'.
#' @param batchsize   If group and/or date are used, size of batches.
#' @param simmat      If softcosine is used, a symmetrical matrix with the similarity scores of terms. If NULL, the cosine similarity of terms in dtm will be used
#' @param simmat_thres If softosine is used, a threshold for the similarity scores of terms
#' @param margin_attr By default, margin attributes are added to meta (see details). This can be turned of for (slightly?) faster computation and less memory usage
#' @param verbose     If TRUE, report progress
#' 
#' @details 
#' For the "igraph" output the meta data is stored as vertex attributes; for the "matrix" output as 
#' the attributes "row_meta" and "col_meta"; for the "edgelist" output as the attributes "from_meta" and "to_meta". Note that
#' attributes are removed if you perform certain operations on a matrix or data.frame, so if you want to use this information it is
#' best to assign it immediately. 
#' 
#' Margin attributes can be added to the meta data with the margin_attr argument.
#' The reason for including this is that some values that are normally available in a similarity matrix are missing if certain filter options are used.
#' If group or date is used, we don't know how many columns a rows has been compared to (normally this is all columns).
#' If a min/max or top_n filter is used, we don't know the true row sums (and thus row means).
#' margin_attr adds the "row_n", "row_sum", "col_n", and "col_sum" data to the meta data.
#' In addition, there are "lag_n" and "lag_sum". this is a special case where row_n and row_sum are calculated for only matches where the column date < row date (lag).
#' This can be used for more refined calculations of edge probabilities before and after (row_n - lag_n) a row document, which is for instance usefull for event matching.
#' 
#' @return A network/graph in the \link[igraph]{igraph} class, or an edgelist data.frame, or a sparse matrix.
#' @export
#' 
#' @examples 
#' rnewsflow_dfm 
#' 
#' dtm = quanteda::dfm_tfidf(rnewsflow_dfm)
#' g = newsflow.compare(dtm, hour.window = c(0.1, 36))
#' 
#' vcount(g) # number of documents, or vertices
#' ecount(g) # number of document pairs, or edges
#' 
#' head(igraph::get.data.frame(g, 'vertices'))
#' head(igraph::get.data.frame(g, 'edges'))
newsflow.compare <- function(dtm, dtm.y=NULL, meta=NULL, meta.y=NULL, date.var='date', hour.window=c(-24,24), group.var=NULL, measure=c('cosine','overlap_pct','overlap','crossprod','softcosine','query_lookup','query_lookup_pct'), 
                             min.similarity=0, n.topsim=NULL, only.from=NULL, only.to=NULL, only.complete.window=TRUE, pvalue = c("none", "normal", "lognormal", "nz_normal", "nz_lognormal", "disparity"), return_as = c('igraph','edgelist','matrix'), 
                             batchsize=1000, simmat=NULL, simmat_thres=NULL, margin_attr=T, verbose=FALSE){
  
  if (margin_attr) {
    ## these are actually 3 options, but we only expose margin_attr for sake of simplicity
    row_attr = T
    col_attr = T
    lag_attr = T
  }
  
  ########### prepare dtm
  if (is.null(meta)) {
    if (!methods::is(dtm, 'dfm')) stop('meta can only be NULL if dtm is a quanteda dfm class')
    meta = quanteda::docvars(dtm)
  }
  dtm = quanteda::as.dfm(dtm)
  meta$document_id = rownames(dtm)
  measure = match.arg(measure)
  return_as = match.arg(return_as)

  if (!date.var %in% colnames(meta)) stop('The name specified in date.var is not a valid dfm docvar')
  date = meta[[date.var]]
  
  if (!methods::is(date, 'POSIXct')) stop("Date has to be of type POSIXct (use as.POSIXct)")
  if (any(sapply(meta, methods::is, 'POSIXlt'))) stop('meta data cannot contain a column of the POSIXlt class')
  
  if (!is.null(group.var)) {
    if (!group.var %in% colnames(meta)) stop('The name specified in group.var is not a valid dfm docvar')
    group = meta[[group.var]]
  } else group = NULL
  
  ########### prepare dtm.y
  if (!is.null(dtm.y)) {
    if (!is.null(only.from)) stop('Cannot use only.from if dtm.y is used')
    if (!is.null(only.to)) stop('Cannot use only.to if dtm.y is used')

    if (is.null(meta.y)) {
      if (!methods::is(dtm.y, 'dfm')) stop('meta.y can only be NULL if dtm.y is a quanteda dfm class')
      meta.y = quanteda::docvars(dtm.y)
    }

    dtm.y = quanteda::as.dfm(dtm.y)
    if (!identical(quanteda::featnames(dtm), quanteda::featnames(dtm.y))) 
      dtm.y <- quanteda::dfm_match(dtm.y, quanteda::featnames(dtm))

    meta.y$document_id = rownames(dtm.y)
    measure = match.arg(measure)
    return_as = match.arg(return_as)
    
    if (!date.var %in% colnames(meta)) stop('The name specified in date.var is not a valid dfm docvar')
    date.y = meta.y[[date.var]]
    if (!methods::is(date.y, 'POSIXct')) stop("Date.y has to be of type POSIXct (use as.POSIXct)")
    if (any(sapply(meta.y, methods::is, 'POSIXlt'))) stop('meta data cannot contain a column of the POSIXlt class')
    
    if (!is.null(group.var)) {
      if (!group.var %in% colnames(meta.y)) stop('The name specified in group.var is not a valid dfm docvar')
      group.y = meta.y[[group.var]]
    } else group.y = NULL
  } else {
    dtm.y = NULL
    date.y = NULL
    group.y = NULL
    if (!is.null(only.from)) {
      if(!class(only.from) == 'logical') only.from = rownames(dtm) %in% only.from
      dtm.y = dtm
      date.y = date
      group.y = group
      dtm = dtm[only.from,]
      date = date[only.from]
      if (!is.null(group)) group = group[only.from]
    }
    if (!is.null(only.to)) {
      if(!class(only.to) == 'logical') only.to = rownames(dtm) %in% only.to
      dtm.y = dtm.y[only.to,]
      date.y = date.y[only.to]
      if (!is.null(group)) group.y = group.y[only.to]
    }
  }
  
  
  ####### inspect and filter date range
  
  mindate = if (!is.null(date.y)) min(date.y) else min(date)
  maxdate = if (!is.null(date.y)) max(date.y) else max(date)
  if (only.complete.window) {
    mindate = mindate - as.difftime(hour.window[1], units = 'hours')
    maxdate = maxdate - as.difftime(hour.window[2], units = 'hours')
  }
  left_out = date < mindate
  right_out = date > maxdate
  if (any(left_out) || any(right_out)) {
    message(sprintf('NOTE: only the rows in dtm that occur within the comparison window range are used.
      dtm date range:          %s - %s 
      comparison window range: %s - %s\n',
                    min(date), max(date), mindate, maxdate))
    dtm = dtm[!left_out & !right_out,]
    date = date[!left_out & !right_out]
    if (!is.null(group)) group = group[!left_out & !right_out]
  }

  
  ############ compare
  
  diag = !is.null(dtm.y)
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, pvalue=pvalue, normalize='l2', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2=group.y,
                                                  date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                  row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                       date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                       row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'overlap') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                       date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                   row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'crossprod') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                   date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                   row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'query_lookup') {
    if (is.null(dtm.y)) {
      dtm.y=dtm
      date.y = date
      group.y = group
    }
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                           date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                           row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  }
  if (measure == 'query_lookup_pct') {
    if (is.null(dtm.y)) {
      dtm.y=dtm
      date.y = date
      group.y = group
    }
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                           date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                           row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  }
  
  if (measure == 'softcosine') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, normalize='softl2', crossfun = 'softprod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                       date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                       row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  
  
  ## meta data is returned as data.table 
  meta = data.table::as.data.table(meta)
  if (!is.null(meta.y)) {
    meta.y = data.table::as.data.table(meta.y)
  } else {
    if (return_as %in% c('matrix','edgelist')) {
      ## if output is matrix or edgelist it is necessary to split meta and meta.x for matching the margin attributes
      ## however, for igraph is it necessary to keep a single meta if no meta.y exists to prevent creating duplicate vertex meta (with row and col attributes)
      meta.y = meta
    }  
  }
  
  
  #meta = data.table::as.data.table(meta)
  #meta.y = if (is.null(meta.y)) meta else data.table::as.data.table(meta.y)
  
  
  ## add row attributes 
  if (row_attr || col_attr || lag_attr) {
    marvars = attr(cp, 'margin')
    
    rowvars = grep('row\\_|lag\\_', names(marvars), value=T)
    colvars = grep('col\\_', names(marvars), value=T)
    if (length(rowvars) > 0) {
      meta_i = match(meta$document_id, rownames(cp))
      for (rowvar in rowvars) {
        vname = if (return_as == 'matrix') rowvar else gsub('^row', 'from', rowvar)
        meta[[vname]] = marvars[[rowvar]][meta_i]
      }
    }
    if (length(colvars) > 0) {
      if (is.null(meta.y)) {
        meta_i = match(meta$document_id, colnames(cp))
        for (colvar in colvars) {
          vname = if (return_as == 'matrix') colvar else gsub('^col', 'to', colvar)
          meta[[vname]] = marvars[[colvar]][meta_i] 
        }
      } else {
        meta_i = match(meta.y$document_id, colnames(cp))
        for (colvar in colvars) {
          vname = if (return_as == 'matrix') colvar else gsub('^col', 'to', colvar)
          meta.y[[vname]] = marvars[[colvar]][meta_i]
        }
      }
    }
  }
    
  ## return as matrix
  if (return_as == 'matrix') {
    attr(cp, 'row_meta') = meta[match(rownames(cp), meta$document_id),]
    if (is.null(meta.y)) {
      attr(cp, 'col_meta') = meta[match(rownames(cp), meta$document_id),]
    } else {
      attr(cp, 'col_meta') = meta.y[match(colnames(cp), meta.y$document_id),]
    }
    attr(cp, 'margin') = NULL
    return(cp)
  }
  
  cp = methods::as(cp, 'dgTMatrix')
  if (length(cp@i) == 0) return(NULL)
  
  ## return as edgelist
  if (return_as == 'edgelist') {
    date = as.numeric(date)
    if (!is.null(dtm.y)) {
      date.y = as.numeric(date.y)
      hourdiff = round((date.y[cp@j+1] - date[cp@i+1]) / (60*60), 3)
      #hourdiff = round(difftime(date.y[cp@j+1], date[cp@i+1], units = 'hours'),3)
    } else {
      #hourdiff = round(difftime(date[cp@j+1], date[cp@i+1], units = 'hours'),3)
      hourdiff = round((date[cp@j+1] - date[cp@i+1]) / (60*60), 3)
    }
    cp = data.table::data.table(from=rownames(cp)[cp@i+1], to=colnames(cp)[cp@j+1], weight=cp@x, hourdiff = hourdiff)
    cp = data.table::setorderv(cp, cols='weight', order = -1)
    cp = cp[!as.character(cp$from) == as.character(cp$to),]
    attr(cp, 'from_meta') = meta
    attr(cp, 'to_meta') = if (is.null(meta.y)) meta else meta.y
    return(cp)
  } 
  
  ## return as igraph network
  if (return_as == 'igraph') {
    cp = data.table::data.table(x=rownames(cp)[cp@i+1], y=colnames(cp)[cp@j+1], similarity=cp@x)
    cp = cp[!as.character(cp$x) == as.character(cp$y),]
  
    if (verbose) message('Creating network')
    
    if (!is.null(meta.y)) {
      meta = data.table::rbindlist(list(meta, meta.y), use.names = T, fill = T)
      meta = unique(meta)
    }
    g = document.network(cp, as.data.frame(meta), 'document_id', date.var)
    return(g)
  } 
}

#' Delete duplicate (or similar) documents from a document term matrix 
#' 
#' Delete duplicate (or similar) documents from a document term matrix. 
#' Duplicates are defined by: having high content similarity, occuring within a given time distance and being published by the same source.
#' 
#' Note that this can also be used to delete "updates" of articles (e.g., on news sites, news agencies). 
#' This should be considered if the temporal order of publications is relevant for the analysis. 
#' 
#' @param dtm A quanteda \link[quanteda]{dfm}. Alternatively, a DocumentTermMatrix from the tm package can be used, but then the meta parameter needs to be specified manually
#' @param meta If dtm is a quanteda dfm, docvars(meta) is used by default (meta is NULL) to obtain the meta data. Otherwise, the meta data.frame has to be given by the user, with the rows of the meta data.frame matching the rows of the dtm (i.e. each row is a document)
#' @param date.var The name of the column in meta that specifies the document date. default is "date". The values should be of type POSIXlt or POSIXct
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param group.var Optionally,  The name of the column in meta that specifies a group (e.g., source, sourcetype). If given, only documents within the same group will be compared.
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param similarity a threshold for similarity. Documents of which similarity is equal or higher are deleted
#' @param keep A character indicating whether to keep the 'first' or 'last' published of duplicate documents.
#' @param tf.idf if TRUE, weight the dtm with tf.idf before comparing documents. The original (non-weighted) DTM is returned.
#' @param dup_csv Optionally, a path for writing a csv file with the duplicates edgelist. For each duplicate pair it is noted if "from" or "to" is the duplicate, or if "both" are duplicates (of other documents)
#' @param verbose if TRUE, report progress
#' 
#' @return A dtm with the duplicate documents deleted
#' @export
#' 
#' @examples
#' ## example with very low similarity threshold (normally not recommended!)
#' dtm2 = delete.duplicates(rnewsflow_dfm, similarity = 0.5, keep='first', tf.idf = TRUE)
delete.duplicates <- function(dtm, meta=NULL, date.var='date', hour.window=c(-24,24), group.var=NULL, measure=c('cosine','overlap_pct'), similarity=1, keep='first', tf.idf=FALSE, dup_csv=NULL, verbose=F){
  measure = match.arg(measure)
  if (is.null(meta)) {
    if (!methods::is(dtm, 'dfm')) stop('meta can only be NULL if dtm is a quanteda dfm class')
    meta = quanteda::docvars(dtm)
  }
  dtm = quanteda::as.dfm(dtm)
  if (!nrow(dtm) == nrow(meta)) stop('Number of rows in dtm and meta is not the same (each row represents a document)')
  if(tf.idf) dtm = quanteda::dfm_tfidf(dtm)
  
  d = newsflow.compare(dtm, meta=meta, date.var=date.var, hour.window=hour.window, group.var=group.var, measure=measure, 
                       min.similarity=similarity, only.complete.window = F, return_as = 'edgelist', verbose=verbose)
  
  #e = igraph::get.edges(g, igraph::E(g))
  #d = igraph::get.data.frame(g, 'edges')  
  
  duplicates = c()
  if(keep == 'first') {
    duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff > 0])))
    duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff < 0])))
  }
  if(keep == 'last') {
    duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff > 0])))
    duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff < 0])))
  }
  
  ## if there are duplicate articles that occured simultaneously, delete first match to dtm rows
  ds = d[!d$from %in% duplicates & !d$to %in% duplicates,] ## duplicates that occur simultaneously
  ds$fromi = match(ds$from, rownames(dtm)) ## makes unique match to all ids in d
  ds$toi = match(ds$to, rownames(dtm))
  duplicates = unique(c(duplicates, 
                        as.character(ds$from[ds$fromi < ds$toi]),
                        as.character(ds$to[ds$fromi > ds$toi])))

  if (length(duplicates) == 0) {
    message("There are no duplicates")
    return(dtm)
  }
  
  message('Deleting ', length(duplicates), ' duplicates')
 
  if (!is.null(group.var)) {
    duplicates.med = meta[[group.var]][match(duplicates, rownames(dtm))]
    counts.med = table(duplicates.med)
    for(source in names(counts.med)){
        message('\t',source, ': ', counts.med[source])
    }
  } 
  
  d$is_duplicate = NA
  is_from = d$from %in% duplicates
  is_to = d$to %in% duplicates
  d$is_duplicate[is_from & !is_to] = 'from'
  d$is_duplicate[!is_from & is_to] = 'to'
  d$is_duplicate[is_from & is_to] = 'both'
  if (anyNA(d$is_duplicate)) warning(sprintf("There are %s document pairs for which neither document is marked as a duplicate. This shouldn't happen, so please report as a bug", sum(is.na(d$is_duplicate))))
  if (!is.null(dup_csv)) {
    d = d[order(d$from, d$hourdiff),]
    d$weight = round(d$weight,4)
    utils::write.csv(d, dup_csv, row.names=F)
  }

  dtm[!rownames(dtm) %in% duplicates,]
}

############
