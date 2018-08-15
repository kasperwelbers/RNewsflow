## FUNCTIONS FOR COMPARING DOCUMENTS


cosineSimilarity <- function(m1, m2=NULL){
  norm = sqrt(Matrix::colSums(m1^2))
  m1@x = m1@x / norm[m1@j+1]  
  if(!is.null(m2)){
    norm = sqrt(Matrix::colSums(m2^2))
    m2@x = m2@x / norm[m2@j+1]
    cp = Matrix::crossprod(m1,m2) 
  } else cp = Matrix::crossprod(m1)
  cp
}

termOverlap <- function(m1, m2=m1){
  m2@x[Matrix::which(m2@x > 0)] = 1
  Matrix::crossprod(m1,m2)
}

termOverlap_pct <- function(m1, m2=m1, reverse=FALSE){
  totalterms = if(!reverse) Matrix::colSums(methods::as(m1, 'dgCMatrix')) else Matrix::colSums(methods::as(m2, 'dgCMatrix'))
  m2@x[Matrix::which(m2@x > 0)] = 1
  Matrix::crossprod(m1,m2) / totalterms
}

termIndex <- function(m1, m2=m1){
  m2@x[Matrix::which(m2@x > 0)] = 1
  totalterms = Matrix::colSums(methods::as(m1, 'dgCMatrix'))
  Matrix::crossprod(m1,m2) / totalterms
}

Nth.max <- function(x, N){
  N = min(N, length(x)) 
  -sort(-x, partial=N)[N]
}

filterResults <- function(results, min.similarity, n.topsim){
  if(!is.null(min.similarity)) results@x[Matrix::which(results@x < min.similarity)] = 0
  if(!is.null(n.topsim)) {
    simthres = apply(results, 1, Nth.max, N=n.topsim)
    results@x[Matrix::which(results < simthres)] = 0
  }
  results
}
  
calculate.similarity <- function(m.x, m.y, measure){
  if(measure == 'cosine') results = cosineSimilarity(m.x, m.y)
  if(measure == 'percentage.from') results = termOverlap_pct(m.x, m.y)
  if(measure == 'percentage.to') results = termOverlap_pct(m.x, m.y, reverse = TRUE)
  results
}


reindexTerms <- function(dtm, terms){
  dtm = as(dtm, 'dgTMatrix')
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
#' @param dtm A document-term matrix in the tm \link[tm]{DocumentTermMatrix} class. It is recommended to weight the DTM beforehand, for instance using \link[tm]{weightTfIdf}.
#' @param dtm.y Optional. If given, documents from dtm will only be compared to the documents in dtm.y
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine", for cosine similarity. Also supports assymetrical measures "percentage.from" and "percentage.to" for the percentage of overlapping terms (term scores taken into account). Here "percentage.from" gives the percentage of the document that is compared to the other, whereas "percentage.to" gives the percentage of the document to which is compared.
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' 
#' @return A data frame with pairs of documents and their similarities. 
#' @export
#' 
#' @import tm
#' 
#' @examples
#' data(dtm)
#' 
#' comp = documents.compare(dtm, min.similarity=0.4)
#' head(comp)
documents.compare <- function(dtm, dtm.y=NULL, measure=c('cosine','overlap_pct'), min.similarity=0, n.topsim=NULL) {  
  measure = match.arg(measure)
  
  out = vector('list')

  if(!is.null(dtm.y)){
    if(!all(colnames(dtm) == colnames(dtm.y))){
      ## if colnames do not match, reindex them.
      terms = unique(c(colnames(dtm), colnames(dtm.y)))
      dtm = reindexTerms(dtm, terms)
      dtm.y = reindexTerms(dtm.y, terms)
    }
  }
  dtm = as(dtm, 'dgCMatrix')
  if (!is.null(dtm.y)) dtm.y = as(dtm.y, 'dgCMatrix')
  
  diag = !is.null(dtm.y)
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, l2norm = T, min_value = min.similarity, top_n = n.topsim, diag=diag)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag)
  cp = as(cp, 'dgTMatrix')
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
#' @param dtm A quanteda \link[quanteda{dfm}]
#' @param date.var The name of the dfm docvar containing the document date. default is "date". The values should be of type POSIXlt or POSIXct
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param group.var Optionally, the name of a dfm docvar to specify groups (e.g., document of the same source). If given, only documents within the same group will be compared.
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0.1 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param only.from A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare only these documents to other documents. 
#' @param only.to A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare other documents to only these documents.
#' @param only.complete.window if True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' @param verbose If TRUE, report progress
#' 
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#' 
#' @examples 
#' data(dtm)
#' data(meta)
#' 
#' dtm = tm::weightTfIdf(dtm)
#' g = newsflow.compare(dtm, meta, hour.window = c(0.1, 36))
#' 
#' vcount(g) # number of documents, or vertices
#' ecount(g) # number of document pairs, or edges
#' 
#' head(igraph::get.data.frame(g, 'vertices'))
#' head(igraph::get.data.frame(g, 'edges'))
newsflow.compare <- function(dtm, date.var='date', hour.window=c(-24,24), group.var=NULL, measure=c('cosine','overlap_pct'), 
                             min.similarity=0, n.topsim=NULL, only.from=NULL, only.to=NULL, only.complete.window=TRUE, verbose=FALSE){
  measure = match.arg(measure)
  meta = quanteda::docvars(dtm)

  if (!date.var %in% colnames(meta)) stop('The name specified in date.var is not a valid dfm docvar')
  date = quanteda::docvars(dtm, date.var)
  if (!is(date, 'POSIXlt') && !is(date, 'POSIXct')) stop("Date has to be of type POSIXlt or POSIXct (use as.POSIXlt or strptime)")
  
  if (!is.null(group.var)) {
    if (!group.var %in% colnames(quanteda::docvars(dtm))) stop('The name specified in group.var is not a valid dfm docvar')
    group = quanteda::docvars(dtm, group.var)
  } else group = NULL
  
  if (is.null(only.from) && is.null(only.to)) {
    dtm.y = NULL
    date.y = NULL
    group.y = NULL
  }
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
    if (!is.null(group)) group.y = group.y[only.from]
  }
  
  if (only.complete.window) {
    left_out = date + as.difftime(hour.window[1], units = 'hours') < min(date)
    right_out = date + as.difftime(hour.window[2], units = 'hours') > max(date)
    dtm = dtm[!left_out & !right_out,]
    date = date[!left_out & !right_out]
  }
  
  diag = !is.null(dtm.y)
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, l2norm = T, min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2=group.y,
                                                  date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', verbose=verbose)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                       date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', verbose=verbose)
  cp = as(cp, 'dgTMatrix')
  cp = data.frame(x=rownames(cp)[cp@i+1], y=colnames(cp)[cp@j+1], similarity=cp@x)
  cp = cp[!as.character(cp$x) == as.character(cp$y),]
  
  if (verbose) message('Matching document meta')
  meta$document_id = rownames(meta)
  rownames(meta) = NULL
  g = document.network(cp, meta, 'document_id', date.var)
}

#' Delete duplicate (or similar) documents from a document term matrix 
#' 
#' Delete duplicate (or similar) documents from a document term matrix. 
#' Duplicates are defined by: having high content similarity, occuring within a given time distance and being published by the same source.
#' 
#' Note that this can also be used to delete "updates" of articles (e.g., on news sites, news agencies). 
#' This should be considered if the temporal order of publications is relevant for the analysis. 
#' 
#' @param dtm A document-term matrix in the tm \link[tm]{DocumentTermMatrix} class. It is recommended to weight the DTM beforehand, for instance using \link[tm]{weightTfIdf}.
#' @param date.var Optionally, the name of the dfm docvar containing the document date. If given, only documents within the given hour.window will be compared. default is "date". The values should be of type POSIXlt or POSIXct
#' @param group.var Optionally, the name of the dfm docvar containing a group of documents (e.g. of the same source). If given, only documents of the same group are compared.
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. By default c(-24,24), which compares each document to all other documents within a 24 hour time distance.
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param similarity a threshold for similarity. Documents of which similarity is equal or higher are deleted
#' @param keep A character indicating whether to keep the 'first' or 'last' published of duplicate documents.
#' @param tf.idf if TRUE, weight the dtm with tf.idf before comparing documents. The original (non-weighted) DTM is returned.
#' 
#' @return A dtm with the duplicate documents deleted
#' @export
#' 
#' @examples
#' data(dtm)
#' 
#' ## example with very low similarity threshold (normally not recommended!)
#' dtm2 = delete.duplicates(dtm, similarity = 0.5, keep='first', tf.idf = TRUE)
delete.duplicates <- function(dtm, date.var='date', group.var=NULL, hour.window=c(-24,24), measure=c('cosine','overlap_pct'), similarity=1, keep='first', tf.idf=FALSE){
  if(tf.idf) dtm = quanteda::dfm_tfidf(dtm)
  g = newsflow.compare(dtm, date.var, hour.window, group.var, measure=measure, min.similarity=similarity, only.complete.window = F)
  
  e = igraph::get.edges(g, igraph::E(g))
  d = igraph::get.data.frame(g, 'edges')  
  
  duplicates = c()
  if(keep == 'first') {
    duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff > 0])))
    duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff < 0])))
  }
  if(keep == 'last') {
    duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff > 0])))
    duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff < 0])))
  }
  d = d[!d$from %in% duplicates & !d$to %in% duplicates,]
  
  ## if there are identical articles that occured simultaneously, delete randomly
  d = d[sample(1:nrow(d), nrow(d)),]
  d$fromi = match(d$from, unique(d$from, d$to))
  d$toi = match(d$to, unique(d$from, d$to))
  d = d[d$fromi < d$toi,]
  duplicates = unique(c(duplicates, as.character(d$from)))
  
  message('Deleting ', length(duplicates), ' duplicates')
 
  duplicates.med = meta$source[match(duplicates, rownames(dtm))]
  counts.med = table(duplicates.med)
  for(source in names(counts.med)){
      message('\t',source, ': ', counts.med[source])
  }
  
  dtm[!rownames(dtm) %in% duplicates,]
}

############

