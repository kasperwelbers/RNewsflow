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

termProduct <- function(m1, m2=m1){
  #m2@x[Matrix::which(m2@x > 0)] = 1
  Matrix::crossprod(m1,m2)
}

termOverlap_pct <- function(m1, m2=m1){
  totalterms = Matrix::colSums(as(m1, 'dgCMatrix'))
  m2@x[Matrix::which(m2@x > 0)] = 1
  Matrix::crossprod(m1,m2) / totalterms
}

termIndex <- function(m1, m2=m1){
  m2@x[Matrix::which(m2@x > 0)] = 1
  totalterms = Matrix::colSums(as(m1, 'dgCMatrix'))
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
  if(measure == 'overlap') results = termOverlap(m.x, m.y)
  if(measure == 'overlap_pct') results = termOverlap_pct(m.x, m.y)
  if(measure == 'product') results = termProduct(m.x, m.y)
  results
}

reindexTerms <- function(dtm, terms){
  dtm = dtmToSparseMatrix(dtm)
  documents = rownames(dtm)
  dtm = spMatrix(nrow(dtm), length(terms), dtm@i+1, match(colnames(dtm)[dtm@j+1], terms), dtm@x)
  dimnames(dtm) = list(documents, terms)
  as.DocumentTermMatrix(dtm, weighting = weightTf)
}

#' Compare the documents in two corpora/dtms
#' 
#' Compare the documents in corpus dtm.x with reference corpus dtm.y. 
#' 
#' @param dtm a document-term matrix in format of the tm package.
#' @param dtm.y Optional. If given, documents from dtm will only be compared to the documents in dtm.y
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently only cosine is supported
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0.1 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param return.zeros If true, all comparison results are returned, including those with zero similarity (quite possibly the worst thing to do with large data)
#' @return A data frame with sets of documents and their similarities. 
#' @export
documents.compare <- function(dtm, dtm.y=NULL, measure='cosine', min.similarity=0.1, n.topsim=NULL, return.zeros=F) {  
  if(!is.null(dtm.y)){
    if(mean(colnames(dtm) == colnames(dtm.y)) < 1){
      ## if colnames do not match, reindex them.
      terms = unique(c(colnames(dtm), colnames(dtm.y)))
      dtm = reindexTerms(dtm, terms)
      dtm.y = reindexTerms(dtm.y, terms)
    }
    m.x = Matrix::t(dtmToSparseMatrix(dtm))
    m.y = Matrix::t(dtmToSparseMatrix(dtm.y))
  } else {
    m.x = m.y = Matrix::t(dtmToSparseMatrix(dtm))
  }
    
  results = calculate.similarity(m.x, m.y, measure)
  results = filterResults(results, min.similarity, n.topsim)
  
  results = as(results, 'dgTMatrix')
  if(return.zeros) {
    results = Matrix(which(!is.na(results), arr.ind=T))
    results = data.frame(x=colnames(m.x)[results[,1]], y=colnames(m.y)[results[,2]], similarity=as.vector(results))
  } else{
    if(sum(results) == 0) return(NULL)
    results = data.frame(x=colnames(m.x)[results@i+1], y=colnames(m.y)[results@j+1], similarity=results@x)
    results = results[results$similarity > 0 & !is.na(results$similarity),]
  }
  results[!as.character(results$x) == as.character(results$y),]
}

#' @export
unlistWindow <- function(list_object, i, window){
  indices = i + window
  indices = indices[indices > 0 & indices <= length(list_object)]
  unlist(list_object[indices], use.names=F)
}

getDateIds <- function(date, row_filter=NULL){
  if(is.null(row_filter)) row_filter = rep(T, length(date))
  
  datetime = as.Date(date)
  datetimeseq = seq.Date(min(datetime), max(datetime), by='days')
  
  nonempty = which(datetimeseq %in% unique(datetime))
  nonempty_datetime_ids = llply(datetimeseq[nonempty], function(dtime) which(datetime == dtime & row_filter))
  datetime_ids = vector("list", length(datetimeseq))
  datetime_ids[nonempty] = nonempty_datetime_ids
  datetime_ids
}


#' Compare the documents in a dtm per time frame
#' 
#' This is an extension of the documents.compare function to only compare documents within a given range of days. 
#' 
#' Note that in this function it is not possible to use dtm.y to compare only certain documents (e.g., news and PR messages). Instead, the only.from and only.to parameters can be used for this purpose. 
#' 
#' @param dtm a document-term matrix in format of the tm package.
#' @param date a vector of date class, of the same length and order as the documents (rows) of the dtm.
#' @param window.size the timeframe in days within which articles must occur in order to be compared. e.g., if 0, articles are only compared to articles of the same day. If 1, articles are compared to all articles of the previous, same or next day.
#' @param window.direction For a more specific selection of which articles in the window to compare to. This is given with a combination of the symbols '<' (before x) '=' (simultanous with x) and '>' (after x). default is '<=>', which means all articles. '<>' means all articles before or after the [time.unit] of an article itself. '<' means all previous articles, and '<=' means all previous and simultaneous articles. etc.  
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently only cosine is supported
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0.1 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param only.from A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare only these documents to other documents. 
#' @param only.to A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare other documents to only these documents.
#' @param return.zeros If true, all comparison results are returned, including those with zero similarity (quite possibly the worst thing to do with large data)
#' @param return.date If true, the dates for x and y are given in the output
#' @param only.complete.window if True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' @return A data frame with columns x, y and similarity. If return.date == T, date.x and date.y are returned as well.
#' @export
documents.window.compare <- function(dtm, meta, id.var='document_id', date.var='date', hour.window=c(-24,24), measure='cosine', min.similarity=0, n.topsim=NULL, only.from=NULL, only.to=NULL, return.zeros=F, only.complete.window=T){
  confirm.dtm.meta(meta, id.var, date.var)
  meta = match.dtm.meta(dtm, meta, id.var)
  
  message('Indexing articles by date/time')
  if(is.null(only.from) & is.null(only.to)){
    dateids.x = dateids.y = getDateIds(meta[,date.var])
  } else{ 
    if(is.null(only.from)) only.from = rep(T, nrow(dtm))
    if(is.null(only.to)) only.to = rep(T, nrow(dtm))
    if(!class(only.from) == 'logical') only.from = rownames(dtm) %in% only.from
    if(!class(only.to) == 'logical') only.to = rownames(dtm) %in% only.to
    dateids.x = getDateIds(meta[,date.var], only.from)
    dateids.y = getDateIds(meta[,date.var], only.to)
  }
  dateindex = which(lapply(dateids.x, length) > 0)
  
  window = floor(hour.window[1]/24):ceiling(hour.window[2]/24)
  if(only.complete.window){
    if(window[1] < 0) dateindex = dateindex[dateindex > window[1]]
    if(rev(window)[1] > 0) dateindex = dateindex[dateindex <= length(dateids.x) - rev(window)[1]]  
  }
  
  message('Comparing documents')
  output = ldply(dateindex, function(i) ldply_documents.compare(i, dtm, dateids.x, dateids.y, window, measure, min.similarity, n.topsim, return.zeros), .progress='text')
  output = output[,!colnames(output) == '.id']
  
  message('Matching document meta')
  g = document.network(output, meta, id.var, date.var)
  
  delete.pairs = which(E(g)$hourdiff < hour.window[1] | E(g)$hourdiff > hour.window[2])
  g = delete.edges(g, delete.pairs)
  g
}

ldply_documents.compare <- function(i, dtm, dateids.x, dateids.y, window, measure, min.similarity, n.topsim, return.zeros){
  ## special function to be used in ldply in document.window.compare
  dtm.x_indices = unique(dateids.x[[i]])
  dtm.y_indices = unique(unlistWindow(dateids.y,i,window))
  if(length(dtm.y_indices) == 0) return(NULL)

  documents.compare(dtm[dtm.x_indices,], dtm[dtm.y_indices,], measure, min.similarity, n.topsim, return.zeros)
}

#' Delete duplicate (or similar) documents from a document term matrix 
#' 
#' Delete duplicate (or similar) documents from a document term matrix 
#' 
#' @param dtm a document-term matrix in format of the tm package.
#' @param date a vector of date class, of the same length and order as the documents (rows) of the dtm.
#' @param source a vector of the same length and order as the documents (rows) of the dtm, indicating what the source of the document is.
#' @param window.size the timeframe in days within which articles must occur in order to be compared. e.g., if 0, articles are only compared to articles of the same day. If 1, articles are compared to all articles of the previous, same or next day.
#' @param similarity a threshold for similarity. lower values are deleted
#' @param keep A character indicating whether to keep the 'first' or 'last' of duplicate documents.
#' @param tf.idf if True, weight the dtm with tf.idf before comparing documents
#' @return A dtm with the duplicate documents deleted
#' @export
delete.duplicates <- function(dtm, meta, id.var='document_id', date.var='date', source.var='source', hour.window=c(-24,24), similarity=1, keep='first', tf.idf=F){
  if(tf.idf) {
    g = documents.window.compare(weightTfIdf(dtm), meta, min.similarity = similarity, hour.window=hour.window)
  } else {
    g = documents.window.compare(dtm, meta, min.similarity = similarity, hour.window=hour.window)
  }
  
  e = get.edges(g, E(g))
  d = get.data.frame(g, 'edges')  
  d$med.x = V(g)$source[e[,1]]
  d$med.y = V(g)$source[e[,2]]
  d = d[d$med.x == d$med.y,]
  
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

