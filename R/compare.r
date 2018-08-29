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
#' @param dtm A quanteda \link[quanteda]{dfm}. Alternatively, a DocumentTermMatrix from the tm package can be used.
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
#' @param dtm A quanteda \link[quanteda]{dfm}. Alternatively, a DocumentTermMatrix from the tm package can be used, but then the meta parameter needs to be specified manually
#' @param meta If dtm is a quanteda dfm, docvars(meta) is used by default (meta is NULL) to obtain the meta data. Otherwise, the meta data.frame has to be given by the user, with the rows of the meta data.frame matching the rows of the dtm (i.e. each row is a document)
#' @param date.var The name of the column in meta that specifies the document date. default is "date". The values should be of type POSIXlt or POSIXct
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param group.var Optionally,  The name of the column in meta that specifies a group (e.g., source, sourcetype). If given, only documents within the same group will be compared.
#' @param measure the measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param min.similarity a threshold for similarity. lower values are deleted. Set to 0.1 by default.
#' @param n.topsim An alternative or additional sort of threshold for similarity. Only keep the [n.topsim] highest similarity scores for x. Can return more than [n.topsim] similarity scores in the case of duplicate similarities.
#' @param only.from A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare only these documents to other documents. 
#' @param only.to A vector with names/ids of documents (dtm rownames), or a logical vector that matches the rows of the dtm. Use to compare other documents to only these documents.
#' @param only.complete.window if True, only compare articles (x) of which a full window of reference articles (y) is available. Thus, for the first and last [window.size] days, there will be no results for x.
#' @param return_as Detemine whether output is returned as an "edgelist" or "igraph" network.
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
newsflow.compare <- function(dtm, meta=NULL, date.var='date', hour.window=c(-24,24), group.var=NULL, measure=c('cosine','overlap_pct'), 
                             min.similarity=0, n.topsim=NULL, only.from=NULL, only.to=NULL, only.complete.window=TRUE, return_as = c('igraph','edgelist'), verbose=FALSE){
  if (is.null(meta)) {
    if (!is(dtm, 'dfm')) stop('meta can only be NULL if dtm is a quanteda dfm class')
    meta = quanteda::docvars(dtm)
  }
  dtm = quanteda::as.dfm(dtm)
  meta$document_id = rownames(dtm)
  measure = match.arg(measure)
  return_as = match.arg(return_as)
  #meta = quanteda::docvars(dtm)
  
  #d = quanteda::as.dfm(dtm)
  #m = quanteda::docvars(d)
  #m$document_id = rownames(m)
  if (!date.var %in% colnames(meta)) stop('The name specified in date.var is not a valid dfm docvar')
  date = meta[[date.var]]
  if (!is(date, 'POSIXlt') && !is(date, 'POSIXct')) stop("Date has to be of type POSIXlt or POSIXct (use as.POSIXlt or strptime)")
  
  if (!is.null(group.var)) {
    if (!group.var %in% colnames(meta)) stop('The name specified in group.var is not a valid dfm docvar')
    group = meta[[group.var]]
  } else group = NULL
  
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
  
  if (only.complete.window) {
    left_out = date + as.difftime(hour.window[1], units = 'hours') < min(date)
    right_out = date + as.difftime(hour.window[2], units = 'hours') > max(date)
    dtm = dtm[!left_out & !right_out,]
    date = date[!left_out & !right_out]
    if (!is.null(group)) group = group[!left_out & !right_out]
  }
  
  diag = !is.null(dtm.y)
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, l2norm = T, min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2=group.y,
                                                  date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', verbose=verbose)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                       date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', verbose=verbose)
  cp = as(cp, 'dgTMatrix')
  
  if (return_as == 'edgelist') {
    if (!is.null(dtm.y)) 
      hourdiff = round(difftime(date.y[cp@j+1], date[cp@i+1], units = 'hours'),3)
    else
      hourdiff = round(difftime(date[cp@j+1], date[cp@i+1], units = 'hours'),3)
    cp = data.frame(from=rownames(cp)[cp@i+1], to=colnames(cp)[cp@j+1], weight=cp@x, hourdiff = hourdiff)
    return(cp[!as.character(cp$from) == as.character(cp$to),])
  } 
  
  if (return_as == 'igraph') {
    cp = data.frame(x=rownames(cp)[cp@i+1], y=colnames(cp)[cp@j+1], similarity=cp@x)
    cp = cp[!as.character(cp$x) == as.character(cp$y),]
  
    if (verbose) message('Creating network')
    g = document.network(cp, meta, 'document_id', date.var)
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
#' data(dtm)
#' 
#' ## example with very low similarity threshold (normally not recommended!)
#' dtm2 = delete.duplicates(dtm, similarity = 0.5, keep='first', tf.idf = TRUE)
delete.duplicates <- function(dtm, meta=NULL, date.var='date', hour.window=c(-24,24), group.var=NULL, measure=c('cosine','overlap_pct'), similarity=1, keep='first', tf.idf=FALSE, dup_csv=NULL, verbose=F){
  if (is.null(meta)) {
    if (!is(dtm, 'dfm')) stop('meta can only be NULL if dtm is a quanteda dfm class')
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
    write.csv(d, dup_csv, row.names=F)
  }

  dtm[!rownames(dtm) %in% duplicates,]
}

############

