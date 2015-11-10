library(tm)

#' Transform a dtm into a sparse matrix.
#' 
#' @param dtm a document-term matrix
#' @return a sparse matrix
dtmToSparseMatrix <- function(dtm){
  sm = spMatrix(nrow(dtm), ncol(dtm), dtm$i, dtm$j, dtm$v)
  rownames(sm) = rownames(dtm)
  colnames(sm) = colnames(dtm)
  sm
}

confirm.dtm.meta <- function(meta, id.var, date.var){
  if(!id.var %in% colnames(meta)) stop(sprintf('Meta data.frame should contain a column that matches the id.var parameter (currently set to "%s")', id.var))
  if(!date.var %in% colnames(meta)) stop(sprintf('Meta data.frame should contain a column that matches the id.var parameter (currently set to "%s")', date.var))
}

match.dtm.meta <- function(dtm, meta, id.var){
  if(mean(rownames(dtm) %in% meta[,id.var]) < 1) stop('Not all documents in DTM match with a document in the meta data.frame')
  meta[match(rownames(dtm), meta[,id.var]),]
}

term.day.dist <- function(dtm, meta, id.var='document_id', date.var='date'){
  confirm.dtm.meta(meta, id.var, date.var)
  meta = match.dtm.meta(dtm, meta, id.var)
  
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  
  cs = Matrix::colSums(dtm)
  if(sum(cs == 0) > 0) {
    message("dtm contains empty columns/terms. These will be ignored (and won't appear in the output)")
    dtm = dtm[,col_sums(dtm) > 0]
  } 
  
  document.date = as.Date(meta[,date.var])
  dateseq = seq.Date(min(document.date), max(document.date), by='days')
  i = document.date[dtm@i+1]
  i = match(i, dateseq)
  m = spMatrix(length(dateseq), ncol(dtm), i, dtm@j+1, dtm@x)
  
  m = as(m, 'dgCMatrix')
  days.entropy = columnEntropy(m)
  days.n = Matrix::colSums(m>0)

  d = data.frame(term=colnames(dtm),
                 freq = Matrix::colSums(dtm),
                 doc.freq = Matrix::colSums(dtm > 0),
                 days.n = days.n, 
                 days.pct = days.n / length(dateseq),
                 days.entropy = days.entropy, 
                 days.entropy.norm=days.entropy / length(dateseq))
  d$term = as.character(d$term)
  d
}
