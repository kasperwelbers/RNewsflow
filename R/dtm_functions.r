dtmToSparseMatrix <- function(dtm){
  sm = Matrix::spMatrix(nrow(dtm), ncol(dtm), dtm$i, dtm$j, dtm$v)
  rownames(sm) = rownames(dtm)
  colnames(sm) = colnames(dtm)
  sm
}

pad_dfm <- function(dfm, cnames) {
  ## not pretty, but quanteda:::pad_dfm is not exported
  dv = quanteda::docvars(dfm)
  rnames = rownames(dfm)
  dfm = methods::as(dfm, 'dgTMatrix')
  colindex = match(colnames(dfm), cnames)
  isna = is.na(colindex[dfm@j+1])
  out = Matrix::spMatrix(nrow(dfm), length(cnames), i=dfm@i[!isna]+1, j=colindex[dfm@j[!isna]+1], x = dfm@x[!isna])
  out = quanteda::as.dfm(out)
  rownames(out) = rnames
  colnames(out) = cnames
  if (nrow(dv) > 0) quanteda::docvars(out) = dv
  out
}

match.dtm.meta <- function(dtm, meta, id.var){
  if(!all(rownames(dtm) %in% meta[,id.var])) stop('Not all documents in DTM match with a document in the meta data.frame')
  meta[match(rownames(dtm), meta[,id.var]),]
}

#' Calculate statistics for term occurence across days
#'
#' @param dtm A quanteda \link[quanteda]{dfm}. Alternatively, a DocumentTermMatrix from the tm package can be used, but then the meta parameter needs to be specified manually
#' @param meta If dtm is a quanteda dfm, docvars(meta) is used by default (meta is NULL) to obtain the meta data. Otherwise, the meta data.frame has to be given by the user, with the rows of the meta data.frame matching the rows of the dtm (i.e. each row is a document)
#' @param date.var The name of the meta column specifying the document date. default is "date". The values should be of type POSIXlt or POSIXct
#'
#' @return A data.frame with statistics for each term.
#' \itemize{
#'  \item{freq:}{ The number of times a term occurred}
#'  \item{doc.freq:}{ The number of documents in which a term occured}
#'  \item{days.n:}{ The number of days on which a term occured}
#'  \item{days.pct:}{ The percentage of days on which a term occured}
#'  \item{days.entropy:}{ The entropy of the distribution of term frequency across days}
#'  \item{days.entropy.norm:}{ The normalized days.entropy, where 1 is a discrete uniform distribution}
#' }
#' @export
#'
#' @examples
#' tdd = term.day.dist(rnewsflow_dfm, date.var='date')
#' head(tdd)
#' tail(tdd)
term.day.dist <- function(dtm, meta=NULL, date.var='date'){
  if (is.null(meta)) meta = quanteda::docvars(dtm)
  dtm = quanteda::as.dfm(dtm)
  dtm = methods::as(dtm, 'dgTMatrix')
  if (!date.var %in% colnames(meta)) stop('The name specified in date.var is not a valid dfm docvar')
  
  document.date = as.Date(meta[[date.var]])
  
  if (any(is.na(document.date))) {
    message("date contains NA. These documents will be ignored")
    filter = !is.na(document.date)
    dtm = dtm[filter,]
    document.date = document.date[filter]
  }
  
  cs = Matrix::colSums(dtm)
  if(any(cs == 0)) {
    message("dtm contains empty columns/terms. These will be ignored (and won't appear in the output)")
    filter = Matrix::colSums(dtm) > 0
    dtm = dtm[,filter]
  } 
  
  
  
  dateseq = seq.Date(min(document.date), max(document.date), by='days')
  i = document.date[dtm@i+1]
  i = match(i, dateseq)
  m = Matrix::spMatrix(length(dateseq), ncol(dtm), i, dtm@j+1, dtm@x)
  
  m = methods::as(m, 'dgCMatrix')
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
  rownames(d) = NULL
  d
}

