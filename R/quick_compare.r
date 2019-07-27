match_by <- function(x_text, y_text, stem=T, remove=NULL, weight=T, ...) {
  c(list(x_text=x_text, y_text=y_text, stem=stem, remove=remove, weight=weight), 
    list(...))
}

weight_dtm_x <- function(dfm_x, dfm_y) {
  ts = colSums(dfm_y > 0)
  ts = ts[match(colnames(dfm_x), names(ts))]
  ts[is.na(ts)] = 0
  idf = log(1 + (nrow(dfm_y) / (ts+1)))
  t(t(dfm_x) * idf)
}

match_texts <- function(x, y, date.var='date', date.var.x=date.var, date.var.y=date.var, ...,
                         hour.window=c(-24,24), measure='overlap_pct', min.similarity=0, verbose=T) {
  mb = list(...)
  d = NULL
  out = list()
  y[[date.var]] = as.POSIXct(y[[date.var.x]])
  x[[date.var]] = as.POSIXct(x[[date.var.y]])
  for (i in seq_along(mb)) {
    label = names(mb)[i]
    message(paste('Matching', label))
    
    dfm_args = mb[-(1:2)]

    if (!identical(d[-2], mb[[i]][-2])) {
      text = y[[mb[[i]]$y_text]]
      dfm_y = do.call(quanteda::dfm, args = c(list(x=as.character(text)), dfm_args))
    }
    text = x[[mb[[i]]$x_text]]
    dfm_x = do.call(quanteda::dfm, args = c(list(x=as.character(text)), dfm_args))
    
    if (mb[[i]]$weight) dfm_x = weight_dtm_x(dfm_x, dfm_y)
    if (mb)
    rownames(dfm_x) = 1:nrow(dfm_x)
    rownames(dfm_y) = 1:nrow(dfm_y)
    out[[label]] = RNewsflow::newsflow.compare(dfm_x, dfm_y, meta = x, meta.y = y, date.var = date.var, only.complete.window = F, return_as = 'edgelist', 
                                           measure = mb[[i]]$measure, hour.window = mb[[i]]$hour.window, min.similarity = min.similarity, verbose=verbose)
  }
  out = data.table::rbindlist(out, id='match_by')
  data.table::setnames(out, old=c('from','to'), c('x','y'))
  data.table::dcast(out, x + y + hourdiff ~ match_by, value.var='weight')
}
