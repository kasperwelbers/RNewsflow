#' Experimental: Convert dtm scores to a term innovation score, based on changes in term use over time
#'
#' For each term in m, the usage before and after the document date is compared (with a chi2 test) to see whether usage increased. 
#'
#' @param m A dgCMatrix 
#' @param date a character vector that specifies a date for each row in m. If given, only pairs of rows within a given date range (see lwindow, rwindow and date_unit) are calculated. 
#' @param m2 Optionally, use a different matrix for calculating the innovation scores. For example, if m is a DTM of press releases, m2 can be a DTM of news articles, to see if term usage increased in the news after the press release.
#' @param date2 If m2 is used, date2 has to be used to specify the date for the rows in m2 (otherwise date will be ignored)
#' @param lwindow If date (and date2) are used, lwindow determines the left side of the date window. e.g. -10 means that rows are only matched with rows for which date is within 10 [date_units] before.
#' @param rwindow Like lwindow, but for the right side. e.g. an lwindow of -1 and rwindow of 1, with date_unit is "days", means that only rows are matched for which the dates are within a 1 day distance
#' @param date_unit The date unit used in lwindow and rwindow. Supports "days", "hours", "minutes" and "seconds". Note that refers to the time distance between two rows ("days" doesn't refer to calendar days, but to a time of 24 hours)
#' @param min_chi  The minimum chi-square value
#' @param min_ratio The minimum ratio (rwindow score / lwindow score)
#' @param smooth    The smoothing factor (prevents -Inf/Inf ratio)
#'
#' @return A dgCMatrix
term_innovation <- function(m, date, m2=NULL, date2=NULL, lwindow=-7, rwindow=7, date_unit=c('days','hours','minutes','seconds'), min_chi=5.024, min_ratio=2, smooth=1) {
  date_unit = match.arg(date_unit)
  
  m = quanteda::as.dfm(m)
  if (is.null(m2)) {
    m2 = m
    date2 = date
  } else {
    m2 = quanteda::as.dfm(m2)
    if (!identical(quanteda::featnames(m), quanteda::featnames(m2))) 
      m2 = dfm_match(m2, quanteda::featnames(m))
  }
  if (!methods::is(date, 'POSIXct')) stop("Date has to be of type POSIXct (use as.POSIXct)")
  if (!methods::is(date2, 'POSIXct')) stop("Date has to be of type POSIXct (use as.POSIXct)")
  
  if (is.null(colnames(m))) colnames(m) = 1:ncol(m) ##     mainly for testing (normally there should really be column names)
  if (is.null(colnames(m2))) colnames(m2) = 1:ncol(m2) ## 

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
  
  ti = window_corp_comp(m,m2,order1,order2, lwindow, rwindow, min_chi=min_chi, min_ratio=min_ratio, smooth=smooth)
  dimnames(ti) = dimnames(m)
  
  ti = quanteda::as.dfm(ti)
  quanteda::docvars(ti) = quanteda::docvars(m)
  ti
}

