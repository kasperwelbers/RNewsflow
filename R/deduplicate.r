
#' Delete duplicate (or similar) documents from a document term matrix 
#' 
#' Delete duplicate (or similar) documents from a document term matrix. 
#' Duplicates are defined by: having high content similarity, occuring within a given time distance and being published by the same source.
#' 
#' Note that this can also be used to delete "updates" of articles (e.g., on news sites, news agencies). 
#' This should be considered if the temporal order of publications is relevant for the analysis. 
#' 
#' @param dtm         A quanteda \link[quanteda]{dfm}. 
#' @param date_var    The name of the column in docvars(dtm) that specifies the document date. The values should be of type POSIXlt or POSIXct
#' @param hour_window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param group_var   Optionally,  column name in docvars(dtm) that specifies a group (e.g., source, sourcetype). If given, only documents within the same group will be compared.
#' @param measure     The measure that should be used to calculate similarity/distance/adjacency. Currently supports the symmetrical measure "cosine" (cosine similarity), and the assymetrical measures "overlap_pct" (percentage of term scores in the document that also occur in the other document).
#' @param similarity  A threshold for similarity. Documents of which similarity is equal or higher are deleted
#' @param keep        A character indicating whether to keep the 'first' or 'last' published of duplicate documents.
#' @param tf_idf      If TRUE, weight the dtm with tf_idf before comparing documents. The original (non-weighted) DTM is returned.
#' @param dup_csv     Optionally, a path for writing a csv file with the duplicates edgelist. For each duplicate pair it is noted if "from" or "to" is the duplicate, or if "both" are duplicates (of other documents)
#' @param verbose     If TRUE, report progress
#' 
#' @return A dtm with the duplicate documents deleted
#' @export
#' 
#' @examples
#' ## example with very low similarity threshold (normally not recommended!)
#' dtm2 = delete_duplicates(rnewsflow_dfm, similarity = 0.5, keep='first', tf_idf = TRUE)
delete_duplicates <- function(dtm, date_var=NULL, hour_window=c(-24,24), group_var=NULL, measure=c('cosine','overlap_pct'), similarity=1, keep='first', tf_idf=FALSE, dup_csv=NULL, verbose=F){
  if (!methods::is(dtm, 'dfm')) stop('dtm has to be a quanteda dfm')
  measure = match.arg(measure)
  if(tf_idf) dtm = quanteda::dfm_tfidf(dtm)
  
  d = compare_documents(dtm, date_var=date_var, hour_window=hour_window, group_var=group_var, measure=measure, 
                             min_similarity=similarity, only_complete_window = F, verbose=verbose)
  d = d$d
  #e = igraph::get.edges(g, igraph::E(g))
  #d = igraph::get.data.frame(g, 'edges')  
  
  duplicates = c()
  if (is.null(date_var)) {
    if(keep == 'first') {
      duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff > 0])))
      duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff < 0])))
    }
    if(keep == 'last') {
      duplicates = c(duplicates, as.character(unique(d$from[d$hourdiff > 0])))
      duplicates = c(duplicates, as.character(unique(d$to[d$hourdiff < 0])))
    }
  } else {
    duplicates = unique(c(as.character(d$to), as.character(d$from)))
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
  
  if (!is.null(group_var)) {
    duplicates.med = quanteda::docvars(dtm)[[group_var]][match(duplicates, rownames(dtm))]
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
    d = d[order(d$from),]
    d$weight = round(d$weight,4)
    utils::write.csv(d, dup_csv, row.names=F)
  }
  
  dtm[!rownames(dtm) %in% duplicates,]
}