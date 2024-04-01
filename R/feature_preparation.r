#' Automatically infer queries from combinations of terms in a dtm
#'
#' This function was designed for the task of matching short event descriptions to news articles, but can more generally
#' be used for document matching tasks. However, it should be noted that it will require exponentially more
#' memory for dtms with more unique terms, which is why it is less suitable for matching larger documents. This only applies to
#' the dtm, not the ref_dtm. Thus, if your goal is to match smaller documents such as event descriptions to news, this function
#' might be usefull. 
#'
#' The main purpose of the function is that it intersects the terms in a dtm based to increase sparsity. 
#' This can improve certain document matching tasks, but at the cost of creating a bigger dtm. 
#' If all terms are combined this would be a quadratic increase of columns.
#' However, only term combinations that occur in dtm (not ref_dtm) will be used.
#' This is not a problem as long as the similarity of the documents in dtm to documents in dtm_y is calculated
#' as an assymetric similarity measure (i.e. in which the sum of terms in dtm_y is not used). 
#' 
#' To emphasize that this feature preparation step is geared towards the task of 'looking up' documents,
#' we use the terminolog of a 'query'. The output of the function is a list of two dtm: query_dtm and ref_dtm.
#' Both dtms have the exact same columns that contain the query terms.
#' The values in query_dtm are by default tfidf weighted, and the values in ref_dtm are binary.
#'  
#' Several options are given to only create term combinations that are informative. Firstly, a minimum and maximum document frequency of term combinations can be defined. 
#' Secondly, a minimum observed/expected ratio can be given. The expected probability of a combination of term A and term B
#' is the joint probability. If the observed probability is not higher, the combination is not more informative than chance.
#' Thirdly, before intersecting terms, one can first cluster very similar terms together as single columns to reduce the number
#' of possible combinations. 
#'
#' @param dtm          A quanteda \link[quanteda]{dfm}
#' @param ref_dtm      Optionally, another quanteda \link[quanteda]{dfm}. If given, the ref_dtm will be used to calculate the docfreq/docprob scores.
#' @param min_docfreq  The minimum frequency for terms or combinations of terms
#' @param max_docprob  The maximum probability (document frequency / N) for terms or combinations of terms
#' @param weight       Determine how to weight the queries (if ref_dtm is used, uses the idf of the ref_dtm, or of both the dtm and ref dtm if use_dtm_and_ref is T). 
#'                     Default is "binary" (does/does not occur). "tfidf" uses common tf-idf weighting (actually just idf, since scores are binary). 
#' @param norm_weight  Normalize the weight score so that the highest value is 1. If "max" is used, max is the highest possible value. "doc_max" uses the highest value within each document, and "dtm_max" uses the highest observed value in the dtm.
#' @param min_obs_exp  The minimum ratio of the observed and expected frequency of a term combination
#' @param union_sim_thres If given, a number between 0 and 1, used as the cosine similarity threshold for combining clusters of terms 
#' @param combine_all  If True, combine all terms. If False (default), terms that are included as unigrams (i.e. that are within the min_docfreq and max_docprob) are not combined with other terms.
#' @param only_dtm_combs Only include term combinations that occur in dtm. This makes sense (and saves a lot of memory) if you are only interested in assymetric similarity measures based on the query
#' @param use_dtm_and_ref if a ref_dtm is used, the weight is computed based only on the document frequencies in the ref dtm. If use_dtm_and_ref is set to TRUE, both the dtm and ref_dtm are used.
#' @param verbose      If true, report progress
#'
#' @return a list with a query dtm and ref_dtm. Designed for use in \code{\link{compare_documents}} using the special `query_lookup` measure
#' @export
#'
#' @examples
#'  q = create_queries(rnewsflow_dfm, min_docfreq = 2, union_sim_thres = 0.9, 
#'                     max_docprob = 0.05, verbose = FALSE)
#'  head(colnames(q$query_dtm),100)
create_queries <- function(dtm, ref_dtm=NULL, min_docfreq=2, max_docprob=0.01, weight=c('tfidf', 'binary'), norm_weight=c('max','doc_max','dtm_max','none'), min_obs_exp=NA, union_sim_thres=NA, combine_all=T, only_dtm_combs=T, use_dtm_and_ref=F, verbose=F) {
  weight = match.arg(weight)
  norm_weight = match.arg(norm_weight)
  if (!methods::is(dtm, 'dfm')) stop('dtm has to be a quanteda dfm')
  if (!is.null(ref_dtm) && !methods::is(dtm, 'dfm')) stop('ref_dtm has to be a quanteda dfm')
  
  if (!is.null(ref_dtm)) {
    voc = colnames(ref_dtm)[Matrix::colSums(ref_dtm) >= min_docfreq]
    voc = intersect(colnames(dtm), voc)
    m = ref_dtm[,voc]
    if (use_dtm_and_ref) {
      m = rbind(m, dtm[,voc])                                  ## add dtm rows
      orig_ref = c(rep(T, nrow(ref_dtm)), rep(F, nrow(dtm)))   ## remember which, to remove later
    }
  } else {
    m = dtm[,Matrix::colSums(dtm > 0) > min_docfreq]
  }
  if (min_docfreq < 0) stop('min docfreq must be zero or higher')
  if (max_docprob <= 0 || max_docprob > 1) stop('max_docprob must be a value between 0 and 1')
  
  if (!is.na(union_sim_thres)) {
    if (union_sim_thres <= 0 || union_sim_thres > 1) stop('Union sim thres must be a value between 0 and 1')
    if (verbose) message('Computing clusters of similar terms')
    simmat1 = term_occur_sim(m, union_sim_thres, verbose=verbose)
    m = term_union(m, methods::as(simmat1, 'CsparseMatrix'), as_dfm=F)
  }
  
  if (verbose) message('Computing term combinations')
  simmat2 = term_cooccurence_docprob(m, max_docfreq = max_docprob * nrow(m), min_obs_exp=min_obs_exp,
                                    min_docfreq=min_docfreq)


  if (verbose) message('Building new dtm')
  if (!is.null(ref_dtm)) {
    m_ref = m ## remember dtm_ref
    
    m = dtm[,voc]
    if (!is.na(union_sim_thres)) {
      m = term_union(m, methods::as(simmat1,'CsparseMatrix'), as_dfm = F, verbose = F)
      voc = intersect(colnames(m), colnames(simmat2))
      m = m[,voc]
      simmat2 = simmat2[voc,voc]
    }
    
    if (only_dtm_combs) {
      simmat_filter = term_cooccurence_docprob(m, max_docfreq = nrow(m), min_obs_exp=min_obs_exp,
                                       min_docfreq=1)
      nz = which(simmat2 > 0)
      dropval = nz[simmat_filter[nz] == 0]
      simmat2@x[dropval] = 0
      simmat2 = Matrix::drop0(simmat2)
    }
    
    if (!combine_all) simmat2 = rm_comb_if_diag(simmat2)
    m = term_intersect(methods::as(m, 'dMatrix'), methods::as(simmat2, 'dMatrix'), as_dfm=F, verbose=F)
    
    if (verbose) message('Building new reference dtm')
    m_ref = m_ref[,voc]
    m_ref = term_intersect(methods::as(m_ref, 'dMatrix'), methods::as(simmat2, 'dMatrix'), as_dfm=F, verbose=F)
    
  } else {
    if (!combine_all) simmat2 = rm_comb_if_diag(simmat2)
    m = term_intersect(methods::as(m, 'dMatrix'), methods::as(simmat2, 'dMatrix'), as_dfm=F)
    m_ref = NULL
  }
   
  if (!ncol(m) == 0){
    m = weight_queries(m, m_ref, weight, norm_weight)
    if (!is.null(m_ref)) m_ref = weight_queries(m_ref, m_ref, weight, norm_weight)
  }

  m = quanteda::as.dfm(m)
  quanteda::docvars(m) = quanteda::docvars(dtm)
  
  if (!is.null(m_ref)) {
    if (use_dtm_and_ref) m_ref = m_ref[orig_ref,]   ## remove dtm rows
    m_ref = quanteda::as.dfm(m_ref)
    quanteda::docvars(m_ref) = quanteda::docvars(ref_dtm)
  } else {
    m_ref = m
  }
  
  l = list(query_dtm=m, ref_dtm=m_ref)
  class(l) = c('RNewflow_queries', class(l))
  l
}

weight_queries <- function(dfm_x, dfm_y=NULL, weight, norm_weight) {
  if (weight == 'binary') return(dfm_x > 0)
  if (is.null(dfm_y)) dfm_y = dfm_x
  ts = Matrix::colSums(dfm_y > 0)  
  ts = ts[match(colnames(dfm_x), names(ts))]
  ts[is.na(ts)] = 0
  
  if (weight == 'docprob') w = (ts / nrow(dfm_y))
  if (weight == 'tfidf') w = log((nrow(dfm_y) - ts + 0.5) / (ts+0.5), base=2)
  if (weight == 'tfidf_sq') w = log((nrow(dfm_y) - ts + 0.5) / (ts+0.5), base=2)^2
  
  dfm_x = t(t(dfm_x > 0) * w)
  
  if (norm_weight == 'doc_max') dfm_x = dfm_x / apply(dfm_x, MARGIN = 1, max)
  if (norm_weight == 'dtm_max') dfm_x = dfm_x / max(dfm_x)
  if (norm_weight == 'max') {
    if (weight == 'docprob') max_possible = (1 / nrow(dfm_y))
    if (weight == 'tfidf') max_possible = log((nrow(dfm_y) - 1 + 0.5) / (1+0.5), base=2)
    if (weight == 'tfidf_sq') max_possible = log((nrow(dfm_y) - 1 + 0.5) / (1+0.5), base=2)^2
    dfm_x = dfm_x / max_possible
  }
  return(dfm_x)
}


match_simmat_terms <- function(dtm, simmat) {
  if (!all(colnames(dtm) %in% colnames(simmat))) stop('not all terms in dtm are in  simmat')
  if (!ncol(dtm) == ncol(simmat)) {
    simmat = simmat[colnames(dtm), colnames(dtm)]
  } else {
    if (!all(colnames(dtm) == colnames(simmat))) {
      simmat = simmat[colnames(dtm), colnames(dtm)]
    }
  }
  simmat
}

#' Combine terms in a dtm
#' 
#' Given a dtm and a similarity (adjacency) matrix, group clusters of similar terms (simmat > 0) into a single column.
#' Column names will be concatenated, with a "|" seperator (read as OR)
#'
#' @param dtm          A quanteda \link[quanteda]{dfm} or a CsparseMatrix.
#' @param simmat       A similarity matrix in CsparseMatrix format. For instance, created with \link{term_char_sim}
#' @param as_dfm       If True, return as quanteda dfm
#' @param verbose      If True, report progress
#' @param sep          The separator used for pasting the terms
#' @param par          If TRUE, add parentheses to colnames before combining. This is mainly for internal use, as it allows
#'                     specification if OR (term_union) and AND (term_intersect) operations are combined. 
#'                     If NA, this is based on whether parenthese are present. 
#'
#' @return  A CsparseMatrix or quanteda dfm
#' @export
#' 
#' @examples 
#' dfm = quanteda::tokens(c('That guy Gadaffi','Do you mean Kadaffi?',
#'                          'Nah more like Gadaffel','Not Kadaffel?')) |>
#'   quanteda::dfm()
#' simmat = term_char_sim(colnames(dfm), same_start=0)
#' term_union(dfm, simmat, verbose = FALSE)
term_union <- function(dtm, simmat, as_dfm=T, verbose=F, sep='|', par=NA) {
  if (methods::is(dtm, "DocumentTermMatrix")) stop('this function does not work for tm DocumentTermMatrix class')
  dtm = quanteda::dfm_match(quanteda::as.dfm(dtm), colnames(simmat))
  
  parentheses = if (is.na(par)) grepl('[&]', colnames(dtm)) else par
  
  ml = term_union_cpp(methods::as(methods::as(methods::as(dtm, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
                      methods::as(methods::as(methods::as(simmat, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
                      colnames(dtm), parentheses, verbose, sep)
  colnames(ml$m) = ml$colnames
  rownames(ml$m) = rownames(dtm)
  ml$m = ml$m[,colSums(ml$m) > 0]
  
  if (as_dfm && methods::is(dtm, 'dfm')) {
    m = quanteda::as.dfm(ml$m > 0)
    quanteda::docvars(m) = quanteda::docvars(dtm)
    return(m)
  } else {
    return(ml$m > 0)
  }
}

#' Combine terms in a dtm
#' 
#' Given a dtm and a similarity (adjacency) matrix, create a new column for each nonzero cell in the
#' similarity matrix. For the term combinations  (everything except the diagonal) the column names will be
#' pasted together with a "&" separator (read as AND)
#'
#' @param dtm          A quanteda \link[quanteda]{dfm} or a CsparseMatrix.
#' @param simmat       A similarity matrix in CsparseMatrix format. For instance, created with \link{term_char_sim}
#' @param as_dfm       If True, return as quanteda dfm
#' @param verbose      If True, report progress
#' @param sep          The separator used for pasting the terms
#' @param par          If TRUE, add parentheses to colnames before combining. This is mainly for internal use, as it allows
#'                     specification if OR (term_union) and AND (term_intersect) operations are combined. 
#'                     If NA, this is based on whether parenthese are present.
#'
#' @return  A CsparseMatrix or quanteda dfm
#' @export
term_intersect <- function(dtm, simmat, as_dfm=T, verbose=F, sep=' & ', par=NA) {
  if (methods::is(dtm, "DocumentTermMatrix")) stop('this function does not work for tm DocumentTermMatrix class')
  #simmat = match_simmat_terms(dtm, simmat)
  dtm <- quanteda::dfm_match(quanteda::as.dfm(dtm), colnames(simmat))

  parentheses = if (is.na(par)) grepl('[|]', colnames(dtm)) else par
  ml = term_intersect_cpp(methods::as(methods::as(methods::as(dtm, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
                          methods::as(methods::as(methods::as(simmat, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
                          colnames(dtm), parentheses, verbose, sep)
  colnames(ml$m) = ml$colnames
  rownames(ml$m) = rownames(dtm)
  ml$m = ml$m[,Matrix::colSums(ml$m) > 0]
  
  if (as_dfm && methods::is(dtm, 'dfm')) {
    m = quanteda::as.dfm(ml$m > 0)
    quanteda::docvars(m) = quanteda::docvars(dtm)
    return(dtm)
  } else {
    return(ml$m > 0)
  }
}

term_occur_sim <- function(m, min_cos, verbose=F) {
  simmat = tcrossprod_sparse(t(m), min_value = min_cos, normalize='l2', verbose=verbose) > 0
  methods::as(simmat, 'CsparseMatrix')
}

## create a matrix with document probabilities (docfreq / n) for all column combinations
## max_docfreq is used to only keep sufficiently rare combinations
## typically, min_docfreq is used to drop very sparse terms, and max_docfreq is used to drop terms that are too common to be informative.
term_cooccurence_docprob <- function(m, max_docfreq, min_docfreq=NULL, min_obs_exp=NA) {
  #simmat = tcrossprod_sparse(methods::as(t(m > 0), 'CsparseMatrix'), min_value=min_docfreq, max_value = max_docfreq, verbose=verbose)
  simmat = crossprod(m>0)
  simmat@x[simmat@x < min_docfreq] = 0
  simmat@x[simmat@x > max_docfreq] = 0
  simmat = Matrix::drop0(simmat)
  
  if (!is.na(min_obs_exp)) {
    simmat = methods::as(methods::as(simmat, 'generalMatrix'), 'TsparseMatrix')
    prob = Matrix::colMeans(m > 0)
    exp = prob[simmat@i+1] * prob[simmat@j+1] * nrow(m)
    simmat@x[(simmat@x / exp) < min_obs_exp] = 0
    simmat = Matrix::drop0(simmat)
  }
  
  methods::as(simmat, 'CsparseMatrix')
}

## if the docfreq of a term is lower than max_docfreq
## its combinations with other terms are ignored. This allows us to only include combinations
## for terms that are not informative enough on their own
rm_comb_if_diag <- function(simmat) {
  d = diag(simmat)
  has_diag = d > 0
  if (any(has_diag)) {
    simmat[,has_diag] = 0
    simmat[has_diag,] = 0
    simmat = Matrix::drop0(simmat)
    diag(simmat) = d
  }
  simmat
}

char_grams <- function(x, type=c('tri','bi'), pad=T, drop_non_alpha=T, min_length=3) {
  type = match.arg(type)
  voc = x
  if (!pad && type == 'tri' && min_length < 3) stop('cannot use trigrams if length < 3 and pad = F')
  if (!pad && type == 'bi' && min_length < 2) stop('cannot use bigrams if length < 3 and pad = F')
  if (drop_non_alpha) x[!grepl('[a-zA-Z]', x)] = ''
  uni = stringi::stri_split_boundaries(x, type='character')
  uni = lapply(uni, FUN=function(x) {
    if (length(x) < min_length) return(c())
    if (pad) {
      if (type == 'bi') return(paste(c('lb', x), c(x, 'rb'), sep='_'))
      if (type == 'tri') return(paste(c('lb','lb', x), c('lb',x,'rb'), c(x, 'rb','rb'), sep='_'))
    } else {
      if (type == 'bi') return(paste(c(x[-1]), x[-length(x)], sep='_'))
      if (type == 'tri') return(paste(x[-(1:2)], x[-length(x)][-1], x[-(length(x) - c(1,0))], sep='_'))
    }
  })
  n = sapply(uni, length)
  bi = data.frame(bigram=unlist(uni), i = rep(1:length(uni), times=n))
  bi_voc = unique(bi$bigram)
  bi = Matrix::spMatrix(nrow = length(x), ncol = length(bi_voc), 
                   i = bi$i, j = match(bi$bigram, bi_voc), x = rep(1, nrow(bi)))
  rownames(bi) = voc
  methods::as(bi, 'CsparseMatrix')
}


#' Find terms with similar spelling
#'
#' A quick, language agnostic way for finding terms with similar spelling. 
#' Calculates similarity as percentage of a terms bigram's or trigram's that also occur in the other term. 
#' The percentage has to be above the given threshold for both terms (unless allow_asym = T)  
#'
#' @param voc            A character vector that gives the vocabulary (e.g., colnames of a dtm)
#' @param type           Either "bi" (bigrams) or "tri" (trigrams)
#' @param min_overlap    The minimal overlap percentage. Works together with max_diff to determine required overlap
#' @param max_diff       The maximum number of bi/tri-grams that is different
#' @param pad            If True, pad the left size (ls) and right side (rs) of bi/tri-grams. So, trigrams for "pad" would be: "ls_ls_p", "ls_p_a", "p_a_d", "a_d_rs", "d_rs_rs".
#' @param as_lower       If True, ignore case
#' @param same_start     Should terms start with the same character(s)? Given as a number for the number of same characters. (also greatly speeds up calculation)
#' @param drop_non_alpha If True, ignore non alpha terms (e.g., numbers, punctuation). They will appear in the output matrix, but only with zeros.
#' @param min_length     The minimum number of characters in a term. Terms with fewer characters are ignored. They will appear in the output matrix, but only with zeros.
#' @param allow_asym     If True, the match only needs to be true for at least one term. In practice, this means that "America" would match perfectly with "Southern-America".
#' @param verbose        If True, report progress
#'
#' @return  A similarity matrix in the CsparseMatrix format
#' @export
#'
#' @examples
#' dfm = quanteda::tokens(c('That guy Gadaffi','Do you mean Kadaffi?',
#'                          'Nah more like Gadaffel','What Gargamel?')) |>
#'   quanteda::dfm()
#' simmat = term_char_sim(colnames(dfm), same_start=0)
#' term_union(dfm, simmat, verbose = FALSE)
term_char_sim <- function(voc, type=c('tri','bi'), min_overlap=2/3, max_diff=4, pad=F, as_lower=T, same_start=1, drop_non_alpha=T, min_length=5, allow_asym=F, verbose=T) {
  type = match.arg(type)
  if (!methods::is(voc, 'character')) stop('voc must be a character vector')
  if (as_lower) voc = tolower(voc)
  m = char_grams(voc, type=type, pad=pad, drop_non_alpha = drop_non_alpha, min_length=min_length)  ## sparse matrix of bigrams
  
  max_diff_pct = (Matrix::rowSums(m) - max_diff) / Matrix::rowSums(m)
  min_overlap = ifelse(min_overlap < max_diff_pct, max_diff_pct, min_overlap)
  if (same_start > 0) {
    group = substr(rownames(m), start = 0, stop = same_start)
    simmat = tcrossprod_sparse(m, rowsum_div = T, group = group, crossfun = 'min', diag=F, min_value=min_overlap, verbose=verbose)
  } else {
    simmat = tcrossprod_sparse(m, rowsum_div = T, crossfun = 'min', diag=F, min_value=min_overlap, verbose=verbose)
  }
  if (!allow_asym) simmat = tril(simmat>0) * tril(t(simmat>0)) ## both need to match, because otherwise 
  diag(simmat) = 1
  methods::as(simmat, 'CsparseMatrix')
}

#' View term scores for a given document
#'
#' @param dtm      A quanteda dfm
#' @param docname  name of document to select
#' @param doc_i    alternatively, select document by index
#'
#' @return  A named vector with terms (names) and scores
#' @export
#'
#' @examples
#' get_doc_terms(rnewsflow_dfm, doc_i=1)
get_doc_terms <- function(dtm, docname=NULL, doc_i=NULL) {
  if (is.null(docname) && is.null(doc_i)) stop('either docname or doc_i has to be specified')
  if (!is.null(docname) && !is.null(doc_i)) stop('either (not both) docname or doc_i has to be specified')
  if (!is.null(doc_i)) docname = quanteda::docnames(dtm)[doc_i]
  r = dtm[quanteda::docnames(dtm) == docname,]
  if (nrow(r) == 0) stop('docname is not a document in dtm')
  cs = colSums(r)
  cs[cs > 0]
}


#' View overlapping terms for a given pair of documents
#'
#' @param dtm      A quanteda dfm
#' @param doc.x    The name of the first document in dtm
#' @param doc.y    The name of the second document in dtm (or dtm.y) 
#' @param dtm.y    Optionally, a second dtm (for when the documents occur in separate dtm's) 
#'
#' @return  A character vector
#' @export
#'
#' @examples
#' get_overlap_terms(rnewsflow_dfm, 
#'                   quanteda::docnames(rnewsflow_dfm)[1],
#'                   quanteda::docnames(rnewsflow_dfm)[5])
get_overlap_terms <- function(dtm, doc.x, doc.y, dtm.y=dtm) {
  tx = get_doc_terms(dtm, doc.x)
  ty = get_doc_terms(dtm.y, doc.y)
  intersect(names(tx), names(ty))
}



merge_vocabulary_terms <- function(voc) {
  term = left = right = i = NULL  
  l = stringi::stri_split(voc, fixed = ' #&# ')
  ln = sapply(l, length)
  l = unlist(l)
  lr = stringi::stri_split(l, fixed = ' & ')
  lrn = sapply(lr, length)
  d = data.table::data.table(i = rep(1:length(ln), ln),
                             left = sapply(lr, '[[', 1))
  d[lrn > 1, right := sapply(lr[lrn > 1], '[[', 2)]
  d[is.na(right), right := '']
  
  ds = d[, list(left = paste(left, collapse='|')), by=c('i','right')]
  ds = ds[, list(right = paste(right, collapse='|')), by=c('i','left')]
  ds[, term := left]
  ds[!right == '', term := paste(left, right, sep=' & ')]
  ds = ds[, list(term = paste(term, collapse=' & ')), by=c('i')]
  if (!identical(1:length(voc), ds$i)) stop('test error in merge_vocabulary_terms')
  ds$term
}
