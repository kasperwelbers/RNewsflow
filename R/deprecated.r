is_deprecated <- function(new, old = as.character(sys.call(sys.parent()))[1L]){
  msg <- gettextf("The '%s' function is deprecated, replaced by '%s'. It will be deleted in the next update.\nThis is part of a broader change towards exclusive support for quanteda style DTMs (to improve clarity) and better naming/documentation",old,new)
  .Deprecated(msg=msg)
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
  is_deprecated('term_day_dist')
  
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

#' A wrapper for \link[igraph]{plot.igraph} for visualizing directed networks.
#'
#' This is a convenience function for visualizing directed networks with edge labels using \link[igraph]{plot.igraph}. 
#' It was designed specifically for visualizing aggregated document similarity networks in the RNewsflow package, but works with any network in the \link[igraph]{igraph} class.  
#'
#' @param g A network/graph in the \link[igraph]{igraph} class 
#' @param weight.var The edge attribute that is used to specify the edges
#' @param weight.thres A threshold for weight. Edges below the threshold are ignored 
#' @param delete.isolates If TRUE, isolates (i.e. vertices without edges) are ignored.
#' @param vertex.size The size of the verticex/nodes. Defaults to 30. Can be a vector with values per vertex.
#' @param vertex.color Color of vertices/nodes. Default is lightblue. Can be a vector with values per vertex.
#' @param vertex.label.color Color of labels for vertices/nodes. Defaults to black. Can be a vector with values per vertex. 
#' @param vertex.label.cex Size of the labels for vertices/nodes. Defaults to 0.7. Can be a vector with values per vertex. 
#' @param edge.color Color of the edges. Defaults to grey. Can be a vector with values per edge.
#' @param show.edge.labels Logical. Should edge labels be displayed? Default is TRUE.
#' @param edge.label.color Color of the edge labels. Defaults to black. Can be a vector with values per edge. 
#' @param edge.label.cex Size of the edge labels. Defaults to 0.6. Can be a vector with values per edge. 
#' @param edge.arrow.size Size of the edge arrows. Defaults to 1. Can only be set globally (igraph might update this at some point)
#' @param layout The igraph layout used to plot the network. Defaults to \link[igraph]{layout.davidson.harel}
#' @param ... Arguments to be passed to the \link[igraph]{plot.igraph} function. 
#'
#' @return Nothing
#' @export
#'
#' @examples
#' data(docnet)
#' aggdocnet = network_aggregate(docnet, by='source')
#' directed.network.plot(aggdocnet, weight.var = 'to.Vprop', weight.thres = 0.2)
directed.network.plot <- function(g, weight.var='from.Vprop', weight.thres=NULL, delete.isolates=FALSE, 
                                  vertex.size=30, vertex.color='lightblue', vertex.label.color='black', vertex.label.cex=0.7, 
                                  edge.color = 'grey', show.edge.labels=TRUE, edge.label.color = 'black', edge.label.cex = 0.6, edge.arrow.size=1, 
                                  layout=igraph::layout.davidson.harel, ...){
  is_deprecated('directed_network_plot')
  igraph::E(g)$weight = igraph::get.edge.attribute(g, weight.var)
  if(!is.null(weight.thres)) g = igraph::delete.edges(g, which(igraph::E(g)$weight < weight.thres))
  if(delete.isolates) g = igraph::delete.vertices(g, igraph::degree(g) == 0)
  
  argnames = names(list(...))  
  igraph::E(g)$color = edge.color
  igraph::E(g)$arrow.size = edge.arrow.size
  igraph::E(g)$label.cex = edge.label.cex
  igraph::E(g)$label.color = edge.label.color
  if(!'edge.width' %in% argnames) igraph::E(g)$width = scales::rescale(igraph::E(g)$weight, to = c(1,5))
  if(!'edge.label' %in% argnames & show.edge.labels) igraph::E(g)$label = round(igraph::E(g)$weight, 2)
  if(!'edge.label.font' %in% argnames) igraph::E(g)$label.font = 3
  if(!'edge.curved' %in% argnames) {
    e = igraph::get.edges(g, igraph::E(g))
    twoway = paste(e[,1], e[,2], sep='.') %in% paste(e[,2], e[,1], sep='.')
    igraph::E(g)$curved = ifelse(twoway, 0.35, 0.15)
  } 
  
  igraph::V(g)$size = vertex.size
  igraph::V(g)$label.cex = vertex.label.cex
  igraph::V(g)$color = vertex.color
  igraph::V(g)$label.color = vertex.label.color
  if(!'vertex.label' %in% argnames) igraph::V(g)$label = gsub(' ', '\n', igraph::V(g)$name)  
  if(!'vertex.frame.color' %in% argnames) igraph::V(g)$frame.color = igraph::V(g)$color
  
  g$layout = layout
  igraph::plot.igraph(g, ...)
}



#' Aggregate the edges of a network by vertex attributes 
#' 
#' This function offers a versatile way to aggregate the edges of a network based on the vertex attributes.
#' Although it was designed specifically for document similarity networks, it can be used for any network in the \link[igraph]{igraph} class.
#' 
#' The first argument is the network (in the `igraph` class). 
#' The second argument, for the `by` parameter, is a character vector to indicate one or more vertex attributes based on which the edges are aggregated.
#' Optionally, the `by` parameter can also be specified separately for `by.from` and `by.to`. 
#' 
#' By default, the function returns the aggregated network as an \link[igraph]{igraph} class.
#' The edges in the aggregated network have five standard attributes. 
#' The `edges` attribute counts the number of edges from the `from` group to the `to` group. 
#' The `from.V` attribute shows the number of vertices in the `from` group that matched with a vertex in the `to` group.
#' The `from.Vprop attribute shows this as the proportion of all vertices in the `from` group.
#' The `to.V` and `to.Vprop` attributes show the same for the `to` group.
#' 
#' In addition, one of the edge attributes of the original network can be aggregated with a given function.
#' These are specified in the `edge.attribute` and `agg.FUN` parameters.
#'
#' @param g A network/graph in the \link[igraph]{igraph} class
#' @param by A character string indicating the vertex attributes by which the edges will be aggregated. 
#' @param by.from Optionally, specify different vertex attributes to aggregate the `from` side of edges 
#' @param by.to Optionally, specify different vertex attributes to aggregate the `to` side of edges 
#' @param edge.attribute Select an edge attribute to aggregate using the function specified in `agg.FUN`. Defaults to 'weight'
#' @param agg.FUN The function used to aggregate the edge attribute
#' @param return.df Optional. If TRUE, the results are returned as a data.frame. This can in particular be convenient if by.from and by.to are used.
#' @param keep_isolates if True, also return scores for isolates
#'
#' @return A network/graph in the \link[igraph]{igraph} class, or a data.frame if return.df is TRUE.
#' @export
#'
#' @import data.table
#'
#' @examples
#' data(docnet)
#' aggdocnet = network.aggregate(docnet, by='sourcetype')
#' igraph::get.data.frame(aggdocnet, 'both')
#' 
#' aggdocdf = network.aggregate(docnet, by.from='sourcetype', by.to='source', return.df = TRUE)
#' head(aggdocdf)
network.aggregate <- function(g, by=NULL, by.from=by, by.to=by, edge.attribute='weight', agg.FUN=mean, return.df=FALSE, keep_isolates=T){
  is_deprecated('network_aggregate')
  
  igraph::V(g)$any_vertex = '---'
  if(is.null(by.from)) by.from = 'any_vertex'
  if(is.null(by.to)) by.to = 'any_vertex'
  
  e = data.frame(igraph::get.edges(g, igraph::E(g)))
  v = igraph::get.data.frame(g, 'vertices')  
  
  e = cbind(e, v[e$X1, by.from, drop=FALSE])
  e = cbind(e, v[e$X2, by.to, drop=FALSE])
  colnames(e) = c('from','to', paste('from', by.from, sep='.'), paste('to', by.to, sep='.'))
  aggvars = colnames(e)[3:ncol(e)]
  
  e$from_unique_match = as.numeric(!duplicated(e[,colnames(e)[-2]]))
  e$to_unique_match = as.numeric(!duplicated(e[,colnames(e)[-1]]))
  e$var = igraph::get.edge.attribute(g, edge.attribute)
  
  e = data.table::setDT(e)
  e = e[, list(edges=.N,
               edge.attribute=agg.FUN(get('var')), # note that get('var') is used because the CRAN check would otherwise wrongfully interpret 'var' as a global variable (yielding a NOTE)  
               from.V=sum(get('from_unique_match')),
               to.V=sum(get('to_unique_match'))),
        by= eval(paste(aggvars, collapse=','))]
  e = as.data.frame(e)
  
  colnames(e)[colnames(e) == 'edge.attribute'] = paste('agg',edge.attribute, sep='.')
  
  # match total vertices
  v = data.table::setDT(v)
  from.totals = v[, list(from.N=.N), by= eval(paste(by.from, collapse=','))]
  to.totals = v[, list(to.N=.N), by= eval(paste(by.to, collapse=','))]
  
  e = merge(e, from.totals, by.x=paste('from', by.from, sep='.'), by.y=by.from, all.x=TRUE)
  e = merge(e, to.totals, by.x=paste('to', by.to, sep='.'), by.y=by.to, all.x=TRUE)
  e$from.Vprop = e$from.V / e$from.N
  e$to.Vprop = e$to.V / e$to.N
  
  if(!return.df) {
    e = return_network_aggregate(e, from.totals, to.totals, keep_isolates)
    if('any_vertex' %in% names(igraph::vertex.attributes(e))) e = igraph::delete_vertex_attr(e, 'any_vertex')
    return(e)
  }
  if(return.df) {
    e = e[,!colnames(e) %in% c('from.N', 'to.N', 'from.any_vertex', 'to.any_vertex')]
    
    edge_col = which(colnames(e) == 'edges') # silly but necessary way to order the from and to attribute columns ()
    ord = c(order(colnames(e)[1:(edge_col-1)]), edge_col:ncol(e))
    return(e[,ord]) 
  }
}


#' Filter edges from the document similarity network based on hour difference
#'
#' The `filter.window` function can be used to filter the document pairs (i.e. edges) using the `hour.window` parameter, which works identical to the `hour.window` parameter in the `newsflow.compare` function. 
#' In addition, the `from.vertices` and `to.vertices` parameters can be used to select the vertices (i.e. documents) for which this filter is applied.
#'
#' It is recommended to use the \link[RNewsflow]{show_window} function to verify whether the hour windows are correct according to the assumptions and focus of the study.
#'
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}  
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param to.vertices A filter to select the vertices `to` which an edge is filtered. 
#' For example, if `V(g)$sourcetype == "newspaper"` is used, then the hour.window filter is only applied for edges `to` newspaper documents (specifically, where the sourcetype attribute is "newspaper"). 
#' @param from.vertices A filter to select the vertices `from` which an edge is filtered. 
#' Works identical to `to.vertices`.  
#'
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#'
#' @examples
#' data(docnet)
#' show.window(docnet, to.attribute = 'source') # before filtering
#' 
#' docnet = filter.window(docnet, hour.window = c(0.1,24))
#' 
#' docnet = filter.window(docnet, hour.window = c(6,36), 
#'                        to.vertices = V(docnet)$sourcetype == 'Print NP')
#' 
#' show.window(docnet, to.attribute = 'sourcetype') # after filtering per sourcetype
#' show.window(docnet, to.attribute = 'source') # after filtering per source
filter.window <- function(g, hour.window, to.vertices=NULL, from.vertices=NULL){
  is_deprecated('filter_window')
  
  if(is.null(from.vertices)) from.vertices = rep(TRUE, igraph::vcount(g))
  if(is.null(to.vertices)) to.vertices = rep(TRUE, igraph::vcount(g))
  
  # get vector indicating which edges are selected given the vertex.filter
  from.edge.i = igraph::get.edges(g, igraph::E(g))[,1]
  to.edge.i = igraph::get.edges(g, igraph::E(g))[,2]
  selected_edge = from.vertices[from.edge.i] & to.vertices[to.edge.i]
  
  # filter the selected edges
  delete.i = selected_edge & (igraph::E(g)$hourdiff < hour.window[1] | igraph::E(g)$hourdiff > hour.window[2])
  igraph::delete.edges(g, which(delete.i))
}

#' Show time window of document pairs 
#'
#' This function aggregates the edges for all combinations of attributes specified in `from.attribute` and `to.attribute`, and shows the minimum and maximum hour difference for each combination.
#'
#' The \link[RNewsflow]{filter.window} function can be used to filter edges that fall outside of the intended time window. 
#'
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}  
#' @param to.attribute The vertex attribute to aggregate the `to` group of the edges
#' @param from.attribute The vertex attribute to aggregate the `from` group of the edges
#'
#' @return A data.frame showing the left and right edges of the window for each unique group.
#' @export
#'
#' @examples
#' data(docnet)
#' show.window(docnet, to.attribute = 'source')
#' show.window(docnet, to.attribute = 'sourcetype')
#' show.window(docnet, to.attribute = 'sourcetype', from.attribute = 'sourcetype')
show.window <- function(g, to.attribute=NULL, from.attribute=NULL){
  is_deprecated('show_window')
  
  from.vertices = if(is.null(from.attribute)) rep('any document', igraph::vcount(g)) else igraph::get.vertex.attribute(g, from.attribute)
  to.vertices = if(is.null(to.attribute)) rep('any document', igraph::vcount(g)) else igraph::get.vertex.attribute(g, to.attribute)
  
  e = data.frame(igraph::get.edges(g, igraph::E(g)))
  medwindow = stats::aggregate(igraph::E(g)$hourdiff, by=list(from=from.vertices[e$X1], to=to.vertices[e$X2]), FUN = function(x) cbind(min(x), max(x)))
  d = data.frame(from = medwindow$from, to=medwindow$to, window.left=medwindow$x[,1], window.right=medwindow$x[,2])
  d
}

#' Transform document network so that each document only matches the earliest dated matching document
#'
#' Transforms the network so that a document only has an edge to the earliest dated document it matches within the specified time window[^duplicate].
#' 
#' If there are multiple earliest dated documents (that is, having the same publication date) then edges to all earliest dated documents are kept.
#' 
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}  
#'
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#'
#' @examples
#' data(docnet)
#' 
#' subcomp1 = igraph::decompose.graph(docnet)[[2]]
#' subcomp2 = only.first.match(subcomp1)
#' 
#' igraph::get.data.frame(subcomp1)
#' igraph::get.data.frame(subcomp2)
#' 
#' graphics::par(mfrow=c(2,1))
#' document.network.plot(subcomp1, main='All matches')
#' document.network.plot(subcomp2, main='Only first match')
#' graphics::par(mfrow=c(1,1))
only.first.match <- function(g){
  is_deprecated('only_first_match')
  e = data.frame(igraph::get.edges(g, igraph::E(g)))
  
  # first make from and to columns of same temporal order (from before to) by switching columns with negative hourdiff.
  switch_from.to = igraph::get.data.frame(g, 'edges')$hourdiff < 0
  e = data.frame(from = ifelse(switch_from.to, e$X2, e$X1),
                 to   = ifelse(switch_from.to, e$X1, e$X2))
  
  e$i = 1:nrow(e)
  e$from.date = igraph::V(g)$date[e$from]
  e = e[order(e$from.date),]
  
  first = !duplicated(e$to)
  first = paste(e$to, e$from.date, sep='.') %in% paste(e$to[first], e$from.date[first], sep='.') # also first if same date as first
  
  igraph::delete.edges(g, e$i[!first])
}

#' Visualize (a subcomponent) of the document similarity network 
#'
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}
#' @param date.attribute The label of the vertex/document date attribute. Default is "date"
#' @param source.attribute The label of the vertex/document source attribute. Default is "source"
#' @param subcomp_i Optional. If an integer is given, the network is decomposed into subcomponents (i.e. unconnected components) and only the i-th component is visualized.
#' @param dtm Optional. If a document-term matrix that contains the documents in g is given, a wordcloud with the most common words of the network is added. 
#' @param sources Optional. Use a character vector to select only certain sources
#' @param only.outer.date If TRUE, only the labels for the first and last date are reported on the x-axis
#' @param date.format The date format of the date labels (see \link[base]{format.POSIXct})
#' @param margins The margins of the network plot. The four values represent bottom, left, top and right margin. 
#' @param isolate.color Optional. Set a custom color for isolates
#' @param source.loops If set to FALSE, all edges between vertices/documents of the same source are ignored.
#' @param ... Additional arguments for the network plotting function \link[igraph]{plot.igraph} 
#'
#' @return Nothing. 
#' @export
#'
#' @examples
#' docnet = docnet
#' dtm = rnewsflow_dfm
#' 
#' docnet_comps = igraph::decompose.graph(docnet) # get subcomponents
#' 
#' # subcomponent 1
#' document.network.plot(docnet_comps[[1]]) 
#' 
#' # subcomponent 2 with wordcloud
#' document.network.plot(docnet_comps[[2]], dtm=dtm) 
#' 
#' # subcomponent 3 with additional arguments passed to plot.igraph 
#' document.network.plot(docnet_comps[[3]], dtm=dtm, vertex.color='red') 
document.network.plot <- function(g, date.attribute='date', source.attribute='source', subcomp_i=NULL, dtm=NULL, sources=NULL, only.outer.date=FALSE, date.format='%Y-%m-%d %H:%M', margins=c(5,8,1,13), isolate.color=NULL, source.loops=TRUE, ...){
  is_deprecated('document_network_plot')
  g = igraph::set.vertex.attribute(g, 'date', value= igraph::get.vertex.attribute(g, date.attribute))
  g = igraph::set.vertex.attribute(g, 'source', value= igraph::get.vertex.attribute(g, source.attribute))
  
  if(is.null(subcomp_i)){
    cluster = g
  } else {
    clusters = igraph::decompose.graph(g)
    cluster = clusters[[subcomp_i]]
  }
  
  if(is.null(sources)) sources = sort(unique(igraph::V(cluster)$source))
  cluster = igraph::delete.vertices(cluster, which(!igraph::V(cluster)$source %in% sources))
  
  if(!is.null(dtm)){
    dtm = quanteda::as.dfm(dtm)
    graphics::layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE), widths = c(1.5, 2.5), heights = c(1, 2))  
    topwords = Matrix::colSums(dtm[rownames(dtm) %in% igraph::V(cluster)$name,])
    topwords = topwords[topwords > 0, drop=FALSE]
    topwords = head(topwords[order(-topwords)], 10)
    graphics::par(mar = c(0,0,0,0))
    suppressWarnings(
      wordcloud::wordcloud(names(topwords), sqrt(as.numeric(topwords)), scale=c(2,0.7))
    )
  } 
  
  plot_doc_net(cluster, sources, only.outer.date, date.format, margins, isolate.color=NULL, source.loops=TRUE, ...)
  if(!is.null(dtm)) graphics::layout(matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
}


#' Create a document similarity network
#'
#' Combines document similarity data (d) with document meta data (meta) into an \link[igraph]{igraph} network/graph.
#' 
#' This function is mainly offered to mimic the output of the \link[RNewsflow]{newsflow.compare} function when using imported document similarity data.
#' This way the functions for transforming, aggregating and visualizing the document similarity data can be used.
#'
#' @param d A data.frame with three columns, that represents an edgelist with weight values. 
#' The first and second column represent the names/ids of the 'from' and 'to' documents/vertices. 
#' The third column represents the similarity score. Column names are ignored
#' @param meta A data.frame where rows are documents and columns are document meta information. 
#' Should at least contain 2 columns: the document name/id and date. 
#' The name/id column should match the document names/ids of the edgelist, and its label is specified in the `id.var` argument. 
#' The date column should be intepretable with \link[base]{as.POSIXct}, and its label is specified in the `date.var` argument.            
#' @param id.var The label for the document name/id column in the `meta` data.frame. Default is "document_id"
#' @param date.var The label for the document date column in the `meta` data.frame . default is "date"
#' @param min.similarity For convenience, ignore all edges where the weight is below `min.similarity`. 
#'
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#'
#' @examples
#' d = data.frame(x = c(1,1,1,2,2,3),
#'                y = c(2,3,5,4,5,6),
#'                v = c(0.3,0.4,0.7,0.5,0.2,0.9))
#' 
#' meta = data.frame(document_id = 1:8,
#'                   date = seq.POSIXt(from = as.POSIXct('2010-01-01 12:00:00'), 
#'                          by='hour', length.out = 8),
#'                   medium = c(rep('Newspapers', 4), rep('Blog', 4)))
#' 
#' g = document.network(d, meta)
#' 
#' igraph::get.data.frame(g, 'both')
#' igraph::plot.igraph(g)
document.network <- function(d, meta, id.var='document_id', date.var='date', min.similarity=NA){
  is_deprecated(new='create_document_network')
  if (!date.var %in% colnames(meta)) stop('The name specified in date.var is not a column in meta')
  if (!id.var %in% colnames(meta)) stop('The name specified in id.var is not a column in meta')
  
  if (nrow(d) == 0) d = data.frame(x=numeric(), y=numeric(), similarity=numeric())
  
  colnames(d) = c('x','y','similarity')
  
  if (!is.na(min.similarity)) d = d[d$similarity >= min.similarity, c('x','y','similarity')]
  if (is.data.table(d)) {
    d = data.table::setorderv(d, cols='similarity', order = -1)
  } else {
    d = d[order(-d$similarity),]
  }
  
  g = igraph::graph.data.frame(d[,c('x','y')])
  igraph::E(g)$weight = d$similarity
  if (nrow(d) > 0)
    if (!all(igraph::V(g)$name %in% meta[[id.var]])) stop("Not all documents in d match with an 'id' in the meta information")
  
  ## add documents in meta data.frame that do not appear in the edgelist (in other words, isolates)
  missingmeta = as.character(meta[!meta[,id.var] %in% igraph::V(g)$name, id.var])
  g = igraph::add.vertices(g, nv = length(missingmeta), attr = list(name=missingmeta))
  
  meta = meta[match(igraph::V(g)$name, meta[,id.var]),]
  attribs = colnames(meta)[!colnames(meta) == id.var]
  for(attrib in attribs){
    g = igraph::set.vertex.attribute(g, attrib, value=as.character(meta[,attrib]))
  }
  
  edgelist = igraph::get.edges(g, igraph::E(g))
  dates = as.POSIXct(igraph::V(g)$date) 
  #igraph::E(g)$hourdiff = round(difftime(dates[edgelist[,2]], dates[edgelist[,1]], units = 'hours'),3)
  dates = as.numeric(dates)
  igraph::E(g)$hourdiff = round((dates[edgelist[,2]] - dates[edgelist[,1]]) / (60*60), 3)
  g
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
#' @param max_p       A threshold for maximium p value. 
#' @param pvalue      If max_p < 1, edges are removed based on a p value. For each document in dtm, a p value is calculated over its outward edges. 
#'                    Default is the p-value based on uniform distribution, akin to a "disparity" filter (see \href{https://www.pnas.org/content/106/16/6483.full}{Serrano et al.}) but without filtering on inward edges.
#' @param simmat If softcosine is used, a symmetrical matrix with the similarity scores of terms. If NULL, the cosine similarity of terms in dtm will be used
#' @param simmat_thres If softosine is used, a threshold for the similarity scores of terms
#' 
#' @return A data frame with pairs of documents and their similarities. 
#' @export
#' 
#' @import tm
#' 
#' @examples
#' comp = documents.compare(rnewsflow_dfm, min.similarity=0.4)
#' head(comp)
documents.compare <- function(dtm, dtm.y=NULL, measure=c('cosine','overlap_pct','overlap','crossprod','softcosine','query_lookup','query_lookup_pct'), 
                              min.similarity=0, n.topsim=NULL, max_p=1, pvalue=c("none", "normal", "lognormal", "nz_normal", "nz_lognormal", "disparity"), 
                              simmat=NULL, simmat_thres=NULL) {  
  
  is_deprecated(new = 'compare_documents')
  dtm = quanteda::as.dfm(dtm)
  measure = match.arg(measure)
  pvalue=match.arg(pvalue)
  
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
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, pvalue=pvalue, max_p=max_p, normalize='l2', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, max_p=max_p, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'overlap') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, max_p=max_p, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'crossprod') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, max_p=max_p, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'softcosine') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, max_p=max_p, normalize='softl2', crossfun = 'softprod', min_value = min.similarity, top_n = n.topsim, diag=diag, simmat=simmat, simmat_thres=simmat_thres)
  if (measure == 'query_lookup') {
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, max_p=max_p, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, 
                           simmat=simmat, simmat_thres=simmat_thres)
  }
  if (measure == 'query_lookup_pct') {
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, max_p=max_p, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, 
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
#' @param max_p       A threshold for maximium p value. 
#' @param pvalue      If max_p < 1, edges are removed based on a p value. For each document in dtm, a p value is calculated over its outward edges. 
#'                    Default is the p-value based on uniform distribution, akin to a "disparity" filter (see \href{https://www.pnas.org/content/106/16/6483.full}{Serrano et al.}) but without filtering on inward edges.
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
#' dtm = quanteda::dfm_tfidf(rnewsflow_dfm)
#' g = newsflow.compare(dtm, hour.window = c(0.1, 36))
#' 
#' vcount(g) # number of documents, or vertices
#' ecount(g) # number of document pairs, or edges
#' 
#' head(igraph::get.data.frame(g, 'vertices'))
#' head(igraph::get.data.frame(g, 'edges'))
newsflow.compare <- function(dtm, dtm.y=NULL, meta=NULL, meta.y=NULL, date.var='date', hour.window=c(-24,24), group.var=NULL, measure=c('cosine','overlap_pct','overlap','crossprod','softcosine','query_lookup','query_lookup_pct'), 
                             min.similarity=0, n.topsim=NULL, only.from=NULL, only.to=NULL, only.complete.window=TRUE, pvalue = c("disparity", "normal", "lognormal", "nz_normal", "nz_lognormal"), max_p=1, return_as = c('igraph','edgelist','matrix'), 
                             batchsize=1000, simmat=NULL, simmat_thres=NULL, margin_attr=T, verbose=FALSE){
  
  
  is_deprecated(new='newsflow_compare')
  
  if (margin_attr) {
    row_attr = T
    col_attr = T
    lag_attr = T
  }
  pvalue = match.arg(pvalue)
  
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
  if (measure == 'cosine') cp = tcrossprod_sparse(dtm, dtm.y, pvalue=pvalue, max_p=max_p, normalize='l2', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2=group.y,
                                                  date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                  row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'overlap_pct') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, max_p=max_p, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                       date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                       row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'overlap') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, max_p=max_p, crossfun = 'min', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                   date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                   row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'crossprod') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, max_p=max_p, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                     date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                     row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  if (measure == 'query_lookup') {
    if (is.null(dtm.y)) {
      dtm.y=dtm
      date.y = date
      group.y = group
    }
    if (length(unique(dtm.y@x)) != 1) dtm.y = methods::as(dtm.y>0, 'dgCMatrix')
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = F, pvalue=pvalue, max_p=max_p, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
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
    cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, max_p=max_p, crossfun = 'prod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                           date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                           row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  }
  
  if (measure == 'softcosine') cp = tcrossprod_sparse(dtm, dtm.y, rowsum_div = T, pvalue=pvalue, max_p=max_p, normalize='softl2', crossfun = 'softprod', min_value = min.similarity, top_n = n.topsim, diag=diag, group=group, group2 = group.y,
                                                      date = date, date2 = date.y, lwindow = hour.window[1], rwindow = hour.window[2], date_unit = 'hours', batchsize=batchsize, simmat=simmat, simmat_thres=simmat_thres, 
                                                      row_attr=row_attr, col_attr=col_attr, lag_attr=lag_attr, verbose=verbose)
  
  
  ## meta data is returned as data.table 
  meta = data.table::as.data.table(meta)
  if (!is.null(meta.y)) meta.y = data.table::as.data.table(meta.y)
  if (return_as %in% c('matrix','edgelist')) {
    ## if output is matrix or edgelist it is necessary to split meta and meta.x for matching the margin attributes
    ## however, for igraph is it necessary to keep a single meta if no meta.y exists to prevent creating duplicate vertex meta (with row and col attributes)
    meta.y = meta
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
    diag(cp) = 0
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
#' This function is deprecated, and will at some point be removed. It is replaced by delete_duplicates.
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
  is_deprecated(new='delete_duplicates')
  
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

