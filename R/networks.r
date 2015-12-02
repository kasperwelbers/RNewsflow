library(igraph)
library(data.table)

#' @export
document.network <- function(d, meta, id.var='document_id', date.var='date', min.similarity=0, add.empty.documents=T){
  confirm.dtm.meta(meta, id.var, date.var)
  
  colnames(d) = c('x','y','similarity')
  d = d[d$similarity > min.similarity, c('x','y','similarity')]
  d = d[order(-d$similarity),]
  
  g = graph.data.frame(d[,c('x','y')])
  E(g)$weight = d$similarity
  
  if(mean(V(g)$name %in% meta[,id.var]) < 1) stop("Not all documents in x and y match with an 'id' in the meta information")
  
  if(add.empty.documents){
    missingmeta = as.character(meta[!meta[,id.var] %in% V(g)$name, id.var])
    g = add.vertices(g, nv = length(missingmeta), attr = list(name=missingmeta))
  }
  meta = meta[match(V(g)$name, meta[,id.var]),]
  
  attribs = colnames(meta)[!colnames(meta) == id.var]
  for(attrib in attribs){
    g = set.vertex.attribute(g, attrib, value=as.character(meta[,attrib]))
  }
  
  edgelist = get.edges(g, E(g))
  dates = as.POSIXct(V(g)$date) 
  E(g)$hourdiff = round(difftime(dates[edgelist[,2]], dates[edgelist[,1]], units = 'hours'),3)
  
  g
}

#' @export
plot.document.network <- function(g, date.attribute='date', source.attribute='source', subcomp=NULL, dtm=NULL, fixed.size=F, sources=NULL, only.outer.date=F, date.format='%Y-%m-%d %H:%M', margins=c(5,8,1,13), isolate.color=NULL, source.loops=T, ...){
  g = set.vertex.attribute(g, 'date', value= get.vertex.attribute(g, date.attribute))
  g = set.vertex.attribute(g, 'source', value= get.vertex.attribute(g, source.attribute))
  
  if(is.null(subcomp)){
    cluster = g
  } else {
    clusters = decompose.graph(g)
    cluster = clusters[[subcomp]]
  }
  
  if(is.null(sources)) sources = sort(unique(V(cluster)$source))
  cluster = delete.vertices(cluster, which(!V(cluster)$source %in% sources))
  
  if(!is.null(dtm)){
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE), widths = c(1.5, 2.5), heights = c(1, 2))  
    topwords = col_sums(dtm[rownames(dtm) %in% V(cluster)$name,])
    topwords = topwords[topwords > 0, drop=F]
    topwords = head(topwords[order(-topwords)], 10)
    par(mar = c(0,0,0,0))
    suppressWarnings(
    wordcloud::wordcloud(names(topwords), sqrt(as.numeric(topwords)), scale=c(2,0.7))
    )
  } else {
    layout(matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
  }
  
  plotDocNet(cluster, sources, only.outer.date, date.format, margins, isolate.color=NULL, source.loops=T, ...)
  layout(matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
}

plotDocNet <- function(cluster, sources=NULL, only.outer.date=F, date.format='%Y-%m-%d %H:%M', margins=c(5,8,1,13), isolate.color=NULL, source.loops=T, ...){
  par(mar = margins)
  if(is.null(sources)) sources = sort(unique(V(cluster)$source))
  E(cluster)$width = scales::rescale(E(cluster)$weight, c(2,4), from = c(0,1))
  E(cluster)$arrow.size = 0.45
  E(cluster)$color = grey(scales::rescale(E(cluster)$weight, c(0.8,0.2), from = c(0,1)))
  
  if(!source.loops){
    edgelist = get.edgelist(cluster)
    xmatch = match(edgelist[,1], V(cluster)$name)
    ymatch = match(edgelist[,2], V(cluster)$name)
    sourceloop = which(V(cluster)$source[xmatch] == V(cluster)$source[ymatch])
    cluster = delete.edges(cluster, sourceloop)
  }
  
  V(cluster)$x = as.POSIXct(V(cluster)$date)
  V(cluster)$y = length(sources) + 1 - match(V(cluster)$source, sources)
  
  vertex.atts = names(vertex.attributes(g))
  if(!'size' %in% vertex.atts) V(cluster)$size = 7.5
  if(!'color' %in% vertex.atts) V(cluster)$color = 'lightgrey'   
  if(!'label' %in% vertex.atts) V(cluster)$label = ''
  if(!'label.color' %in% vertex.atts) V(cluster)$label.color = 'black'
  
  if(!is.null(isolate.color)) V(cluster)$color[degree(cluster) == 0] = isolate.color
  
  ymin = 1 - max(V(cluster)$y)*(1/max(V(cluster)$y))
  cluster$layout = layout.norm(cbind(as.numeric(V(cluster)$x), V(cluster)$y), -1, 1, ymin, 1)
  
  if(length(unique(V(cluster)$x)) == 1) cluster$layout[,1] = 0
  plot(cluster, rescale=F, asp=0, ...)
  
  for(source in sources){
    ycoord = cluster$layout[V(cluster)$source == source,2][1]
    text(1.1, ycoord, source, adj=0)
  }
  
  ## determine inner vertical lines (to indicate date)
  inner_date = seq.POSIXt(as.POSIXct(min(V(cluster)$date)), as.POSIXct(max(V(cluster)$date)), by='day')
  inner_date = as.POSIXct(as.Date(inner_date))
  inner_date = inner_date[as.numeric(inner_date) > min(V(cluster)$x) & as.numeric(inner_date) < max(V(cluster)$x)]
  inner_date_i = scales::rescale(as.numeric(inner_date), to=c(-1,1), from=as.numeric(c(min(V(cluster)$x), max(V(cluster)$x))))
  
  clip(-1.05,1,ymin,1.05)
  abline(h=unique(cluster$layout[,2]), col="grey10", lty='dotted')
  
  ## determine outer vertical lines
  outer_date_i = c(-1,1)
  outer_date = c(min(V(cluster)$date), max(V(cluster)$date))
  if(length(unique(outer_date)) == 1){
    outer_date_i = 0
    outer_date = outer_date[1]
  }
  ## draw lines, add text
  abline(v=outer_date_i, col="grey10", lty='dotted')
  abline(v=inner_date_i, col="grey10", lty='dotted')
  
  clip(-2,2,-100,2)
  outer_date_print = format(as.POSIXct(outer_date), date.format)
  if(!only.outer.date & length(inner_date) > 0){
    text(x=outer_date_i, y=ymin-0.2, label=outer_date_print, srt=30, adj=1) 
    inner_date_print = format(as.POSIXct(inner_date), date.format)
    text(x=inner_date_i, y=ymin-0.2, label=inner_date_print, srt=30, adj=1) 
  } else text(x=outer_date_i, y=ymin-0.2, label=outer_date_print, srt=30, adj=1) 
}

#' @export
show.window <- function(g, vertex.attribute=NULL){
  from.edge.i = get.edges(g, E(g))[,1]
  
  if(is.null(vertex.attribute)) {
    edgemed = rep('All vertices', vcount(g))[from.edge.i]
  } else {
    edgemed = get.vertex.attribute(g, vertex.attribute)[from.edge.i]  
  }
  medwindow = aggregate(E(g)$hourdiff, by=list(source=edgemed), FUN = function(x) cbind(min(x), max(x)))
  d = data.frame(source = medwindow$source, window.left=medwindow$x[,1], window.right=medwindow$x[,2])
  colnames(d) = c('vertex.attribute', 'window.left','window.right')
  d
}

#' @export
filter.window <- function(g, hour.window, select.vertices=NULL){
  if(is.null(select.vertices)) select.vertices = rep(T, vcount(g))
  
  # get vector indicating which edges are selected given the vertex.filter
  from.edge.i = get.edges(g, E(g))[,1]
  selected_edge = select.vertices[from.edge.i]
  
  # filter the selected edges
  delete.i = selected_edge & (E(g)$hourdiff < hour.window[1] | E(g)$hourdiff > hour.window[2])
  delete.edges(g, which(delete.i))
}


###
###
#' @export
source.window.settings <- function(g, direct.edit=F, max.window.left=0.5, max.window.right=36, meta.attribute='source'){
  sources = unique(get.vertex.attribute(g, meta.attribute))
  if(class(g) == 'data.frame') sources = unique(g$source)
  
  d = data.frame(source=sources, window.left=max.window.left, window.right=max.window.right)
  colnames(d) = c(meta.attribute, 'window.left','window.right')
  if(direct.edit) {
    d = edit(d)
    message('Created the following source.window.table:')
    print(d)
  }
  d
}

#' @export
filter.window_oud <- function(g, source.window){
  if(is.null(source.window)) source.window = source.window.settings(g, direct.edit=T)
  
  sourcevar = colnames(source.window)[1]
  dsource = get.vertex.attribute(g, sourcevar)[get.edges(g, E(g))[,1]]
  
  delete_i = list()
  for(i in 1:nrow(source.window)){
    source = as.character(source.window[i,1])
    winleft = as.numeric(source.window$window.left[i])
    winright = as.numeric(source.window$window.right[i])
    delete = which(dsource == source & (E(g)$hourdiff < winleft | E(g)$hourdiff > winright))
    delete_i[['']] = delete
  }
  delete.edges(g, unique(unlist(delete_i)))
}
####
###


#' @export
only.first.match <- function(g){
  e = data.frame(get.edges(g, E(g)))
  e$i = 1:nrow(e)
  e$to.date = V(g)$date[e$X2]
  e = e[order(e$to.date),]
  
  first = !duplicated(e$X1)
  first = paste(e$X1, e$to.date, sep='.') %in% paste(e$X1[first], e$to.date[first], sep='.') # also first if same date as first
  
  delete.edges(g, e$i[!first])
}

#' @export
only.strongest.match <- function(g, by.vertex.meta=NULL){
  e = data.frame(get.edges(g, E(g)))
  e$i = 1:nrow(e)

  if(!is.null(by.vertex.meta)) {
    v = get.data.frame(g, 'vertices')
    e$by = apply(v[,by.vertex.meta, drop=F], 1, paste, collapse=';')[e$X2]
  } else {
    e$by = 1
  }
  
  e = e[order(-E(g)$weight),]
  e = e[duplicated(e[,c('X1', 'by')]),]
  delete.edges(g, e$i)
}

#' @export
aggregate.network <- function(g, by=NULL, by.from=by, by.to=by, edge.attribute='weight', agg.FUN=median, return.df=F){
  V(g)$all_vertices = 'all_vertices'
  if(is.null(by.from)) by.from = 'all_vertices'
  if(is.null(by.to)) by.to = 'all_vertices'
  
  e = data.frame(get.edges(g, E(g)))
  v = get.data.frame(g, 'vertices')  
  
  e = cbind(e, v[e$X1, by.from, drop=F])
  e = cbind(e, v[e$X2, by.to, drop=F])
  colnames(e) = c('from','to', paste('from', by.from, sep='.'), paste('to', by.to, sep='.'))
  aggvars = colnames(e)[3:ncol(e)]
  
  e$from.unique_match = as.numeric(!duplicated(e[,colnames(e)[-2]]))
  e$to.unique_match = as.numeric(!duplicated(e[,colnames(e)[-1]]))
  e$var = get.edge.attribute(g, edge.attribute)
  e = setDT(e)[, .(edges=.N,
                   edge.attribute=agg.FUN(var),
                   from.V=sum(from.unique_match),
                   to.V=sum(to.unique_match)),
                   by= eval(paste(aggvars, collapse=','))]
  e = as.data.frame(e)
  colnames(e)[colnames(e) == 'edge.attribute'] = paste('agg',edge.attribute, sep='.')
  
  # match total nodes
  from.totals = setDT(v)[, .(from.N=.N), by= eval(paste(by.from, collapse=','))]
  to.totals = setDT(v)[, .(to.N=.N), by= eval(paste(by.to, collapse=','))]
  
  e = merge(e, from.totals, by.x=paste('from', by.from, sep='.'), by.y=by.from, all.x=T)
  e = merge(e, to.totals, by.x=paste('to', by.to, sep='.'), by.y=by.to, all.x=T)
  e$from.Vprop = e$from.V / e$from.N
  e$to.Vprop = e$to.V / e$to.N
  
  if(!return.df) return(return.aggregate.network(e))  
  if(return.df) {
    e = e[,!colnames(e) %in% c('from.N', 'to.N')]
    return(e) 
  }
}

return.aggregate.network <- function(g.df){
  vnames = colnames(g.df)[1:(grep('^edges$', colnames(g.df))-1)]
  
  # create network
  from.v = g.df[,grep('^from\\.', vnames), drop=F]
  to.v = g.df[,grep('^to\\.', vnames), drop=F]
  from.name = if(ncol(from.v) > 1) apply(from.v, MARGIN = 1, paste, collapse='; ') else from.v
  to.name = if(ncol(to.v) > 1) apply(to.v, MARGIN = 1, paste, collapse='; ') else to.v
  g.agg = graph.data.frame(cbind(from.name,to.name))
  
  aggvar = colnames(g.df)[grep('agg\\.', colnames(g.df))]
  for(edge.attrib in c('edges',aggvar,'from.V','from.Vprop','to.V','to.Vprop')){
    g.agg = set.edge.attribute(g.agg, edge.attrib, value=g.df[,edge.attrib])
  }
  
  # create meta table
  from.v = cbind(from.name, from.v, from.N=g.df$from.N)
  to.v = cbind(to.name, to.v, to.N=g.df$to.N)
  
  colnames(from.v) = gsub('^from\\.', '', colnames(from.v))
  colnames(to.v) = gsub('^to\\.', '', colnames(to.v))
  metavars = unique(c(colnames(from.v), colnames(to.v)))
  for(metavar in metavars){
    if(!metavar %in% colnames(from.v)) from.v[,metavar] = NA
    if(!metavar %in% colnames(to.v)) to.v[,metavar] = NA
  }

  meta = unique(rbind(from.v[,metavars], to.v[,metavars]))
  
  # add meta attributes
  meta_i = match(V(g.agg)$name, meta$name)
  for(metavar in metavars[!metavars == 'name']){
    g.agg = set.vertex.attribute(g.agg, metavar, value=meta[meta_i,metavar])
  }
  
  g.agg
}

graph.plot.presets <- function(g){
  E(g)$curved = 0.2
  E(g)$width = E(g)$weight * 5
  E(g)$arrow.size = 0.7
  E(g)$label = round(E(g)$weight, 2)
  
  V(g)$label.color = 'black'
  V(g)$color = 'white'
  V(g)$size = 10
  g
}

#' @export
plot.aggregate.network <- function(g, weight.var='from.Vprop', weight.thres=NULL, delete.isolates=F, 
                                   vertex.size=30, vertex.color='lightblue', vertex.label.color='black', vertex.label.cex=0.7, 
                                   edge.color = 'grey', show.edge.labels=T, edge.label.color = 'black', edge.label.cex = 0.6, edge.arrow.size=1, 
                                   layout=layout.davidson.harel, ...){
  E(g)$weight = get.edge.attribute(g, 'from.Vprop')
  if(!is.null(weight.thres)) g = delete.edges(g, which(E(g)$weight < weight.thres))
  if(delete.isolates) g = delete.vertices(g, degree(g) == 0)
  
  argnames = names(list(...))  
  E(g)$color = edge.color
  E(g)$arrow.size = edge.arrow.size
  E(g)$label.cex = edge.label.cex
  E(g)$label.color = edge.label.color
  if(!'edge.width' %in% argnames) E(g)$width = scales::rescale(E(g)$weight, to = c(1,5))
  if(!'edge.label' %in% argnames & show.edge.labels) E(g)$label = round(E(g)$weight, 2)
  if(!'edge.label.font' %in% argnames) E(g)$label.font = 3
  if(!'edge.curved' %in% argnames) {
    e = get.edges(g, E(g))
    twoway = paste(e[,1], e[,2], sep='.') %in% paste(e[,2], e[,1], sep='.')
    E(g)$curved = ifelse(twoway, 0.35, 0.15)
  } 
  
  V(g)$size = vertex.size
  V(g)$label.cex = vertex.label.cex
  V(g)$color = vertex.color
  V(g)$label.color = vertex.label.color
  if(!'vertex.label' %in% argnames) V(g)$label = gsub(' ', '\n', V(g)$name)  
  if(!'vertex.frame.color' %in% argnames) V(g)$frame.color = V(g)$color
  
  g$layout = layout
  plot(g, ...)
}

