## ---- include=FALSE------------------------------------------------------
options(digits=3)
library(knitr)

## ---- message=FALSE, warning=FALSE, echo=TRUE----------------------------
library(quanteda)
d = data.frame(id = c(1,2,3),
               text = c('Socrates is human', 'Humans are mortal', 'Therefore, Socrates is mortal'),
               author = c('Aristotle','Aristotle','Aristotle'),
               stringsAsFactors = F)

corp = corpus(d, docid_field = 'id', text_field='text')
dtm = dfm(corp)

dtm
docvars(dtm)

## ---- message=FALSE, warning=FALSE, echo=TRUE----------------------------
library(RNewsflow)

dtm = rnewsflow_dfm ## copy the demo data
dtm

## ------------------------------------------------------------------------
head(docvars(dtm), 3)

## ------------------------------------------------------------------------
tdd = term.day.dist(dtm)
tail(tdd, 3)

## ------------------------------------------------------------------------
select_terms = tdd$term[tdd$days.entropy.norm <= 0.3]
dtm = dtm[,select_terms]

## ------------------------------------------------------------------------
dtm = quanteda::dfm_tfidf(dtm)

## ----results='hide', message=FALSE, warning=FALSE------------------------
g = newsflow.compare(dtm, hour.window = c(0,36), min.similarity = 0.4)

## ------------------------------------------------------------------------
vertex.sourcetype = V(g)$sourcetype
edge.hourdiff = E(g)$hourdiff

head(vertex.sourcetype)
head(edge.hourdiff)

## ---- fig.show='hold'----------------------------------------------------
v = get.data.frame(g, 'vertices')
e = get.data.frame(g, 'edges')

head(v,3)
head(e,3)    

## ---- fig.width = 7, fig.height = 3--------------------------------------
hist(E(g)$hourdiff, main='Time distance of document pairs', 
     xlab = 'Time difference in hours', breaks = 150)

## ------------------------------------------------------------------------
# set window for all vertices
g = filter.window(g,  hour.window = c(0.1, 36))

# set window for print newspapers
g = filter.window(g,  hour.window = c(6, 36), 
           to.vertices = V(g)$sourcetype == 'Print NP')

## ------------------------------------------------------------------------
show.window(g, to.attribute = 'source')

## ------------------------------------------------------------------------
g_subcomps = decompose.graph(g)

## ---- fig.width = 7, fig.height = 3--------------------------------------
gs = g_subcomps[[55]] # select the second sub-component
document.network.plot(gs)

## ---- fig.width = 7, fig.height = 4--------------------------------------
document.network.plot(gs, source.attribute = 'sourcetype', 
                      dtm=dtm)

## ---- fig.width = 7, fig.height = 3--------------------------------------
gs_onlyfirst = only.first.match(gs)
document.network.plot(gs_onlyfirst)

## ------------------------------------------------------------------------
g.agg = network.aggregate(g, by='source', 
                          edge.attribute='hourdiff', 
                          agg.FUN=median)

e = get.data.frame(g.agg, 'edges')
head(e)

## ------------------------------------------------------------------------
adj.m = get.adjacency(g.agg, attr= 'to.Vprop', sparse = FALSE)
round(adj.m, 2) # round on 2 decimals

## ---- fig.align='center', fig.width=7, fig.height=5----------------------
directed.network.plot(g.agg, weight.var = 'to.Vprop',
                       weight.thres = 0.2)

## ---- fig.align='center', fig.width=7, fig.height=5----------------------
g2 = only.first.match(g)
g2.agg = network.aggregate(g2, by='source', 
                           edge.attribute='hourdiff', 
                           agg.FUN=median)

directed.network.plot(g2.agg, weight.var = 'to.Vprop',
                       weight.thres = 0.2)

## ------------------------------------------------------------------------
V(g)$day = format(as.Date(V(g)$date), '%Y-%m-%d')
agg.perday = network.aggregate(g, 
              by.from='sourcetype', by.to=c('sourcetype', 'day'), 
              edge.attribute='hourdiff', agg.FUN=median, 
              return.df=TRUE)


head(agg.perday[agg.perday$to.sourcetype == 'Online NP',  
                c('from.sourcetype', 'to.sourcetype', 'to.day','to.Vprop')])


## ------------------------------------------------------------------------
agg.perdoc = network.aggregate(g, 
                  by.from='name', by.to='sourcetype', 
                  edge.attribute='weight', agg.FUN=max,
                  return.df=TRUE)
docXsource = xtabs(agg.weight ~ from.name + to.sourcetype, 
                   agg.perdoc, sparse = FALSE)
head(docXsource)

## ----eval=FALSE----------------------------------------------------------
#  library(RNewsflow)
#  library(quanteda)
#  
#  # Prepare DTM
#  dtm = rnewsflow_dfm  ## copy demo data
#  
#  tdd = term.day.dist(dtm)
#  dtm = dtm[,tdd$term[tdd$days.entropy.norm <= 0.3]]
#  
#  dtm = dfm_tfidf(dtm)
#  
#  # Prepare document similarity network
#  g = newsflow.compare(dtm, hour.window = c(-0.1,60),
#                       min.similarity = 0.4)
#  g = filter.window(g, hour.window = c(6, 36),
#                    to.vertices = V(g)$sourcetype == 'Print NP')
#  
#  show.window(g, to.attribute = 'source')
#  
#  g_subcomps = decompose.graph(g)
#  document.network.plot(g_subcomps[[2]], dtm=dtm)
#  
#  g.agg = network.aggregate(g, by='source',
#                            edge.attribute='hourdiff', agg.FUN=median)
#  
#  get.adjacency(g.agg, attr='to.Vprop')
#  directed.network.plot(g.agg, weight.var = 'to.Vprop',
#                        weight.thres=0.2)
#  
#  g2 = only.first.match(g)
#  g2.agg = network.aggregate(g2, by='source',
#                             edge.attribute='hourdiff', agg.FUN=median)
#  
#  get.adjacency(g2.agg, attr='to.Vprop')
#  directed.network.plot(g2.agg, weight.var = 'to.Vprop',
#                        weight.thres=0.2)

