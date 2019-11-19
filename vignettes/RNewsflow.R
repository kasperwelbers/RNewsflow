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
tdd = term_day_dist(dtm)
tail(tdd, 3)

## ------------------------------------------------------------------------
select_terms = tdd$term[tdd$days.entropy.norm <= 0.3]
dtm = dtm[,select_terms]

## ------------------------------------------------------------------------
dtm = quanteda::dfm_tfidf(dtm)

## ----results='hide', message=FALSE, warning=FALSE------------------------
g = newsflow_compare(dtm, date_var='date',
                     hour_window = c(0,36), 
                     min_similarity = 0.4)

## ------------------------------------------------------------------------
vertex_sourcetype = V(g)$sourcetype
edge_hourdiff = E(g)$hourdiff

head(vertex_sourcetype)
head(edge_hourdiff)

## ---- fig.show='hold'----------------------------------------------------
v = as_data_frame(g, 'vertices')
e = as_data_frame(g, 'edges')

head(v[,c('name','date','source','sourcetype')],3)
head(e,3)    

## ---- fig.width = 7, fig.height = 3--------------------------------------
hist(E(g)$hourdiff, main='Time distance of document pairs', 
     xlab = 'Time difference in hours', breaks = 150, right=F)

## ------------------------------------------------------------------------
# set window for all vertices
g = filter_window(g,  hour_window = c(0.1, 36))

# set window for print newspapers
g = filter_window(g,  hour_window = c(6, 36), 
           to_vertices = V(g)$sourcetype == 'Print NP')

## ------------------------------------------------------------------------
show_window(g, to_attribute = 'source')

## ------------------------------------------------------------------------
g_subcomps = decompose(g)

## ---- fig.width = 7, fig.height = 3--------------------------------------
gs = g_subcomps[[55]] # select the second sub-component
document_network_plot(gs)

## ---- fig.width = 7, fig.height = 4--------------------------------------
document_network_plot(gs, source_attribute = 'sourcetype', 
                      dtm=dtm)

## ---- fig.width = 7, fig.height = 3--------------------------------------
gs_onlyfirst = only_first_match(gs)
document_network_plot(gs_onlyfirst)

## ------------------------------------------------------------------------
g_agg = network_aggregate(g, by='source', 
                          edge_attribute='hourdiff', 
                          agg_FUN=median)

e = as_data_frame(g_agg, 'edges')
head(e)

## ------------------------------------------------------------------------
adj_m = as_adjacency_matrix(g_agg, attr= 'to.Vprop', sparse = FALSE)
round(adj_m, 2) # round on 2 decimals

## ---- fig.align='center', fig.width=7, fig.height=5----------------------
directed_network_plot(g_agg, weight_var = 'to.Vprop',
                       weight_thres = 0.2)

## ---- fig.align='center', fig.width=7, fig.height=5----------------------
g2 = only_first_match(g)
g2_agg = network_aggregate(g2, by='source', 
                           edge_attribute='hourdiff', 
                           agg_FUN=median)

directed_network_plot(g2_agg, weight_var = 'to.Vprop',
                       weight_thres = 0.2)

## ------------------------------------------------------------------------
V(g)$day = format(as.Date(V(g)$date), '%Y-%m-%d')
agg_perday = network_aggregate(g, 
              by_from='sourcetype', by_to=c('sourcetype', 'day'), 
              edge_attribute='hourdiff', agg_FUN=median, 
              return_df=TRUE)

head(agg_perday[agg_perday$to.sourcetype == 'Online NP',  
                c('from.sourcetype', 'to.sourcetype', 'to.day','to.Vprop')])

## ------------------------------------------------------------------------
agg_perdoc = network_aggregate(g, 
                  by_from='name', by_to='sourcetype', 
                  edge_attribute='weight', agg_FUN=max,
                  return_df=TRUE)
docXsource = xtabs(agg.weight ~ from.name + to.sourcetype, 
                   agg_perdoc, sparse = FALSE)
head(docXsource)

## ----eval=FALSE----------------------------------------------------------
#  library(RNewsflow)
#  library(quanteda)
#  
#  # Prepare DTM
#  dtm = rnewsflow_dfm  ## copy demo data
#  
#  tdd = term_day_dist(dtm)
#  dtm = dtm[,tdd$term[tdd$days.entropy.norm <= 0.3]]
#  
#  dtm = dfm_tfidf(dtm)
#  
#  # Prepare document similarity network
#  g = newsflow_compare(dtm, hour_window = c(-0.1,60),
#                       min_similarity = 0.4)
#  g = filter_window(g, hour_window = c(6, 36),
#                    to_vertices = V(g)$sourcetype == 'Print NP')
#  show_window(g, to_attribute = 'source')
#  
#  g_subcomps = decompose(g)
#  document_network_plot(g_subcomps[[55]], dtm=dtm)
#  
#  g_agg = network_aggregate(g, by='source',
#                            edge_attribute='hourdiff', agg_FUN=median)
#  
#  as_adjacency_matrix(g_agg, attr='to.Vprop')
#  directed_network_plot(g_agg, weight_var = 'to.Vprop',
#                        weight_thres=0.2)
#  
#  g2 = only_first_match(g)
#  g2_agg = network_aggregate(g2, by='source',
#                             edge_attribute='hourdiff', agg_FUN=median)
#  
#  as_adjacency_matrix(g2_agg, attr='to.Vprop')
#  directed_network_plot(g2_agg, weight_var = 'to.Vprop',
#                        weight_thres=0.2)

