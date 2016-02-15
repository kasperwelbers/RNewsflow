## ---- include=FALSE------------------------------------------------------
options(digits=3)
library(knitr)

## ---- message=F, warning=F, echo=T---------------------------------------
doc1 = 'Socrates is human'
doc2 = 'Humans are mortal'
doc3 = 'Therefore, Socrates is mortal'
dtm = RTextTools::create_matrix(
  textColumns = c(doc1,doc2,doc3), 
  minWordLength = 1, removeStopwords = F)

rownames(dtm) = paste('Document', 1:nrow(dtm))
as.matrix(dtm)

## ---- message=F, warning=F, echo=T---------------------------------------
library(RNewsflow)

data(dtm)
as.matrix(dtm[1:3,1:5])

## ------------------------------------------------------------------------
data(meta)
head(meta,3)

## ------------------------------------------------------------------------
tdd = term.day.dist(dtm, meta)
tdd[sample(1:nrow(tdd), 4),] # show 4 random rows

## ------------------------------------------------------------------------
select_terms = tdd$term[tdd$days.entropy.norm <= 0.3]
dtm = dtm[,select_terms]

## ------------------------------------------------------------------------
dtm = weightTfIdf(dtm)

## ----results='hide', message=FALSE, warning=FALSE------------------------
g = newsflow.compare(dtm, meta,
                             hour.window = c(0,36), 
                             min.similarity = 0.4)

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
gs = g_subcomps[[2]] # select the second sub-component
document.network.plot(gs)

## ---- fig.width = 7, fig.height = 4--------------------------------------
document.network.plot(gs, source.attribute = 'sourcetype', 
                      dtm=dtm)

## ---- fig.width = 7, fig.height = 3--------------------------------------
gs_onlyfirst = only.first.match(gs)
document.network.plot(gs_onlyfirst)

