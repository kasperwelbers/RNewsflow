% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.r
\name{only.first.match}
\alias{only.first.match}
\title{Transform document network so that each document only matches the earliest dated matching document}
\usage{
only.first.match(g)
}
\arguments{
\item{g}{A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}}
}
\value{
A network/graph in the \link[igraph]{igraph} class
}
\description{
Transforms the network so that a document only has an edge to the earliest dated document it matches within the specified time window[^duplicate].
}
\details{
If there are multiple earliest dated documents (that is, having the same publication date) then edges to all earliest dated documents are kept.
}
\examples{
data(docnet)

subcomp1 = igraph::decompose.graph(docnet)[[2]]
subcomp2 = only.first.match(subcomp1)

igraph::get.data.frame(subcomp1)
igraph::get.data.frame(subcomp2)

graphics::par(mfrow=c(2,1))
document.network.plot(subcomp1, main='All matches')
document.network.plot(subcomp2, main='Only first match')
graphics::par(mfrow=c(1,1))
}
