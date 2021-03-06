% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_predict.R
\name{cov_predict}
\alias{cov_predict}
\title{Classification Rule for Coverage-based Random Intersection Trees.}
\usage{
cov_predict(rit, testset)
}
\arguments{
\item{rit}{A model produced by \code{cov_RIT}}

\item{testset}{A dataframe containing the instances to classify}
}
\value{
A response vector for the \code{testset} instances
}
\description{
Applies a basic \code{argmax} rule in order to classify new instances.
}
\references{
Ballarini Robin. Random intersection trees for genomic data analysis. Ecole polytechnique de Louvain, Université catholique de Louvain, 2016. Prom. : Dupont, Pierre.
}

