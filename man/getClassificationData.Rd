% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markup_processing.R
\name{getClassificationData}
\alias{getClassificationData}
\title{Get Classification data}
\usage{
getClassificationData(
  cm,
  markers,
  prenormalized = F,
  data.gene.id.type = "SYMBOL",
  marker.gene.id.type = "SYMBOL",
  db = NULL,
  verbose = F
)
}
\arguments{
\item{cm}{gene count matrix with cells by columns and genes by rows. May be in raw, TC-normalized or tf-idf-normalized format (in case of tf-idf, `prenormalized` must be set to `T`)}

\item{markers}{path to the file with marker genes or parsed marker list from `parseMarkerFile` function}

\item{prenormalized}{is `cm` in tf-idf-normalized format? Default: FALSE.}
}
\description{
prepare information neccessary for cell type classification
}
