\name{ClusProc}
\alias{ClusProc}
\title{This function returns the optimal number of clusters and some descriptive statistics for each object.}
\usage{
  ClusProc(signal, N, prop = 0.9, threshold = 1e-05,
    varSelection = c("PCA", "RAW", "MEAN", "PCA1"),
    thres_sil, adjust = TRUE, thres_MAF = 0.01,
    scale = FALSE)
}
\arguments{
  \item{signal}{for specifying the probe intensity
  measures. It can be a numeric vector, matrix, data.frame
  or list.}

  \item{N}{for applying the given number of clusters to the
  data. N needs to be larger than 1 and if it is 1, error
  will be returned.}

  \item{prop}{a numeric value between 0 and 1. This is used
  to select the number of PC scores, and PC scores
  corresponding to the minimum number of eigenvalues whose
  cumulative proportion is larger than prop will be
  selected. The default is 0.9 and it is meaningful if the
  value of {varSelection} is TRUE. Otherwise it is ignored.
  See {varSelection}.}

  \item{threshold}{Convergence threshold. The iteration
  stops if the absolute difference of loglikelihood between
  successive iterations is less than {threshold}.}

  \item{varSelection}{a factor. For specifying how to
  handle the intensity values. It must take value on 'RAW',
  'PCA', 'PCA1' and 'MEAN'. If the value is 'RAW', then the
  raw intensity value will be used. If it is 'PCA', then
  the first several PCA scores which account for certain
  proportion (this value will be defined in parameter
  {prop}) of all the variance will be used. If the value is
  'PCA1', then the first PCA scores will be used.If the
  value is 'MEAN', the mean of all the probes will be
  used.}

  \item{thres_sil}{??}

  \item{adjust}{??}

  \item{thres_MAF}{??}

  \item{scale}{??}
}
\value{
  It returns object of class 'clust'. 'clust' is a list
  containing at least following componets:
}
\description{
  This function chooses the optimal number of clusters. The
  number of cluster which maximizes the silhouette width is
  selected. no further details yet
}
\author{
  Meiling
}

