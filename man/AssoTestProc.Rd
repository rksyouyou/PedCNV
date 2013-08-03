\name{AssoTestProc}
\alias{AssoTestProc}
\title{CNV association testing}
\usage{
  AssoTestProc(signal, ped, envirX, phi, N,
    threshold = 1e-05,
    varSelection = c("PCA", "RAW", "MEAN", "PCA1"),
    H0 = TRUE, thresEM = 0.005, thresAI = 1e-05, itermax)
}
\arguments{
  \item{signal}{TODO}

  \item{ped}{TODO}

  \item{envirX}{TODO}

  \item{phi}{TODO}

  \item{N}{for applying the given number of clusters to the
  data. N needs to be larger than 1 and if it is 1, error
  will be returned.}

  \item{threshold}{Convergence threshold. The iteration
  stops if the absolute difference of loglikelihood between
  successive iterations is less than threshold.\\}

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

  \item{H0}{TODO}

  \item{thresEM}{TODO}

  \item{thresAI}{TODO}

  \item{itermax}{TODO}
}
\value{
  It returns object of class 'asso'. 'asso' is a list
  containing at least following componets:
}
\description{
  This function tests the association of CNV with
  continuous trait based on CNV-linear mixed model.
}
\details{
  no further details yet.
}
\author{
  Meiling
}

