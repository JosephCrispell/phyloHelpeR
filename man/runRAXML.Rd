% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runRAXML.R
\name{runRAXML}
\alias{runRAXML}
\title{Run command line RAxML tool from within R with bootstrapping}
\usage{
runRAXML(
  fastaFile,
  date,
  path,
  nBootstraps = 100,
  nThreads = 6,
  outgroup = NULL,
  model = "GTRCAT"
)
}
\arguments{
\item{fastaFile}{A character string containing the full path to the FASTA file}

\item{date}{A character string corresponding to date e.g. 06-12-19 for 6th December 2019}

\item{nBootstraps}{The number of bootstraps for RAxML to run. Defaults to 100.}

\item{nThreads}{The number of threads for RAxML to use. Defaults to 6.}

\item{outgroup}{A character string identifying a tip in the phylogeny to use as an outgroup. Defaults to \code{NULL}}

\item{model}{A character string identifying the substitution model to used by RAxML}
}
\description{
Runs RAxML from within R with bootstrapping. Note that it requires raxmlHPC to available in the path.
}
\keyword{likelihood}
\keyword{maximum}
\keyword{phylogeny}
\keyword{raxml}
