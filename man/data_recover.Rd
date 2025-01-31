% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_recovery.R
\name{data_recover}
\alias{data_recover}
\alias{data_recover.DGCPCA}
\title{Recover the mapped data (i.e. all discrete features are numerical) using first several PC scores}
\usage{
data_recover(x, npc_used, ...)

\method{data_recover}{DGCPCA}(x, npc_used, ...)
}
\arguments{
\item{x}{An \code{DGCPCA} object.}

\item{npc_used}{The number of PCs used for data recovery, no larger than \code{x$npc}.}
}
\description{
Recover the mapped data (i.e. all discrete features are numerical) using first several PC scores
}
\section{Methods (by class)}{
\itemize{
\item \code{DGCPCA}: S3 method for \code{DGCPCA}
}}

