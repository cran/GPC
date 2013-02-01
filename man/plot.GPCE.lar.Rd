\name{plot.GPCE.lar}
\alias{plot.GPCE.lar}
\alias{plot,GPCE.lar-method}
\title{
Plot results of adaptive sparse generalized polynomial chaos expansion based on least angle regression
}
\description{
Plot results of adaptive sparse generalized polynomial chaos expansion based on least angle regression
}
\usage{
%Out=PCE(PCSpace="Gaussian",InputDim=3,InputDistrib=c("Gaussian","Gaussian","Gaussian"),Q2tgt=1-10^(-6),pmaxi=12,DesignLength=8,SeedSob=27,CoeffEstiMeth="regression")
\method{plot}{GPCE.lar}(x, ylim = c(0, 1), ...)
}


\arguments{
  \item{x}{A list of arguments returned by the \code{PCE} function}
  \item{ylim}{y-coordinate plotting limits}
  \item{...}{Supplmeentary arguments passed to nodeplot}
}

\seealso{
\code{\link{GPCE.lar}}
}
