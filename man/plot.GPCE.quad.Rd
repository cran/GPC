\name{plot.GPCE.quad}
\alias{plot.GPCE.quad}
\alias{plot,GPCE.quad-method}
\title{
Plot results of generalized polynomial chaos expansion based on the quadrature solution
}
\description{
Plot results of generalized polynomial chaos expansion based on the quadrature solution
}
\usage{
%Out=PCE(PCSpace="Gaussian",InputDim=3,InputDistrib=c("Gaussian","Gaussian","Gaussian"),Q2tgt=1-10^(-6),pmaxi=12,DesignLength=8,SeedSob=27,CoeffEstiMeth="regression")
\method{plot}{GPCE.quad}(x,...)
}


\arguments{
  \item{x}{A list of arguments returned by the \code{PCE} function}
  \item{...}{Supplmeentary arguments passed to nodeplot}
}

\seealso{
\code{\link{GPCE.quad}}
}
