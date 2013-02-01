\name{print.GPCE.sparse}
\alias{print.GPCE.sparse}
\alias{print,GPCE.sparse-method}
\title{
Print results of adaptive sparse generalized polynomial chaos expansion.
}
\description{
Print results of adaptive sparse generalized polynomial chaos expansion.
}
\usage{
%Out=PCE(PCSpace="Gaussian",InputDim=3,InputDistrib=c("Gaussian","Gaussian","Gaussian"),Q2tgt=1-10^(-6),pmaxi=12,DesignLength=8,SeedSob=27,CoeffEstiMeth="regression")
\method{print}{GPCE.sparse}(x,...)
}


\arguments{
  \item{x}{A list of arguments returned by the \code{PCE} function}
  \item{...}{Supplementary arguments passed to print()}
}

\seealso{
\code{\link{GPCE.sparse}}
}
