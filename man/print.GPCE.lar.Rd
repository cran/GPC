\name{print.GPCE.lar}
\alias{print.GPCE.lar}
\alias{print,GPCE.lar-method}
\title{
Print results of adaptive lar generalized polynomial chaos expansion based on least angle regression
}
\description{
Print results of adaptive lar generalized polynomial chaos expansion based on least angle regression
}
\usage{
%Out=PCE(PCSpace="Gaussian",InputDim=3,InputDistrib=c("Gaussian","Gaussian","Gaussian"),Q2tgt=1-10^(-6),pmaxi=12,DesignLength=8,SeedSob=27,CoeffEstiMeth="regression")
\method{print}{GPCE.lar}(x,...)
}


\arguments{
  \item{x}{A list of arguments returned by the \code{PCE} function}
  \item{...}{Supplementary arguments passed to print()}
}

\seealso{
\code{\link{GPCE.lar}}
}
