\name{polyNorm}
\alias{polyNorm}
\title{.}
\description{expansion.}

\usage{polyNorm(degree,x,polynomType,alpha,beta)}
\arguments{
  \item{degree}{Number of random inputs}
  \item{x}{A (m x d) matrix that donates the polynomial order for each random variables in the canonical construction of PCE}
  \item{polynomType}{A m-tuple vector containing the PCE coefficients, ordered according to the canonical sequence of the multivariate polynomial expansion}
  \item{alpha}{A m-tuple vector containing the variance of each.}
  \item{beta}{A m-tuple vector containing the variance of each.}
}

\value{
  \item{y}{}
}

\examples{
d <- 3
p <- 2
Index <- indexCardinal(d,p)
}

\author{Jordan Ko}
