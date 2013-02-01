\name{getM}
\alias{getM}
\title{Get canonical polynomial chaos expansion size}
\description{
  The getM is a very simple function that determines the size of the canonical polynomial chaos expansion for d random variables up to and including degree p in the expansion. The size is determined from choose(d+p,d) or factorial(d+p)/factorial(d)/factorial(p).   
} 
\usage{getM(d,p)}
\arguments{
  \item{d}{Number of input random variables}
  \item{p}{Order of the polynomial chaos expansion}
}

\value{
  \item{M}{Number of terms in the canonical d-dimensional polynomial chaos expansion up to and including degree p}
}

\references{
  R. Ghanem and P. Spanos, 1991, \emph{Stochastic Finite Elements: A Spectral Approach}. Berlin: Springer.
}

% \seealso{
%  \code{\link{tell.GPCE.lar}}
% }

\examples{
d <- 10
p <- 15
m <- getM(d,p)
}

\author{Jordan Ko}
