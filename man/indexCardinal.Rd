\name{indexCardinal}
\alias{indexCardinal}
\title{Get the order of canonical polynomial chaos expansion}
\description{
  For a d-dimensional polynomial chaos expansion up to order p, there are a total of M polynomials, determine from get(d,p). Each polynomial $phi_m(x)$ is expressed as the product of the polynomial from each random variables, i.e. $phi_m(x) = phi_{m,1}(x_1)phi_{m,2}(x_2)...phi_{m,d}(x_d)$, each with a different polynomial order. We can thus very succintly denote phi_m(x) with a n-tuple vector containing the polynomial order from each dimension and the entire canonical PCE expansion can be express as a (m x n) matrix. 
}
\usage{indexCardinal(d,p,m,index)}
\arguments{
  \item{d}{Number of input random variables}
  \item{p}{Order of the polynomial chaos expansion}
  \item{m}{Pointer to the current n-tuple vector being calculated}
  \item{index}{A dummy variable that stores the polynoial order calculated so far}
}

\value{
  \item{Index}{A (m x d) matrix that donates the polynomial order for each random variables in the canonical construction of PCE}
}

\references{
  R. Ghanem and P. Spanos, 1991, \emph{Stochastic Finite Elements: A Spectral Approach}. Berlin: Springer.
  
  J. Ko, 2009, \emph{Applications of the generalized polynomial chaos to the numerical simulationof stochastic shear flows}, Doctoral thesis, Universit\'e Paris VI.
}

\seealso{
 \code{\link{getM}}
}

\examples{
d <- 3
p <- 5
m <- getM(d,p)
index <- indexCardinal(d,p)
print(t(index))
}

\author{Jordan Ko}
