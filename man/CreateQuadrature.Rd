\name{CreateQuadrature}
\alias{CreateQuadrature}
\title{Create a full or sparse numerical quadrature.}

\usage{CreateQuadrature(N,L,QuadPoly,ExpPoly,QuadType,ParamDistrib)}

\description{Creates a full-tensor multivariate quadrature based on one-dimensional quadrature rules or a Smolyak sparse-tensor multivariate quadrature based on nested one-dimensional quadrature rules. 
}

\arguments{
  \item{N}{Number of random variables}
  \item{L}{Level of quadrature used in the approximation.}
  \item{QuadPoly}{The type of one-dimensional quadrature rule to be used. For \code{SPARSE}, one can use \code{ClenshawCurtis} or \code{Fejer}. For \code{FULL}, one could choose \code{Hermite},\code{Legendre},\code{Jacobi} or \code{Laguerre}}  
  \item{ExpPoly}{The polynomial used in the expansion. Options are \code{Hermite}, \code{Legendre}, \code{Jacobi}, \code{Laguerre}}
  \item{QuadType}{Type of quadrature. Options are \code{SPARSE} or \code{FULL}}
  \item{ParamDistrib}{Shape parameters used for the definition of random variables.}
}

\value{
  \item{QuadSize}{Number of quadrature points}
  \item{QuadNodes}{A (QuadSize x n) matrix containing the QuadSize d-dimensional quadrature points.}
  \item{QuadWeights}{A n-tuple vector containing the quadrature weights associated with each quadrature point}
  \item{PolyNodes}{A (m x n) matrix containing the m d-dimensional multivariate polynomial evaluated at the n quadrature points.}
}

\references{J. Ko, 2009, \emph{Applications of the generalized polynomial chaos to the numerical simulationof stochastic shear flows}, Doctoral thesis, University of Paris VI.}

\seealso{
 \code{\link{tell.GPCE.quad}}
}

\examples{
rm(list = setdiff(ls(), lsf.str()))

# random variable definitgeion
d <- 3              # number of random variables
L <- 4              # quadrature level in each dimention. 
                    # could be anisotropic eg c(3,4,5) for full quadrature 

# PCE definition
QuadType <- "FULL"              # type of quadrature
QuadPoly <- rep("LEGENDRE",d)   # polynomial to use
ExpPoly <- rep("LEGENDRE",d)    # polynomial to use
ParamDistrib <- NULL

# QuadType <- "SPARSE"                                  # type of quadrature
# QuadPoly <- 'ClenshawCurtis'                          # polynomial to use
# ExpPoly <- rep("LEGENDRE",d)    # polynomial to use

Quadrature = CreateQuadrature(d,L,QuadPoly,ExpPoly,QuadType,ParamDistrib) # quadrature

}

\author{Jordan Ko}
