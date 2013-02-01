\name{getSobol}
\alias{getSobol}
\title{Determine Sobol indices from PCE solutions}
\description{Computes the Sobol sensitivity indices from the coefficient of the polynomial chaos expansion.}

\usage{getSobol(d,Index,Coeff,PhiIJ)}
\arguments{
  \item{d}{Number of random inputs}
  \item{Index}{A (m x d) matrix that donates the polynomial order for each random variables in the canonical construction of PCE}
  \item{Coeff}{A m-tuple vector containing the PCE coefficients, ordered according to the canonical sequence of the multivariate polynomial expansion}
  \item{PhiIJ}{A m-tuple vector containing the variance of each.}
}

\value{
  \item{Names}{The }
  \item{Values}{}
}

\references{
  I. M. Sobol', 1993, \emph{Sensitivity estimate for non-linear mathematical models}. Mathematical modelling and computational experiments 1, 407-414.
  
  J. Ko, 2009, \emph{Applications of the generalized polynomial chaos to the numerical simulation of stochastic shear flows}, Doctoral thesis, University of Paris VI.
}

\seealso{
 \code{\link{getM},\link{CreateQuadrature}}
}

\examples{
### Ishigami function definition
Model <- function(x){
  a <- 3
  b <- 7
  res <- sin(pi*x[1,]) + a*sin(pi*x[2,])^2 + b*(pi*x[3,])^4*sin(pi*x[1,])
  return(res)
}

# Find exact solution for the Ishigami test function 
a <- 3
b <- 7
MeanExact <- a/2
VarExact <- a^2/8 + b*pi^4/5 + b^2*pi^8/18 + 1/2
SobolExact <- c(b*pi^4/5 + b^2*pi^8/50 + 1/2, a^2/8, 0, 0, 8*b^2*pi^8/225, 0, 0)

# random variable definition
d <- 3              # number of random variables
L <- 4              # quadrature level in each dimension. 
                    # could be anisotropic eg c(3,4,5) for full quadrature 
P <- L-1            # maximum polynomial expansion (cardinal order)
M <- getM(d,P)      # number of PCE terms
ParamDistrib <- NULL
# ParamDistrib<- list(beta=rep(0.0,d),alpha=rep(0.0,d))

# PCE definition
QuadType <- "FULL"              # type of quadrature
QuadPoly <- rep("LEGENDRE",d)   # polynomial to use
ExpPoly <- rep("LEGENDRE",d)    # polynomial to use

# QuadType <- "SPARSE"                                  # type of quadrature
# QuadPoly <- 'ClenshawCurtis'                          # polynomial to use
# ExpPoly <- rep("LEGENDRE",d)    # polynomial to use

Quadrature = CreateQuadrature(d,L,QuadPoly,ExpPoly,QuadType,ParamDistrib) # quadrature
 
# function sampling
y <- Model(Quadrature$QuadNodes) # Ishigami function d=3 

# generate PCE coefficients
PCE = generatePCEcoeff(M,Quadrature$QuadSize,y,Quadrature$PolyNodes,Quadrature$QuadWeights) # PCE

# Sobol' sensitivity analysis
Index <- indexCardinal(d,P)
Sobol = getSobol(d,Index,PCE$PCEcoeff,PCE$PhiIJ)

cat('PCE Sobol indices are ',Sobol$Values,'\n')
cat('Sobol absolte errors are ',abs(Sobol$Values-SobolExact),'\n')
}

\author{Jordan Ko}
