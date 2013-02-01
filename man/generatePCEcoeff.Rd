\name{generatePCEcoeff}
\alias{generatePCEcoeff}
\title{Determine the PCE coefficients using a numerical quadrature approach}
\description{The inner product used to determine the PCE coefficient is approximated using a numerical quadrature where the integral is approximated with a weighted sum of the function evaluations on some quadrature points.}
\usage{generatePCEcoeff(m,n,y,PhiZn,weights)}
\arguments{
  \item{m}{Number of PCE terms}
  \item{n}{Number of quadrature points}
  \item{y}{A n-tuple vector containing the solution of the function evaluated at the n quadrature points}
  \item{PhiZn}{A (m x n) matrix containing the m d-dimensional multivariate polynomial evaluated at the n quadrature points.}
  \item{weights}{A n-tuple vector containing the quadrature weights associated with each quadrature point}
}

\value{
  \item{ym}{A m-tuple vector containing the PCE coefficients, ordered according to the canonical sequence of the multivariate polynomial expansion}
}

\references{J. Ko, 2009, \emph{Applications of the generalized polynomial chaos to the numerical simulationof stochastic shear flows}, Doctoral thesis, University of Paris VI.}

\seealso{
 \code{\link{getM},\link{CreateQuadrature}}
}

\examples{
rm(list = setdiff(ls(), lsf.str()))

### CASE 1: 
### Ishigami function definition: y= 1 + Phi_1(x1)*Phi_1(x2) + Phi_3(x2)
Model <- function(x){
  a <- 3
  b <- 7
  res <- sin(pi*x[1,]) + a*sin(pi*x[2,])^2 + b*(pi*x[3,])^4*sin(pi*x[1,])
  return(res)
}

# determine the exact solutions for the Ishigami function 
a <- 3
b <- 7
MeanExact <- a/2
VarExact <- a^2/8 + b*pi^4/5 + b^2*pi^8/18 + 1/2
SobolExact <- c(b*pi^4/5 + b^2*pi^8/50 + 1/2, a^2/8, 0, 0, 8*b^2*pi^8/225, 0, 0)

# random variable definitgeion
d <- 3              # number of random variables
L <- 4              # quadrature level in each dimention. 
                    # could be anisotropic eg c(3,4,5) for full quadrature 
P <- L-1            # maximum polynomial expansion (cardinal order)
M <- getM(d,P)      # number of PCE termrs
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


# getting mean and variance
PCEMean = PCE$PCEcoeff[1]
PCEVar = sum(PCE$PCEcoeff[2:M]^2*PCE$PhiIJ[2:M])
# 
# Sobol' sensitivity analysis
Index <- indexCardinal(d,P)
Sobol = getSobol(d,Index,PCE$PCEcoeff,PCE$PhiIJ)

cat('Mean and variance normalized absolute errors are',
     abs(PCEMean-MeanExact)/MeanExact, 'and', abs(PCEVar-VarExact)/VarExact,'\n')
cat('PCE Sobol indices are ',Sobol$Values,'\n')
cat('Sobol absolte errors are ',abs(Sobol$Values-SobolExact),'\n')

}

\author{Jordan Ko}
