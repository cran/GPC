\name{GPCE.quad}
\alias{GPCE.quad}
\title{Solve polynomial chaos expansion with numerical quadrature}
\description{The GPCE.quad function implements a polynomial chaos expansion of a given model or an external model. The strategy for the expansion of the model into a polynomial chaos basis is the Gauss quadrature method where the Gauss quadrature points are used to estimate the integrales corresponding to the coefficients of the expansion. A statistical and a global sensitivity analysis of the model is then carried out.
}

\usage{GPCE.quad(InputDim,PCSpace,InputDistrib,ParamDistrib,
Model,ModelParam,Output,DesignInput,p,ExpPoly,QuadType,QuadPoly,QuadLevel)}
\arguments{
  \item{InputDim}{Number of input random variables}
  \item{PCSpace}{To decide if it should be removed after discussion with Miguel..}
  \item{InputDistrib}{Distribution of the input. Options are \code{Gaussian}, \code{Uniform}, \code{Beta}, \code{Gamma}} 
  \item{ParamDistrib}{Shape parameters of the random variables.}
  \item{Model}{A function defining the model to analyze or NULL if the model is external}
  \item{ModelParam}{Optional parameters for function evaluation.}
  \item{Output}{The results of the model evaluation at the quadrature points.}
  \item{DesignInput}{To decide if it should be removed after discussion with Miguel.}
  \item{p}{The order of the polynomial chaos expansion.}
order) 
  \item{ExpPoly}{The polynomial used in the expansion. Options are \code{Hermite}, \code{Legendre}, \code{Jacobi}, \code{Laguerre}}
  \item{QuadType}{Type of quadrature. Options are \code{SPARSE} or \code{FULL}}
  \item{QuadPoly}{The type of one-dimensional quadrature rule to be used. For \code{SPARSE}, one can use \code{ClenshawCurtis} or \code{Fejer}. For \code{FULL}, one could choose \code{Hermite},\code{Legendre},\code{Jacobi} or \code{Laguerre}}
  \item{QuadLevel}{Level of quadrature used in the approximation.}
}
\value{
  \item{Designs}{A list containing the quadrature size, n, a vector of the n quadrature weights, and a n * d matrix of the quadrature points}
  \item{Output}{Vector of the model output}
  \item{PCEcoeff}{Matrix of the kept sparse polynomial basis. TruncSet_ij is the jth polynomial degree associated to the ith variable}
  \item{Moments}{A list containing the fourth first moments of the output: mean, variance, standard deviation, skweness and kurtosis}
  \item{Sensitivity}{A list containing the sobol sensitivity indices (Values) and the normalized Sobol nominal and total sensitivity indices (S and ST)}
}

\references{
  J. Ko, D. Lucor and P. Sagaut, 2008, \emph{On Sensitivity of Two-Dimensional Spatially Developing Mixing Layer With Respect to Uncertain Inflow Conditions}, Physics of Fluids, 20(7),  07710201-07710220.
  
  J. Ko, 2009, \emph{Applications of the generalized polynomial chaos to the numerical simulationof stochastic shear flows}, Doctoral thesis, University of Paris VI.
  
  J. Ko, D. Lucor and P. Sagaut, 2011, \emph{Effects of base flow uncertainty on Couette flow stability}, Computers and Fluids, 43(1), 82-89.
}

\seealso{
  \code{\link{tell.GPCE.quad}}
}

\examples{

# CASE 1: Validating Legendre PCE with Ishigami function definition
rm(list = setdiff(ls(), lsf.str()))
Model <- function(x){
  a <- 3
  b <- 7
  res <- sin(pi*x[1,]) + a*sin(pi*x[2,])^2 + b*(pi*x[3,])^4*sin(pi*x[1,])
  return(res)
}

# Determine the exact solution for the Ishigami test function 
a <- 3
b <- 7
MeanExact <- a/2
VarExact <- a^2/8 + b*pi^4/5 + b^2*pi^8/18 + 1/2
SobolExact <- c(b*pi^4/5 + b^2*pi^8/50 + 1/2, a^2/8, 0, 0, 8*b^2*pi^8/225, 0, 0)

d = 3
l = 4
ResultObject=GPCE.quad(Model=Model,InputDim=d,PCSpace="Uniform",InputDistrib=rep('Uniform',d),
                       DesignInput=NULL,p=c(l-1),ExpPoly=rep("LEGENDRE",d),QuadType=c("FULL"),
                       QuadPoly=rep("LEGENDRE",d),QuadLevel=c(l),ParamDistrib=NULL,Output=NULL)
names(ResultObject)
cat('Mean and variance normalized absolute errors are', 
  abs(ResultObject$Moments$PCEMean-MeanExact)/MeanExact,
  'and', abs(ResultObject$Moments$PCEVar-VarExact)/VarExact,'\n')
cat('PCE Sobol indices are ',ResultObject$Sensitivity$Values,'\n')
cat('Sobol absolte errors are ',abs(ResultObject$Sensitivity$Values-SobolExact),'\n')

# CASE 2: Validating Hermite PCE with orthogonal Hermite functions
### Model is a R function as a sum of multivariate Hermite polynomials 

Model <- function(x,param){
    d <- param$d
    p <- param$p
    PCETrue <- param$PCETrue
    
    n <- dim(x)[2]
    index <- indexCardinal(d,p)
    PHerm <- hermite.he.polynomials(p, normalized=FALSE)
    y <- rep(0,n)
    
    for (nn in seq(1,n)){
        for (mm in seq(1,getM(d,p))){
            tmp <- 1;
            for (dd in seq(1,d))
            {
            tmp = tmp * unlist(polynomial.values(PHerm[index[dd,mm]+1],x[dd,nn]))
            }     
            y[nn] = y[nn] + PCETrue[mm]*tmp
        }
    }
    return(y)
}

## Problem definition
d = 2;          # random dimension 
l = 3;          # quadrature level
p = l - 1;      # polynomial order of expansion
m = getM(d,p);  # size of polynomial expansion

## Model definition
ModelParam <- NULL
ModelParam$d <- d
ModelParam$p <- p
ModelParam$PCETrue <- sample(seq(1,m),m,replace = FALSE)

## CASE 1: The model is directly evaluated from the GPCE.quad function
ResultObject=GPCE.quad(InputDim=d,PCSpace="Normal",InputDistrib=rep('Gaussian',d),
                       DesignInput=NULL,p=c(p),ExpPoly=rep("HERMITE",d),QuadType=c("FULL"),
                       QuadPoly=rep("HERMITE",d),QuadLevel=c(l),ParamDistrib=NULL,Output=NULL,
                       Model=Model,ModelParam=ModelParam)
cat("The exact PCE coefficients are: \n")
cat(ModelParam$PCETrue,"\n")
cat("The estimated PCE coefficients are: \n")
cat(ResultObject$PCEcoeff,"\n")

## CASE 2: Model is evaluated separately from the GPCE.quad function
# First, the quadrature points are determined from the GPCE.quad function
NoModelResult=GPCE.quad(InputDim=d,PCSpace="Normal",InputDistrib=rep('Gaussian',d),
                       DesignInput=NULL,p=c(p),ExpPoly=rep("HERMITE",d),QuadType=c("FULL"),
                       QuadPoly=rep("HERMITE",d),QuadLevel=c(l),ParamDistrib=NULL,Output=NULL)
cat("The quadrature points can be determined from the Design variable of the output below: \n")
cat(names(NoModelResult),"\n")

# Second, the model is evalauted at the  quadrature points and stored in Output
Output <- Model(NoModelResult$Design$QuadNodes,ModelParam)

# Third, the model output is passed back to GPCE.quad, along with DesignInput and Output
cat("After Design$QuadNodes are evaluated and stored in Output, 
the results is passed back to GPCE.quad:\n")
NoModelResult=GPCE.quad(InputDim=d,PCSpace = "Normal",InputDistrib=rep('Gaussian',d),
                        DesignInput=NoModelResult$Design$QuadNodes,p=c(p),
                        ExpPoly=rep("HERMITE",d),QuadType=c("FULL"),
                        QuadPoly=rep("HERMITE",d),QuadLevel=c(l),
                        ParamDistrib=NULL,Output=Output)
cat("The exact PCE coefficients are:\n")
cat(ModelParam$PCETrue,"\n")
cat("The estimated PCE coefficients are:\n")
cat(NoModelResult$PCEcoeff,"\n")

### Run the algorithm with the lar regression method
# ResultObject=GPCE.lar(Model=Model, PCSpace="Gaussian", InputDim=d, InputDistrib=rep("Gaussian",d))


# CASE 3: Validating Jacobi PCE with orthogonal Jacobi functions 
rm(list = setdiff(ls(), lsf.str()))

# Multivariate Jacobi polynomial test
Model <- function(x){
    a = 2;
    b = 3;
    res <- 2 + 
        3*(2*(a+1)+(a+b+2)*(x[1,]-1))/2 + 
        5*(2*(a+1)+(a+b+2)*(x[2,]-1))/2 + 
        7*(2*(a+1)+(a+b+2)*(x[3,]-1))/2 + 
        11.*(4*(a+1)*(a+2)+4*(a+b+3)*(a+2)*(x[1,]-1)+(a+b+3)*(a+b+4)*(x[1,]-1)^2)/8 + 
        13*(2*(a+1)+(a+b+2)*(x[1,]-1))*(2*(a+1)+(a+b+2)*(x[2,]-1))/4 + 
        17*(2*(a+1)+(a+b+2)*(x[1,]-1))*(2*(a+1)+(a+b+2)*(x[3,]-1))/4 +         
        19.*(4*(a+1)*(a+2)+4*(a+b+3)*(a+2)*(x[2,]-1)+(a+b+3)*(a+b+4)*(x[2,]-1)^2)/8 + 
        23*(2*(a+1)+(a+b+2)*(x[2,]-1))*(2*(a+1)+(a+b+2)*(x[3,]-1))/4 + 
        29.*(4*(a+1)*(a+2)+4*(a+b+3)*(a+2)*(x[3,]-1)+(a+b+3)*(a+b+4)*(x[3,]-1)^2)/8
    return(res)
}

# x is a random variabls following the beta distribution with the pdf
# f(x,a,b) = (1-x)^a(1+x)^b/2^(a+b+1)/Beta(a+1,b+1), where Beta is the beta function 

d = 3
l = 3
ParamDistrib <- NULL
ParamDistrib$alpha <- rep(2,d)  # 1st shape parameters for the random variables
ParamDistrib$beta <- rep(3,d)   # 2nd shape parameters for the random variables

ResultObject=GPCE.quad(Model=Model,InputDim=d,PCSpace="Beta",InputDistrib=rep('Beta',d),
                       DesignInput=NULL,p=c(l-1),ExpPoly=rep("JACOBI",d),QuadType=c("FULL"),
                       QuadPoly=rep("JACOBI",d),QuadLevel=c(l),ParamDistrib,Output=NULL)

cat('Estimated PCE coefficients are', ResultObject$PCEcoeff) 

# CASE 4: Validating Laguerre PCE with orthogonal Laguerre functions 

rm(list = setdiff(ls(), lsf.str()))

# Multivariate Laguerre polynomial test
Model <- function(x){
    a <- 2
    res <- 2 + 
      3*(-x[1,]+a+1) +
      5*(-x[2,]+a+1) + 
      7*(-x[3,]+a+1) +
      11*(x[1,]^2/2-(a+2)*x[1,]+(a+2)*(a+1)/2) + 
      13*(-x[1,]+a+1)*(-x[2,]+a+1) + 
      17*(-x[1,]+a+1)*(-x[3,]+a+1) + 
      19*(x[2,]^2/2-(a+2)*x[2,]+(a+2)*(a+1)/2) + 
      21*(-x[2,]+a+1)*(-x[3,]+a+1) +
      29*(x[3,]^2/2-(a+2)*x[3,]+(a+2)*(a+1)/2) 
    return(res)
}

# x is a random variabls following the gamma distribution with scale parameter theta = 1 and shape 
# parameter k = a + 1 giving the pdf f(x,a) = x^(k-1)*exp(-x)/Gamma(k) = x^a*exp(-x)/Gamma(a+1), 
# where Gamma is the gamma function 

a = 2
d = 3
l = 3

ParamDistrib = NULL
ParamDistrib$alpha = rep(a,d)

ResultObject=GPCE.quad(Model=Model,InputDim=d,PCSpace="Gamma",InputDistrib=rep('Gamma',d),
                       DesignInput=NULL,p=c(l-1),ExpPoly=rep("LAGUERRE",d),QuadType=c("FULL"),
                       QuadPoly=rep("LAGUERRE",d),QuadLevel=c(l),ParamDistrib,Output=NULL)

cat('Estimated PCE coefficients are', ResultObject$PCEcoeff) 
}
\author{Jordan Ko}