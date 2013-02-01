\name{tell.GPCE.quad}
\alias{tell}
\alias{tell.GPCE.quad}
\alias{tell,GPCE.quad-method}
\title{Generalized polynomial chaos expansion based on numerical quadrature}
\description{The function \code{tell.GPCE.quad} is used for generalized polynomial chaos expansion of external models. When a functional model is given to the function \code{GPCE.quad}, this latter automatically call the function \code{tell.GPCE.quad}. 

When in the function \code{GPCE.quad} the \code{Model} option is \code{NULL} then this latter returns a list of \code{Designs} and a list of Arguments, \code{Args}, the ones entered for the previous \code{GPCE.quad} run. 

Then the \code{Designs}, the \code{Args} and the vector of manually calculated outputs need to be given to the function \code{tell.GPCE.quad} to run the algorithm.
}
\usage{
\method{tell}{GPCE.quad}(x,Output,...)
}
\arguments{
The same as the \code{\link{GPCE.quad}} function
  \item{x}{ResultObjectDesign2EvalArgs: This is the argument defined by the user.}
  \item{Output}{The vector of the manually calculated model outputs}
  \item{...}{additional parameters}
}
\value{The same as the \code{\link{GPCE.quad}} function}
\references{
  J. Ko, D. Lucor and P. Sagaut, 2008, \emph{On Sensitivity of Two-Dimensional Spatially Developing Mixing Layer With Respect to Uncertain Inflow Conditions}, Physics of Fluids, 20(7),  07710201-07710220.
  
  J. Ko, 2009, \emph{Applications of the generalized polynomial chaos to the numerical simulationof stochastic shear flows}, Doctoral thesis, University of Paris VI.
  
  J. Ko, D. Lucor and P. Sagaut, 2011, \emph{Effects of base flow uncertainty on Couette flow stability}, Computers and Fluids, 43(1), 82-89.
}
\seealso{
\code{\link{GPCE.quad}}
}
\examples{
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
cat("After Design$QuadNodes are evalauted and stored in Output, 
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
}

\author{Jordan Ko}
