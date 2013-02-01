\name{tell.GPCE.sparse}
\alias{tell.GPCE.sparse}
\alias{tell,GPCE.sparse-method}
\title{
Adaptive sparse generalized polynomial chaos expansion
}
\description{
The function \code{tell.GPCE.sparse} is used for generalized polynomial chaos expansion of external models. When a functional model is given to the function \code{GPCE.sparse}, this latter automatically call the function \code{tell.GPCE.sparse}. When in the function \code{GPCE.sparse} the \code{Model} option is \code{NULL} then this latter returns a list of \code{Designs} and a list of Arguments, \code{Args}, the ones entered for previous the \code{GPCE.sparse} run. Then the \code{Designs}, the \code{Args} and the vector of manually calculated outputs need to be given to the function \code{tell.GPCE.sparse} to run the algorithm.
}
\usage{
%Out=PCE(PCSpace="Gaussian",InputDim=3,InputDistrib=c("Gaussian","Gaussian","Gaussian"),Q2tgt=1-10^(-6),pmaxi=12,DesignLength=8,SeedSob=27,CoeffEstiMeth="regression")
\method{tell}{GPCE.sparse}(x,Output,...)
%tell.GPCE.sparse(ResultObjectDesign2EvalArgs,Output)
}
\arguments{
  \item{x}{ResultObjectDesign2EvalArgs: a list of arguments returned by the \code{GPCE.sparse} function}
  \item{Output}{the vector of the manually calculated model outputs}
  \item{...}{additional parameters}
}
\value{The same as the \code{\link{GPCE.sparse}} function}
\references{
G. Blatman and B. Sudret, 2010, \emph{An adapive algorithm to build up sparse polynomial chaos expansions for stochastic finite element analysis}, Probabilistic Engineering Mechanics, 25, 183--197.
}
\seealso{\code{\link{GPCE.sparse}}}
\examples{

### External model (for the example the function Model defined below will be used externaly)
### Model definition: y= 1 + Phi_1(x1)*Phi_1(x2) + Phi_3(x2)
Model <- function(x){
  PHerm = hermite.he.polynomials(5, normalized=FALSE)
  y=1+unlist(polynomial.values(PHerm[2],x[,1]))*unlist(polynomial.values(PHerm[2],x[,2]))+
  unlist(polynomial.values(PHerm[4],x[,2]))/sqrt(factorial(3))
  return(y)
}

### initialized Output
Output=c()

### Get a first design
ResultObject=GPCE.sparse(PCSpace="Gaussian",InputDim=3,InputDistrib=rep("Gaussian",3))
names(ResultObject)

### Calculate the model output for the given design and concatenate the model output results
### into the output vector 
Output=c(Output,Model(ResultObject$Design2Eval))

### Give the design and the calculated ouput to the tell function
ResultObject=tell(ResultObject,Output)
names(ResultObject)

### If the expansion has been calculated the function tell return the full expansion 
### paramaters, the moments analysis, the sensitivity analysis and the output distribution
### If not the function tell() return an enriched design. 
Output=c(Output,Model(ResultObject$Design2Eval))

### Give the design and the calculated ouput to the tell function
ResultObject=tell(ResultObject,Output)
names(ResultObject)
###
}

\author{
  Munoz Zuniga Miguel
}
