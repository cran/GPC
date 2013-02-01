\name{tell.GPCE.lar}
\alias{tell.GPCE.lar}
\alias{tell,GPCE.lar-method}
\title{
Adaptive sparse generalized polynomial chaos expansion based on least angle regression
}
\description{
The function \code{tell.GPCE.lar} is used for generalized polynomial chaos expansion of external models. When a functional model is given to the function \code{GPCE.lar}, this latter automatically call the function \code{tell.GPCE.lar}. When in the function \code{GPCE.lar} the \code{Model} option is \code{NULL} then this latter returns a list of \code{Designs} and a list of Arguments, \code{Args}, the ones entered for the previous \code{GPCE.lar} run. Then the \code{Designs}, the \code{Args} and the vector of manually calculated outputs need to be given to the function \code{tell.GPCE.lar} to run the algorithm.
}
\usage{
\method{tell}{GPCE.lar}(x,Output,...)
}
\arguments{
  \item{x}{ResultObjectDesign2EvalArgs: A list of arguments returned by the \code{PCE} function}
  \item{Output}{The vector of the manually calculated model outputs}
  \item{...}{additional parameters}
}
\value{
The same as the \code{\link{GPCE.lar}} function
}
\references{
  G. Blatman and B. Sudret, 2011, \emph{Adapive sparse polynomial chaos expansion based on least angle regression}, Journal of Computational Physics, 230, 2345--2367.}
\seealso{
\code{\link{GPCE.lar}}
}
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
ResultObject=GPCE.lar(PCSpace="Gaussian",InputDim=3,InputDistrib=rep("Gaussian",3))
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
### In the later case the user calculate the output 
### externally and give them to the tell function
### with the previous ResultObject for further calculation. 
### See GPCE.sparse documentation for an example.
}
\author{
  Munoz Zuniga Miguel
}
