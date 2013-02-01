\name{GPCE.sparse}
\alias{GPCE.sparse}
\title{Adaptive sparse generalized polynomial chaos expansion}
\description{
The GPCE.sparse function implements a polynomial chaos expansion of a given model or an external model. The strategy for the expansion of the model into a polynomial chaos basis is the adaptive sparse method where the sampled design is iteratively enriched until the expansion coefficients estimation, by regression, is accurate. Meanwile a relevant selection of the polynomial basis to keep in the expansion is carried out. A statistical and a global sensitivity analysis of the model is then carried out.
}
\usage{
GPCE.sparse(Model = NULL, PCSpace = "Uniform", InputDim = 3, InputDistrib = c(),
                 ParamDistrib = NULL, Q2tgt = 1 - 10^(-6), Eps = 10^(-2) * (1 -
                 Q2tgt), EpsForw = Eps, EpsBack = Eps, EnrichStep = 50, jmax
                 = InputDim, pmaxi = 12, DesignLength = 8, SeedSob =
                 sample(1:1000, 1))
%## S3 method for class 'sobol'
%print(x, ...)
%## S3 method for class 'sobol'
%plot(x, ylim = c(0, 1), ...)
}
\arguments{
  \item{Model}{A function defining the model to analyze or NULL if the model is external}

  \item{PCSpace}{The space where the expansion is achieved.
  Options are Gaussian, Uniform and Physic. Physic use the same
  distributions as the input ones for the expansion}
    
  \item{InputDim}{Dimension of the input}

  \item{InputDistrib}{Distribution of the input. Options are \code{Gaussian}, \code{Uniform}, \code{Beta}, \code{Gamma}}

  \item{ParamDistrib}{Parameters of the input distributions}
  
  \item{Q2tgt}{Fix the accuracy of the expansion fitting.
  By default \code{1-10(^-6)}}

  \item{Eps}{Common epsilon for the selection of the basis.
  By default set to \code{10^(-2)*(1-Q2tgt)}}

  \item{EpsForw}{Epsilon used for the forward selection of
  the basis. By default set to \code{Eps}}

  \item{EpsBack}{Epsilon used for the forward selection of
  the basis. By default set to \code{Eps}}

  \item{EnrichStep}{Number of samples to add to the experimental design. By default set to \code{50}}

  \item{DesignLength}{The length of the input design. By default set to \code{8}}
  
  \item{jmax}{The maximum interaction order between the
  input variables}

  \item{pmaxi}{The maximum degree of the polynomial basis}

  \item{SeedSob}{Seed for the Sobol design generation}
}

\value{
  \item{Designs}{A list containing the Sobol design, the input distributions design, the polynomial chaos design and the design length}
  \item{Output}{Vector of the model output}
  \item{TruncSet}{Matrix of the kept sparse polynomial basis. TruncSet_ij is the jth polynomial degree associated to the ith variable}
  \item{CoeffPCE}{Vector of the expansion coefficients associated to the \code{TruncSet}. CoeffPCE_j is the jth coefficient associaed to the jth polynomial basis.}
  \item{R2}{The R2 PCE approximation error}
  \item{Q2}{The Q2 PCE approximation error}
  \item{Moments}{A list containing the fourth first moments of the output: mean, variance, standard deviation, skweness and kurtosis}
  \item{Sensitivity}{A list containing the sobol sensitivity indices and the sobol total sensitivity indices}
  \item{OutputDistrib}{A list containing a kernel estimation of the output distribution and the associated bandwidth}  
}

\references{
G. Blatman and B. Sudret, 2010, \emph{An adapive algorithm to build up sparse polynomial chaos expansions for stochastic finite element analysis}, Probabilistic Engineering Mechanics, 25, 183--197.
}

\seealso{
\code{\link{tell.GPCE.sparse}}
}

\examples{
### CASE 1: model is a R function. 
### Model definition: y= 1 + Phi_1(x1)*Phi_1(x2) + Phi_3(x2)
Model <- function(x){
  PHerm = hermite.he.polynomials(7, normalized=FALSE)
  y=1+unlist(polynomial.values(PHerm[2],x[,1]))*unlist(polynomial.values(PHerm[2],x[,2]))+
  unlist(polynomial.values(PHerm[4],x[,2]))/sqrt(factorial(3))
  return(y)
}

### Run the algorithm with the sparse regression method
ResultObject=GPCE.sparse(Model=Model, PCSpace="Gaussian", InputDim=3, 
                         InputDistrib=rep("Gaussian",3))
names(ResultObject)
###

### CASE 2: external model (for the example the function Model will be used externaly).

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

