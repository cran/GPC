\name{GPC-package}
\alias{GPC-package}
\alias{GPC}
\docType{package}
\title{
Generalized Polynomial Chaos
}
\description{
Generalized polynomial chaos expansion of a model taking independent random input distribution. Three different methods are available for the expansion: the Gauss quadrature method, the adaptive sparse basis method and the least angle regression based adaptive sparse basis method. A statistical and a global sensitivity analysis of the model is then carried out.
}

\details{
The PCE package implements a polynomial chaos expansion of a model. The main option is the choice between three strategy to expand the model into a polynomial chaos basis. 
\itemize{
\item The Gauss quadrature method where the the Gauss quadrature points are used to estimate the integrales giving the coefficients of the expansion. 
\item The adaptive sparse method where the sampled design is iteratively enriched until the expansion coefficients estimation, by regression, is accurate. Meanwile a relevant selection of the polynomial basis to keep in the expansion is carried out.
\item The adaptive sparse method based on the least angle regression.
}
}

\section{Model managing}{
The PCE package has been designed to work with either models written in R than external models
such as heavy computational codes. This is achieved with the input argument model present in all
functions of this package.
The argument model is expected to be either a funtion or NULL.
%either a funtion or a predictor (i.e. an object with a predict function such as lm).
If model = m where m is a function, it will be invoked once by y <- m(X).
%. If model = m where m is a predictor, it will be invoked once by y <- predict(m, X).
X is the design of experiments, i.e. a data.frame with p columns (the input factors) and n lines
(each, an experiment), and y is the vector of length n of the model responses.
The model in invoked once for the whole design of experiment.
The argument model can be left to NULL. This is refered to as the decoupled approach and used with
external computational codes that rarely run on the statistician's computer. See decoupling.
}

\references{
G. Blatman and B. Sudret, 2010, \emph{An adapive algorithm to build up sparse polynomial chaos expansions for stochastic finite element analysis}, Probabilistic Engineering Mechanics, 25, 183--197.

G. Blatman and B. Sudret, 2011, \emph{Adapive sparse polynomial chaos expansion based on least angle regression}, Journal of Computational Physics, 230, 2345--2367.

J. Ko, D. Lucor and P. Sagaut, 2008, \emph{On Sensitivity of Two-Dimensional Spatially Developing Mixing Layer With Respect to Uncertain Inflow Conditions}, Physics of Fluids, 20(7),  07710201-07710220.

M. Munoz Zuniga ang S. Kucherenko, 2013, \emph{Metamodelling with independent and dependent inputs}, Computer Physics Communications,xx ,xx.
}

\keyword{ polynomial chaos} 
\keyword{ sensitivity } 
\keyword{ Gauss quadrature } 
\keyword{ sparse basis} 
\keyword{ least angle regression} 

\seealso{
The model management of the PCE package has some similarity with the one used in the package "sensitivity"
%\code{\link[<sensitivity>:<sensitivity>-package]} 
%{\link[<sensitivity>]}
}
