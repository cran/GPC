\name{predict.GPCE.lar}
\alias{predict.GPCE.lar}
\alias{predict,GPCE.lar-method}
\title{
predict results of adaptive lar generalized polynomial chaos expansion.
}
\description{
predict results of adaptive lar generalized polynomial chaos expansion.
}
\usage{
\method{predict}{GPCE.lar}(object,x,Selection=1:nrow(object$TruncSet),...)
}
\arguments{
  \item{object}{A list of arguments returned by the \code{PCE} function}
  \item{x}{A matrix with the points where the polynomial chaos expansion needs to be evaluated}
  \item{Selection}{A matrix indicating which polynoms of the expansion will be use for the prediction}
  \item{...}{additional parameters}
}
\seealso{
\code{\link{GPCE.lar}}
}
