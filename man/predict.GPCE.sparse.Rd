\name{predict.GPCE.sparse}
\alias{predict.GPCE.sparse}
\alias{predict,GPCE.sparse-method}
\title{
predict results of adaptive sparse generalized polynomial chaos expansion.
}
\description{
predict results of adaptive sparse generalized polynomial chaos expansion.
}
\usage{
\method{predict}{GPCE.sparse}(object,x,Selection=1:nrow(object$TruncSet),...)
}
\arguments{
  \item{object}{A list of arguments returned by the \code{PCE} function}
  \item{x}{A matrix with the points where the polynomial chaos expansion needs to be evaluated}
  \item{Selection}{A matrix indicating which polynoms of the expansion will be use for the prediction}
  \item{...}{additional parameters}
}
\seealso{
\code{\link{GPCE.sparse}}
}

