\name{predict.GPCE.quad}
\alias{predict.GPCE.quad}
% \alias{predict,GPCE.quad-method}
\title{
predict results of adaptive quad generalized polynomial chaos expansion.
}
\description{
predict results of adaptive quad generalized polynomial chaos expansion.
}
\usage{
\method{predict}{GPCE.quad}(object,n,x,...)
}
\arguments{
  \item{object}{A list of arguments returned by the \code{PCE} function}
  \item{n}{Number of points to evaluate}
  \item{x}{A matrix with the points where the polynomial chaos expansion needs to be evaluated}
  \item{...}{additional parameters}
}
\seealso{
\code{\link{GPCE.quad}}
}
