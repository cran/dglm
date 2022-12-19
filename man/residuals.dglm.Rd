\name{residuals.dglm}
\alias{residuals.dglm}

\title{Extract Residuals from Double Generalized Linear Model Fit}

\usage{
\method{residuals}{dglm}(object, ...)
}

\arguments{
\item{object}{a \code{"dglm"} object}
\item{...}{any other parameters are passed to \code{residuals.glm}}
}

\value{
Numeric vector of residuals from the mean submodel.
}

\description{
This implements the 'residuals' generic for the dglm object
}

\author{
Robert W. Corty and Gordon Smyth
}
