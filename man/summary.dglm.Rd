\name{summary.dglm}
\title{Summarizing double generalized linear model fits}
\alias{print.summary.dglm}
\alias{summary.dglm}
\description{These functions are all \code{methods} for class \code{dglm} or
  \code{summary.glm} objects.
}
\usage{
\method{summary}{dglm}(object, dispersion=NULL, correlation = FALSE, ...)
}
\arguments{
   \item{object}{an object of class \code{"dglm"},
   usually, a result of a call to \code{glm}.}
   \item{dispersion}{the dispersion parameter for the fitting family. 
   By default it is obtained from \code{object}.}
   \item{correlation}{logical; if \code{TRUE}, the correlation matrix of 
   the estimated parameters is returned and printed.}
   \item{...}{further arguments to be passed to \code{summary.glm}}
}

\value{
   \code{summary.dglm} returns an object of class 
   \code{"summary.dglm"}, a list with components
	
   \item{call}{the component from \code{object}}
   \item{terms}{the component from \code{object}}
   \item{family}{the component from \code{object}}
   \item{deviance}{the component from \code{object}}
   \item{aic}{\code{NULL} here}
   \item{constrasts}{(where relevant) the contrasts used. NOT WORKING??}
   \item{df.residual}{the component from \code{object}}
   \item{null.deviance}{the component from \code{object}}
   \item{df.null}{the residual degrees of freedom for the null model.}
   \item{iter}{the component from \code{object}}
   \item{deviance.resid}{the deviance residuals: see \code{residuals.glm}}
   \item{coefficients}{the matrix of coefficients, standard errors, 
   \eqn{z}{z}-values and \eqn{p}{p}-values.  Aliased coefficients are omitted.}
   \item{aliased}{named logical vector showing if the original coefficients are
   aliased.}
   \item{dispersion}{either the supplied argument or the estimated dispersion 
   if the latter in \code{NULL}}
   \item{df}{a 3-vector of the rank of the model and the number of
   residual degrees of freedom, plus number of non-aliased coefficients.}
   \item{cov.unscaled}{the unscaled (\code{dispersion = 1}) estimated covariance
   matrix of the estimated coefficients.}
   \item{cov.scaled}{ditto, scaled by \code{dispersion}}
   \item{correlation}{(only if \code{correlation} is true.)  
   The estimated correlations of the estimated coefficients.}
   \item{dispersion.summary}{the summary of the fitted dispersion model}
   \item{outer.iter}{the number of outer iteration of the alternating iterations}
   \item{m2loglik}{minus twice the log-likelihood of the fitted model}
}


\details{
   For more details, see \code{\link{summary.glm}}.
  
  If more than one of \code{etastart}, \code{start} and \code{mustart}
  is specified, the first in the list will be used.
}

\author{Gordon Smyth,
ported to \R\ by Peter Dunn (\email{pdunn2@usc.edu.au})}

\references{
Smyth, G. K. (1989). Generalized linear models with varying dispersion. 
\emph{J. R. Statist. Soc. B}, \bold{51}, 47--60.

Smyth, G. K., and Verbyla, A. P. (1999). 
Adjusted likelihood methods for modelling dispersion in generalized linear models.
\emph{Environmetrics}, \bold{10}, 696-709.

Verbyla, A. P., and Smyth, G. K. (1998). Double generalized linear models:
approximate residual maximum likelihood and diagnostics. 
Research Report, Department of Statistics, University of Adelaide. 
}

\seealso{\code{\link{dglm.object}}, \code{\link{dglm.control}}, 
\code{\link{anova.dglm}},
\code{\link{summary.glm}} }


\note{
The anova method is questionable when applied to an \code{dglm} object with
\code{method="reml"} (stick to \code{method="ml"}). 
}

\examples{
# Continuing the example from  glm, but this time try
# fitting a Gamma double generalized linear model also.
clotting <- data.frame(
      u = c(5,10,15,20,30,40,60,80,100),
      lot1 = c(118,58,42,35,27,25,21,19,18),
      lot2 = c(69,35,26,21,18,16,13,12,12))
         
# The same example as in  glm: the dispersion is modelled as constant
out <- dglm(lot1 ~ log(u), ~1, data=clotting, family=Gamma)
summary(out)

# Try a double glm 
out2 <- dglm(lot1 ~ log(u), ~u, data=clotting, family=Gamma)

summary(out2)
anova(out2)

# Summarize the mean model as for a glm
summary.glm(out2)
    
# Summarize the dispersion model as for a glm
summary(out2$dispersion.fit)
} 
\keyword{models}

