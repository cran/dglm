\name{summary.dglm}
\title{Summarize Double Generalized Linear Model Fit}
\alias{summary.dglm}
\alias{print.summary.dglm}
\description{Summarize objects of class \code{"dglm"}.}
\usage{
\method{summary}{dglm}(object, dispersion=NULL, correlation = FALSE, ...)
}
\arguments{
   \item{object}{an object of class \code{"dglm"}.}
   \item{dispersion}{the dispersion parameter for the fitting family. 
   By default it is obtained from \code{object}.}
   \item{correlation}{logical; if \code{TRUE}, the correlation matrix of 
   the estimated parameters is returned and printed.}
   \item{...}{further arguments to be passed to \code{summary.glm}}
}

\value{
   An object of class \code{"summary.dglm"}, which is a list with the following components:
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
ported to R by Peter Dunn (\email{pdunn2@usc.edu.au})}

\references{
Smyth, G. K. (1989). Generalized linear models with varying dispersion. 
\emph{J. R. Statist. Soc. B}, \bold{51}, 47--60.
\doi{10.1111/j.2517-6161.1989.tb01747.x}

Smyth, G. K., and Verbyla, A. P. (1999). 
Adjusted likelihood methods for modelling dispersion in generalized linear models.
\emph{Environmetrics}, \bold{10}, 696-709.
\doi{10.1002/(SICI)1099-095X(199911/12)10:6<695::AID-ENV385>3.0.CO;2-M}
\url{https://gksmyth.github.io/pubs/Ties98-Preprint.pdf}

Smyth, G. K., and Verbyla, A. P. (1999).
Double generalized linear models: approximate REML and diagnostics.
In \emph{Statistical Modelling: Proceedings of the 14th International Workshop on Statistical Modelling}, Graz, Austria, July 19-23, 1999, H. Friedl, A. Berghold, G. Kauermann (eds.), Technical University, Graz, Austria, pages 66-80.
\url{https://gksmyth.github.io/pubs/iwsm99-Preprint.pdf}
}

\seealso{
\code{\link{dglm}},
\code{\link{dglm-class}},
\code{\link{summary.glm}}
}

\note{
The anova method is questionable when applied to an \code{dglm} object created with
\code{method="reml"} (stick to \code{method="ml"}). 
}

\keyword{models}

