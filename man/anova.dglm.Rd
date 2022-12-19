\name{anova.dglm}
\title{Analysis of Deviance for Double Generalized Linear Model Fits}
\usage{
\method{anova}{dglm}(object, ...)
}
\alias{anova.dglm}
\description{
  Compute an analysis of deviance table for one or more double generalized
  linear model fits.
}
\arguments{
  \item{object}{objects of class \code{dglm}, typically
    the result of a call to \code{\link{dglm}}.}
  \item{...}{Not used.}
}
\details{
  Specifying a single object gives sequential and adjusted likelihood ratio
  tests for the mean and dispersion model components of the fit.
  The aim is to test overall significance for the mean and dispersion components
  of the double generalized linear model fit.
  The sequential tests (i) set both mean and dispersion models constant,
  add the mean model and (ii) sequentially add the dispersion model.
  The adjusted tests determine whether the mean and dispersion models can be set
  constant separately.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The anova method is questionable when applied to an \code{"dglm"} object with 
  \code{method="reml"} (stick to \code{method="ml"}).
}

\author{Gordon Smyth,
ported to \R\ by Peter Dunn (\email{pdunn2@usc.edu.au})}

\references{
Hastie, T. J. and Pregibon, D. (1992)
\emph{Generalized linear models.}
Chapter 6 of \emph{Statistical Models in S},
edited by J. M. Chambers and T. J. Hastie,
Wadsworth and Brooks/Cole.

Smyth, G. K. (1989). Generalized linear models with varying dispersion.  
\emph{J. R. Statist. Soc. B}, \bold{51}, 47--60.

Smyth, G. K., and Verbyla, A. P. (1999). 
Adjusted likelihood methods for modelling dispersion in generalized linear models.
\emph{Environmetrics}, \bold{10}, 696-709.
\url{http://www.statsci.org/smyth/pubs/Ties98-Preprint.pdf}

Smyth, G. K., and Verbyla, A. P. (1999).
Double generalized linear models: approximate REML and diagnostics.
In \emph{Statistical Modelling: Proceedings of the 14th International Workshop on Statistical Modelling}, Graz, Austria, July 19-23, 1999, H. Friedl, A. Berghold, G. Kauermann (eds.), Technical University, Graz, Austria, pages 66-80.
\url{http://www.statsci.org/smyth/pubs/iwsm99-Preprint.pdf}
}

\seealso{
  \code{\link{dglm}}, \code{\link{anova}}.
}
\examples{
# Continuing the example from  glm, but this time try
# fitting a Gamma double generalized linear model also.
clotting <- data.frame(
      u = c(5,10,15,20,30,40,60,80,100),
      lot1 = c(118,58,42,35,27,25,21,19,18),
      lot2 = c(69,35,26,21,18,16,13,12,12))
         
# The same example as in  glm: the dispersion is modeled as constant
out <- dglm(lot1 ~ log(u), ~1, data=clotting, family=Gamma)
summary(out)

# Try a double glm 
out2 <- dglm(lot1 ~ log(u), ~u, data=clotting, family=Gamma)

summary(out2)
anova(out2)
}
\keyword{models}
\keyword{regression}

