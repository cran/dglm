\name{dglm}
\alias{dglm}
\alias{dglm.constant}
\title{DOUBLE GENERALIZED LINEAR MODELS}
\description{Fits a generalized linear model with a link-linear model for the dispersion as well as for the mean.}
\usage{
dglm(formula=formula(data), dformula = ~ 1, family = gaussian, dlink = "log", 
data = sys.parent(), subset = NULL, weights = NULL, contrasts = NULL, 
method = "ml", mustart = NULL, betastart = NULL, phistart = NULL, 
control = dglm.control(...), ykeep = TRUE, xkeep = FALSE, zkeep = FALSE, ...)

dglm.constant(y,family,weights=1)
}
\arguments{
\item{formula}{a symbolic description of the model to be fit. 
The details of model specification are found in \code{\link{dglm}}.}

\item{dformula}{a formula expression of the form  
\code{~ predictor}, the response being ignored. 
This specifies the linear predictor for modelling the dispersion. 
A term of the form \code{offset(expression)} is allowed.}

\item{family}{a description of the error distribution and link function to
be used in the model. 
See \code{\link{glm}} for more information.}

\item{dlink}{link function for modelling the dispersion. 
Any link function accepted by the \code{quasi}
family is allowed, including \code{power(x)}. See details below.}

\item{data}{an optional data frame containing the variables in the model.
See \code{\link{glm}} for more information.}

\item{subset}{an optional vector specifying a subset of observations to be
used in the fitting process.}

\item{weights}{an optional vector of weights to be used in the fitting
process.}

\item{contrasts}{an optional list. See the \code{contrasts.arg} of
\code{\link{model.matrix.default}}.}

\item{method}{the method used to estimate the dispersion parameters; 
the default is \code{"reml"} for restricted maximum likelihood 
and the alternative is \code{"ml"} for maximum likelihood. 
Upper case and partial matches are allowed.}

\item{mustart}{numeric vector giving starting values for the fitted values 
or expected responses. 
Must be of the same length as the response, 
or of length 1 if a constant starting vector is desired. 
Ignored if \code{betastart} is supplied.}

\item{betastart}{numeric vector giving starting values for the 
regression coefficients in the link-linear model for the mean.}

\item{phistart}{numeric vector giving starting values for the dispersion
parameters.}

\item{control}{a list of iteration and algorithmic constants. 
See \code{\link{dglm.control}} for their names and default values. 
These can also be set as arguments to \code{dglm} itself.}

\item{ykeep}{logical flag: if \code{TRUE}, the vector of responses is returned.}

\item{xkeep}{logical flag: if \code{TRUE}, the \code{model.matrix} for the mean
model is returned.}

\item{zkeep}{logical flag: if \code{TRUE}, the \code{model.matrix} for the
dispersion model is returned.}

\item{...}{further arguments passed to or from other methods.}

\item{y}{numeric response vector}
}
\details{
Write \eqn{\mu_i = \mbox{E}[y_i]}{m_i = E(y_i)} for the expectation of the 
\eqn{i}{i}th response. 
Then \eqn{\mbox{Var}[Y_i] = \phi_i V(\mu_i)}{Var(y_i) = s_iV(m_i)} where \eqn{V}{V}
is the variance function and \eqn{\phi_i}{s_i} is the dispersion of the 
\eqn{i}{i}th response 
(often denoted as the Greek character `phi'). 
We assume the link linear models
\eqn{g(\mu_i) = \mathbf{x}_i^T \mathbf{b}}{g(m_i) = x_i^T b} and
\eqn{h(\phi_i) = \mathbf{z}_i^T \mathbf{z}}{h(s_i) = z_i^T a},
where \eqn{\mathbf{x}_i}{x_i} and \eqn{\mathbf{z}_i}{z_i} are vectors of covariates,
and \eqn{\mathbf{b}}{b} and \eqn{\mathbf{a}}{a} are vectors of regression
cofficients affecting the mean and dispersion respectively. 
The argument \code{dlink} specifies \eqn{h}{h}. 
See \code{\link{family}} for how to specify \eqn{g}{g}. 
The optional arguments \code{mustart}, \code{betastart} and \code{phistart}
specify starting values for \eqn{\mu_i}{m_i}, \eqn{\mathbf{b}}{b}
and \eqn{\phi_i}{s_i} respectively.

The parameters \eqn{\mathbf{b}}{b} are estimated as for an ordinary \acronym{glm}.
The parameters \eqn{\mathbf{a}}{a} are estimated by way of a dual \acronym{glm}
in which the deviance components of the ordinary \acronym{glm} appear as responses.
The estimation procedure alternates between one iteration for the mean submodel 
and one iteration for the dispersion submodel until overall convergence.

The output from \code{dglm}, \code{out} say, consists of two \code{glm} objects
(that for the dispersion submodel is \code{out$dispersion.fit}) with a few more
components for the outer iteration and overall likelihood. 
The \code{summary} and \code{anova} functions have special methods for \code{dglm}
objects. 
Any generic function which has methods for \code{glm}s or \code{lm}s will work on
\code{out}, giving information about the mean submodel. 
Information about the dispersion submodel can be obtained by using
 \code{out$dispersion.fit} as argument rather than out itself. 
In particular \code{drop1(out,scale=1)} gives correct score statistics for 
removing terms from the mean submodel, 
while \code{drop1(out$dispersion.fit,scale=2)} gives correct score 
statistics for removing terms from the dispersion submodel.

The dispersion submodel is treated as a gamma family unless the original 
reponses are gamma, in which case the dispersion submodel is digamma. 
(Note that the digamma and trigamma functions are required to fit a digamma
family.) This is exact if the original \acronym{glm} family is \code{gaussian},
\code{Gamma} or \code{inverse.gaussian}. In other cases it can be 
justified by the saddle-point approximation to the density of the responses. 
The results will therefore be close to exact ML or REML when the dispersions 
are small compared to the means. In all cases the dispersion submodel as prior
weights 1, and has its own dispersion parameter which is 2.
}
\value{an object of class \code{dglm} is returned, 
which inherits from \code{glm} and \code{lm}. 
See \code{\link{dglm.object}} for details.}
}

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
\author{Gordon Smyth, ported to \R\ by Peter Dunn (\email{dunn@usq.edu.au})}

\seealso{\code{\link{dglm.object}}, \code{\link{dglm.control}}, 
\code{Digamma family}, \code{Polygamma}}


\note{
The anova method is questionable when applied to an \code{dglm} object with
\code{method="reml"} (stick to \code{method="ml"}). 
}

\examples{
# Continuing the example from  glm, but this time try
# fitting a Gamma double generalized linear model also.
library(statmod)
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

# Examine goodness of fit of dispersion model by plotting residuals
plot(fitted(out2$dispersion.fit),residuals(out2$dispersion.fit)) 
} 
\keyword{models}
\keyword{regression}
