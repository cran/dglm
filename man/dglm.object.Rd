\name{dglm-class}
\docType{class}
\alias{dglm-class}
\title{Double Generalized Linear Model - class}

\description{Class of objects returned by fitting double generalized linear models.}

\section{Generation}{
    This class of objects is returned by the \code{\link{dglm}} function 
	 to represent a fitted double generalized linear model. 
	 Class \code{"dglm"} inherits from class \code{"glm"}, 
	 since it consists of two coupled generalized linear models, 
	 one for the mean and one for the dispersion. 
	 Like \code{glm}, 
	 it also inherits from \code{lm}. 
	 The object returned has all the components of a \code{glm} object. 
	 The returned component \code{object$dispersion.fit} is also a 
	 \code{glm} object in its own right, 
	 representing the result of modelling the dispersion.
}
\section{Methods}{
    Objects of this class have methods for the functions 
	 \code{print}, \code{plot}, \code{summary}, \code{anova}, \code{predict}, 
	 \code{fitted}, \code{drop1}, \code{add1}, and \code{step}, amongst others.
	 Specific methods (not shared with \code{glm}) exist for 
	 \code{summary} and \code{anova}. 
}

\section{Structure}{
    A \code{dglm} object consists of a \code{glm} object with the following
	 additional components:
	 
	 \itemize{
		\item{dispersion.fit}{the dispersion submodel: a \code{glm} object 
		representing the fitted model for the dispersions. 
		The responses for this model are the deviance components from the original
		generalized linear model. 
		The prior weights are 1 and the dispersion or scale of this model is 2.}
		
		\item{iter}{this component now represents the number of outer iterations 
		used to fit the coupled mean-dispersion models. 
		At each outer iteration, one IRLS is done for each of the mean and dispersion
		submodels.}
		
		\item{method}{fitting method used: \code{"ml"} if maximum likelihood 
		was used or \code{"reml"} if adjusted profile likelihood was used.}
		
		\item{m2loglik}{minus twice the log-likelihood or 
		adjusted profile likelihood of the fitted model.}
	}
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

\seealso{\code{\link{dglm}}, \code{Digamma family}, \code{Polygamma}}


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

# Examine goodness of fit of dispersion model by plotting residuals
plot(fitted(out2$dispersion.fit),residuals(out2$dispersion.fit)) 
}

\keyword{classes}
\keyword{regression}
