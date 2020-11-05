\name{dglm.control}
\alias{dglm.control}
\title{Auxiliary for controlling double glm fitting}
\description{Auxiliary function as user interface for fitting double
generalized linear models. 
Typically only used when calling \code{dglm}.}
\usage{dglm.control(epsilon = 1e-007, maxit = 50, trace = FALSE, ...)
}
\arguments{
\item{epsilon}{positive convergence tolerance epsilon; the iterations
   converge when 
   \eqn{(|L_o - L|)/(|L_o| + 1) > \epsilon}{abs(L.old-L)/(abs(L.old)+1) > epsilon},
   where \eqn{L_o}{L.old} is minus twice the values of log-likelihood
   on the previous iteration, and
   \eqn{L}{L} is minus twice the values of log-likelihood
   on the current.
}

\item{maxit}{integer giving the maximal number of outer iterations
of the alternating iterations.}
\item{trace}{logical indicating if (a small amount of) output should 
be produced for each iteration.}
\item{...}{not currently implemented}
}

\details{
When 'trace' is true, calls to 'cat' produce the output for each
outer iteration. Hence, 'options(digits = *)' can be used to
increase the precision; see the example for \code{\link{glm.control}}.}

\author{Gordon Smyth,
ported to \R\ by Peter Dunn (\email{pdunn2@usc.edu.au})
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

\seealso{\code{\link{dglm.object}}, \code{\link{dglm}}}


\examples{
### A variation on  example(dglm) :
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
oo <- options()
options(digits=12) # See more details in tracing
out2 <- dglm(lot1 ~ log(u), ~u, data=clotting, family=Gamma,
   control=dglm.control(epsilon=0.01, trace=TRUE))
   # With this value of epsilon, convergence should be quicker
   # and the results less reliable (compare to example(dglm) )

summary(out2)
options(oo)
} 
\keyword{optimize}
\keyword{models}

