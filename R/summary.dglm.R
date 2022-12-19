summary.dglm <- function(object, dispersion=NULL, correlation=FALSE, ...)
#  Summarize a double glm
#  Gordon Smyth, 7 Jan 1998
{
  sm <- summary.glm(object, dispersion=dispersion, correlation=correlation, ...)
  sd <- summary.glm(object$dispersion.fit, dispersion=2, correlation=correlation, ...)
  ans <- c(sm, list(dispersion.summary = sd, outer.iter = object$iter, m2loglik = object$m2loglik))
  class(ans) <- c("summary.dglm", "summary.glm")
  ans
}
