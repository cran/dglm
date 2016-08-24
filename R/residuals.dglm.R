#' @title residuals.dglm function
#' @name residuals.dglm
#' @author Robert W. Corty
#'
#' @param object The fitted dglm object whose residuals we want
#' @param ... additional parameters to residuals.glm
#' 
#' @description This implements the 'residuals' generic for the dglm object
#'
#' @return the residuals from the mean portion of the dglm object
#' @export
#'
residuals.dglm <- function(object, ...) {
  return(stats::residuals.glm(object = object, ...))
}
