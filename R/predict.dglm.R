# predict.dglm <- function(object) {
#   
#   
#   
#   return(list(mean = predict.glm(object, se.fit = TRUE),
#               var = predict.glm(object$dispersion.fit, se.fit = TRUE)))
#   
#   return(list(mean = object$linear.predictors,
#               var = object$dispersion.fit$linear.predictors))
#   
# }
