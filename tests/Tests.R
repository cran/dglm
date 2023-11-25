library(dglm)
options(warnPartialMatchArgs=TRUE,warnPartialMatchAttr=TRUE,warnPartialMatchDollar=TRUE,width=120)

set.seed(0); u <- runif(100)

n <- 100
y <- rnorm(n)
x <- rnorm(n)
z <- rnorm(n)

fit <- dglm(y~x,~z,method="ml")
summary(fit)
fit <- dglm(y~x,~z,method="reml")
summary(fit)
anova(fit)

y2 <- y^2
fit <- dglm(y2~x,~z,family=Gamma(link="log"),method="ml")
summary(fit)
fit <- dglm(y2~x,~z,family=Gamma(link="log"),method="reml")
summary(fit)

# Continuing the example from glm, but this time try
# fitting a Gamma double generalized linear model also.
clotting <- data.frame(
      u = c(5,10,15,20,30,40,60,80,100),
      lot1 = c(118,58,42,35,27,25,21,19,18),
      lot2 = c(69,35,26,21,18,16,13,12,12))
         
# The same example as in  glm: the dispersion is modelled as constant
# However, dglm uses  ml  not  reml,  so the results are slightly different:
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
