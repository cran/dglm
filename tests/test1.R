n <- 500

a <- runif(n)
b <- runif(n)
c <- runif(n)

betas <- runif(6)

y <- betas[1] + a*betas[2] + b*betas[3] + rnorm(exp(betas[4] + b*betas[5] + c*betas[6]))

rwc.df <- data.frame(y, a, b, c)

mean.formula <- as.formula('y ~ a + b')
var.formula <- as.formula('~ b + c')


my.fit <- dglm(formula = mean.formula,
               dformula = var.formula,
               data = rwc.df)

predict(my.fit, se.fit = TRUE)
