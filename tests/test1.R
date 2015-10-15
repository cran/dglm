library(dglm)
library(CortyKit)

n <- 400

a <- runif(n)
b <- runif(n)
c <- runif(n)
d <- runif(n)

# betas <- c(runif(n = 3, min = 10, max = 20), 2, 1, 1)
betas <- c(runif(n = 3, min = 10, max = 20), runif(n = 3, min = 1, max = 2))

mu <- betas[1] + a*betas[2] + b*betas[3]
sd <- sqrt(exp(betas[4] + c*betas[5] + d*betas[6]))
y <- mu + rnorm(n = n, sd = sd)

rwc.df <- data.frame(y, a, b, c, d)

mean.formula <- as.formula('y ~ a + b')
var.formula <- as.formula('~ c + d')

l <- TestFunc(mean.formula = mean.formula, var.formula = var.formula, model.df = rwc.df)

par(mfrow = c(1, 2))
plot(mu, l$dglm.mean$fit); abline(0, 1)
legend(x = 'topleft', legend = paste(c('intercept  ', '\t\t\t\t\ta', '\t\t\t\t\tb'), round(betas[1:3], 1)))
plot(sd, exp(l$dglm.var$fit/l$dglm.var$residual.scale)); abline(0, 1)
legend(x = 'topleft', legend = paste(c('intercept  ', '\t\t\t\t\tc', '\t\t\t\t\td'), round(betas[4:6], 2)))

# pdf(file = 'heterosked_fit_all_gray.pdf', width = 15, height = 8)
# plot(x = my.fit$linear.predictors, y = y, 
#      main = 'DGLM illustration',
#      xlab = 'Predicted response',
#      ylab = 'True response',
#      col = gray(0.5, 0.8),
#      pch = 19, cex = 1.5)
# dev.off()
# 
# pdf(file = 'heterosked_fit_dark_high_var.pdf', width = 15, height = 8)
# plot(x = my.fit$linear.predictors, y = y, 
#      main = 'DGLM illustration',
#      xlab = 'Predicted response',
#      ylab = 'True response',
#      col = gray(1 - rescale(sd), 0.8),
#      pch = 19, cex = 1.5)
# dev.off()
# 
# pdf(file = 'heterosked_fit_dark_low_var.pdf', width = 15, height = 8)
# plot(x = my.fit$linear.predictors, y = y, 
#      main = 'DGLM illustration',
#      xlab = 'Predicted response',
#      ylab = 'True response',
#      col = gray(rescale(sd), 0.8),
#      pch = 19, cex = 1.5)
# dev.off()
