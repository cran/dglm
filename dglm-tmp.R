dglm <- function(formula = formula(data),
		dformula = ~1,
		family = gaussian(),
		dlink = "log",
		data = sys.parent(),
		subset = NULL,
		weights = NULL,
		contrasts = NULL,
		method = "ml",
		mustart = NULL,
		betastart = NULL,
		phistart = NULL,
		control = dglm.control(...),
		ykeep = T,
		xkeep = F,
		zkeep = F,
		...)
{
#
#  Double generalized linear models
#	Manual page www.maths.uq.edu.au/~gks/s/dglm.html
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	8 Dec 1997. Last revised 22 Oct 1999.
#
#  Set up mean submodel: 
#               y   response
#               X   design matrix
#               w   prior weights
#          offset   offset in linear predictor
#          family   response family
#         mustart   starting values for mu (optional)
#       betastart   starting values for coefficients (optional)
#
#	Save call for future reference
	scall <- match.call()
#
#	Evaluate the model frame
	mnames <- c("", "formula", "data", "weights", "subset")
	cnames <- names(scall)
	cnames <- cnames[match(mnames,cnames,0)]
	mcall <- scall[cnames]
	mcall[[1]] <- as.name("model.frame")
	mframe <- eval(mcall, sys.parent())
#
#	Now extract the glm components
	y <- model.extract(mframe, response)
	if(is.null(dim(y))) N <- length(y) else N <- dim(y)[1]
	mterms <- attr(mframe, "terms")
	X <- model.matrix(mterms, mframe, contrasts)
	w <- model.extract(mframe, weights)
	if(is.null(w)) w <- rep(1,N)
	offset <- model.extract(mframe, offset)
	# family <- as.family(family)
#
#  Set up dispersion submodel:
#               Z   design matrix 
#         doffset   offset in linear predictor
#         dfamily   family for unit deviances
#        phistart   starting values for phi (optional)
#
#	Setup dformula with y on left hand side
#	(to ensure that mean and dispersion submodels have same length)
	mcall$formula <- formula
	mcall$formula[3] <- switch(match(length(dformula),c(0,2,3)),
		1,dformula[2],dformula[3])
#
#	Evaluate the model frame and extract components
	mframe <- eval(mcall, sys.parent())
	dterms <- attr(mframe, "terms")
	Z <- model.matrix(dterms, mframe, contrasts)
	doffset <- model.extract(mframe, offset)
#
#	Parse the dispersion link and evaluate as list
	name.dlink <- substitute(dlink)
cat("name.dlink= ", name.dlink,"\n")
cat("dlink= ", dlink,"\n")
	if(is.name(name.dlink))
		if(is.character(dlink))  # link is character variable
			name.dlink <- dlink
		else                     # link is name without quotes
			dlink <- name.dlink <- as.character(name.dlink)
	else
		if(is.call(name.dlink))  # power link
			name.dlink <- deparse(name.dlink)
	#if(match(name.dlink, dimnames(glm.links)[[2]], F)) # If not match...
	#	dlink <- glm.links[, name.dlink]
	if(!is.null(dlink$name))
		name.dlink <- dlink$name
	cat("name.dlink=",name.dlink,"\n")
#
#	Construct the dispersion variance function
	Digamma <- family$family=="Gamma" || (family$family=="Tweedie" && family$variance$p==2)
print(family)
print(family$family)
names(family)
cat("Digamma=",Digamma,"\n")
	if (Digamma) {
		linkinv <- make.link(name.dlink)$linkinv
		linkfun <- make.link(name.dlink)$linkfun
		mu.eta <- make.link(name.dlink)$mu.eta
		valid.eta <- make.link(name.dlink)$valid.eta
		init <- expression({
		    if (any(y <= 0)) 
      		  stop(paste("Non-positive values not", "allowed for the Digamma family"))
		    n <- rep.int(1, length(y) )
		    mustart <- y
			}
		)
		dfamily <- structure(list(family="Digamma", variance=varfun.digamma, 
				dev.resids=function(y,mu){sign(y - mu) * sqrt(deviance.digamma(y,mu))},
				aic=NA, link=name.dlink, linkfun=linkfun, linkinv=linkinv, 
				mu.eta=mu.eta, initialize=init, 
				validmu=function(mu) { all(mu>0) }, valideta=valid.eta) )

#		dfamily <- make.family("Digamma", dlink, var, name.dlink, "2*[1/theta + trigamma(-theta)]")
	} else
		dfamily <- make.family("Gamma", dlink, glm.variances[, "mu^2"], name.dlink, "Square: mu^2")

# structure(list(family = "Tweedie", variance = variance, dev.resids = dev.resids, 
#        aic = aic, link = paste("mu^", as.character(lambda), 
#            sep = ""), linkfun = linkfun, linkinv = linkinv, 
#        mu.eta = mu.eta, initialize = initialize, validmu = validmu, 
#        valideta = valideta), class = "family")




#
#	Remember if log-link for use with initial values
	dlink <- as.character(dfamily$link)
	logdlink <- dlink=="log"
#
#  Match method (ml or reml)
#
	if(!is.null(scall$method)) {
		name.method <- substitute(method)
		if(!is.character(name.method))
			name.method <- deparse(name.method)
		list.methods <- c("ml","reml","ML","REML","Ml","Reml")
		i.method <- pmatch(method,list.methods,nomatch=0)
		if(!i.method) stop("Method must be ml or reml")
		method <- switch(i.method,"ml","reml","ml","reml","ml","reml")
	}
	reml <- method=="reml"
#
#  Starting values.  If explicit starting values are not supplied,
#  regress link(y) on X and dlink(d) on Z by ordinary least squares.
# 

## This line causes probelms for some reason
#	eval(family$initialize)
	nobs <- length(y)
	mu <- y+(1/6)
	bob <- family
	print(bob$initialize)
	cat(bob$family,"\n")

	if(!is.null(betastart)) {
		eta <- X %*% betastart
		mu <- family$inverse(eta+offset)}
	else if(!is.null(mustart)) {
		mu <- mustart
		eta <- family$link(mu)-offset}
	else {
		eta <- lm.fit(X,family$linkfun(mu)-offset,singular.ok=T)$fitted.values
		mu <- family$inverse(eta+offset)}
	d <- family$deviance(mu, y, w, residuals=T)^2
	if(!is.null(phistart)) {
		phi <- phistart
		deta <- dfamily$link(phi)-doffset}
	else {
		deta <- lm.fit(Z,dfamily$link(d+(d==0)/6)-doffset,singular.ok=T)$fitted.values
		if(logdlink) deta <- deta + 1.27036
		phi <- dfamily$inverse(deta+offset)}
#
#  First iteration of mean submodel
#  (necessary to get adjusted profile likelihood correct)
#
	zm <- eta + (y - mu) * family$deriv(mu)
	wm <- eval(family$weight)/phi
	mfit <- lm.wfit(X,zm,wm,method="qr",singular.ok=T,qr=reml)
	eta <- mfit$fitted.values
	mu <- family$inverse(eta+offset)
	d <- family$deviance(mu, y, w, residuals=T)^2
#
#  Initial (minus twice log) likelihood or adjusted profile likelihood
#
	const <- glm.constant(y,family,w)
	if(Digamma)
		h <- 2*(lgamma(w/phi)+(1+log(phi/w))*w/phi)
	else
		h <- log(phi/w)
	m2loglik <- const + sum(h + d/phi)
	if(reml)
		m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))
	m2loglikold <- m2loglik+1
#
#  Estimate model by alternate iterations
#
	epsilon <- control$epsilon
	maxit <- control$maxit
	trace <- control$trace
	iter <- 0
	while ( abs(m2loglikold-m2loglik)/(abs(m2loglikold)+1) > epsilon && iter < maxit )
	{
#
#		dispersion submodel
		hdot <- dfamily$deriv(phi)
		if(Digamma) {
			delta <- 2*w*( log(w/phi) - digamma(w/phi) )
			u <- 2*w^2*( trigamma(w/phi) - phi/w)
			fdot <- phi^2 / u * hdot
			}
		else {
			delta <- phi
			u <- phi^2
			fdot <- hdot
			}
		wd <- 1 / (fdot^2 * u)
#		wd <- eval(dfamily$weight,local=list(mu=phi,w=rep(1,N),family=dfamily))}
		if(reml) {
			h <- hat(mfit$qr)
			delta <- delta - phi*h
			wd <- wd - 2*(h/hdot^2/phi^2) + h^2
		}
		if(any(wd<0)) wd[wd<0] <- 0
		zd <- deta + (d - delta) * fdot
		dfit <- lm.wfit(Z,zd,wd,method="qr",singular.ok=T)
		deta <- dfit$fitted.values
		phi <- dfamily$inverse(deta+doffset)
#
#		mean submodel
		zm <- eta + (y - mu) * family$deriv(mu)
		wm <- eval(family$weight)/phi
		mfit <- lm.wfit(X,zm,wm,method="qr",singular.ok=T,qr=reml)
		eta <- mfit$fitted.values
		mu <- family$inverse(eta+offset)
		d <- family$deviance(mu, y, w, residuals=T)^2
#
#		overall likelihood
		m2loglikold <- m2loglik
		if(Digamma)
			h <- 2*(lgamma(w/phi)+(1+log(phi/w))*w/phi)
		else
			h <- log(phi/w)
		m2loglik <- const + sum(h + d/phi)
		if(reml)
			m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))
		iter <- iter+1
		if(trace)
			cat("DGLM iteration ", iter, ": -2*log-likelihood = ",
			format(round(m2loglik, 4)), " \n", sep = "")
	}
#
#  Output for mean model:
#  Exactly as for glm.object.  As for lm.object except that
#  linear.predictors and prior.weights are new components, and fitted.values
#  has a new definition.
#
	mfit$family <- family$family
	mfit$linear.predictors <- mfit$fitted.values+offset
	mfit$fitted.values <- mu
	if(!is.null(scall$weights)) mfit$prior.weights <- w
	mfit$terms <- mterms
	mfit$call <- scall
	mfit$deviance <- sum(d/phi)
	mfit$null.deviance <- glm.null(X, y, w/phi, offset, family)
	if(length(mfit$null.deviance)>1) mfit$null.deviance <- mfit$null.deviance$null.deviance
	if(ykeep) mfit$y <- y
	if(xkeep) mfit$x <- X
	mfit$formula <- as.vector(attr(mterms, "formula"))
#
#  Output for dispersion model:
#  As for glm.object except that prior.weights are not relevant.  Is
#  nested in one output component.
#
	dfit$family <- dfamily$family
	dfit$linear.predictors <- dfit$fitted.values+doffset
	dfit$fitted.values <- phi
	dfit$terms <- dterms
	scall$formula <- scall$dformula
	scall$dformula <- NULL
	scall$family <- call(dfamily$family["name"],link=as.name(dlink))
	dfit$call <- scall
	dfit$deviance <- dfamily$deviance(phi, d, w=rep(1/2,N))
	dfit$null.deviance <- glm.null(Z, d, rep(1/2,N), doffset, dfamily)
	if(length(dfit$null.deviance)>1) dfit$null.deviance <- dfit$null.deviance$null.deviance
	if(ykeep) dfit$y <- d
	if(zkeep) dfit$z <- Z
	dfit$formula <- as.vector(attr(dterms, "formula"))
	dfit$iter <- iter
	class(dfit) <- c("glm","lm")
	out <- c(mfit, list(dispersion.fit = dfit, iter=iter, method=method, m2loglik=m2loglik))
	class(out) <- c("dglm","glm","lm")
	out
}

dglm.control <- function(epsilon = 1e-007, maxit = 50, trace = F, ...)
{
#  Control iteration for dglm
#  As for glm.control but with different defaults
#  GKS  9 Jan 98
#
	if(epsilon <= 0) {
		warning("the value of epsilon supplied is zero or negative;
the default value of 1e-7 was used instead")
		epsilon <- 1e-007}
	if(maxit < 1) {
		warning("the value of maxit supplied is zero or negative;
the default value of 50 was used instead")
		maxit <- 50}
	list(epsilon = epsilon, maxit = maxit, trace = as.logical(trace)[1])
}

glm.constant <- function(y,family,w=1)
{
#  Constant term appearing in glm log-likelihood
#  Used by dglm
#  GKS  6 Jan 98, 4 Jul 98, 23 Sep 99.
#
#  Exact cases (in binomial case, exact for phi near 1)
#
	const <- switch(family$family[1],
		"Gaussian" = length(y)*log(2*pi),
		"Poisson" = 2*sum(y-y*ifelse(y>0,log(y),0)+lgamma(y+1)),
		"Gamma" = 2*sum(log(y)),
		"Inverse Gaussian" = sum(log(2*pi*y^3)),
		"Tweedie" = switch(match(family$variance[[2]],c(0,1,2,3),nomatch=0),
			length(y)*log(2*pi),
			2*sum(y-y*ifelse(y>0,log(y),0)+lgamma(y+1)),
			2*sum(log(y)),
			sum(log(2*pi*y^3)) ),
		"Binomial" = -2*sum(lgamma(w+1)-lgamma(w*y+1)-lgamma(w*(1-y)+1)+w*(y*ifelse(y>0,log(y),0)+(1-y)*ifelse(1-y>0,log(1-y),0)))+sum(log(w))
	)
#
#  Saddle-point approximation
#
	if(is.null(const)) {
		V <- family$variance(y)
		if(any(V==0)) V[V==0] <- family$variance(y[V==0]+1/6)
		const <- sum(log(2*pi*V))
		if(length(V)==1 && length(y)>1) const <- length(y)*const
	}
	const
}

summary.dglm <- function(object, correlation = F)
{
#  Summarize a double glm
#  GKS  7 Jan 98
	sm <- summary.glm(object,dispersion=1,correlation=correlation)
	sd <- summary.glm(object$dispersion.fit,dispersion=2,correlation=correlation)
	structure( c(sm, list(dispersion.summary=sd, outer.iter=object$iter,
		m2loglik=object$m2loglik)), class=c("summary.dglm","summary.glm"))
}

print.summary.dglm <- function(x, digits = NULL, quote = T, prefix = "", residuals = F)
{
#  Print summary of double glm
#  GKS  7 Jan 98
#
	xd <- x$dispersion.summary
	x$dispersion.summary <- NULL
	if(is.null(digits))
		digits <- options()$digits
	else {
		old.digits <- options(digits = digits)
		on.exit(options(old.digits))
	}
	cat("\nCall: ")
	print(x$call)
#
#  Mean submodel
#
#	cat("\nMEAN MODEL")
	nas <- x$nas
	coef <- x$coef
	correl <- x$correl
	if(any(nas)) {
		nc <- length(nas)
		cnames <- names(nas)
		coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
		coef1[!nas,  ] <- coef
		coef <- coef1
		if(!is.null(correl)) {
			correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, cnames)
				)
			correl1[!nas, !nas] <- correl
			correl <- correl1
		}
	}
	dresid <- x$deviance.resid
	n <- length(dresid)
	rdf <- x$df[2]
	if(residuals) {
		if(rdf > 5) {
			cat("Deviance Residuals:\n")
			rq <- quantile(as.vector(dresid))
			names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
			print(rq, digits = digits)
		}
		else if(rdf > 0) {
			cat("Deviance Residuals:\n")
			print(dresid, digits = digits)
		}
	}
	if(any(nas))
		cat("\nMean Coefficients: (", sum(nas), 
			" not defined because of singularities)\n", sep = "")
	else cat("\nMean Coefficients:\n")
	print(coef, digits = digits)
	cat(paste("(Dispersion Parameters for", names(x$dispersion), 
		"family estimated as below", ")\n"))
	int <- attr(x$terms, "intercept")
	if(is.null(int))
		int <- 1
	cat("\n    Scaled Null Deviance:", format(round(x$null.deviance, digits)), "on", n - 
		int, "degrees of freedom\n")
	cat("Scaled Residual Deviance:", format(round(x$deviance, digits)), "on", round(
		rdf, digits), "degrees of freedom\n")
#       cat("\nNumber of Fisher Scoring Iterations:", format(trunc(x$iter)), "\n")
	if(!is.null(correl)) {
		p <- dim(correl)[2]
		if(p > 1) {
			cat("\nCorrelation of Coefficients:\n")
			ll <- lower.tri(correl)
			correl[ll] <- format(round(correl[ll], digits))
			correl[!ll] <- ""
			print(correl[-1,  - p, drop = F], quote = F, digits = digits)
		}
	}
#
#  Dispersion submodel
#
	nas <- xd$nas
	coef <- xd$coef
	correl <- xd$correl
	if(any(nas)) {
		nc <- length(nas)
		cnames <- names(nas)
		coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
		coef1[!nas,  ] <- coef
		coef <- coef1
		if(!is.null(correl)) {
			correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, cnames)
				)
			correl1[!nas, !nas] <- correl
			correl <- correl1
		}
	}
	dresid <- xd$deviance.resid
	n <- length(dresid)
	rdf <- xd$df[2]
	if(residuals) {
		if(rdf > 5) {
			cat("Deviance Residuals:\n")
			rq <- quantile(as.vector(dresid))
			names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
			print(rq, digits = digits)
		}
		else if(rdf > 0) {
			cat("Deviance Residuals:\n")
			print(dresid, digits = digits)
		}
	}
	if(any(nas))
		cat("\nDispersion Coefficients: (", sum(nas), 
			" not defined because of singularities)\n", sep = "")
	else cat("\nDispersion Coefficients:\n")
	print(coef, digits = digits)
	cat(paste("(Dispersion Parameter for", names(xd$dispersion), 
		"family taken to be", format(round(xd$dispersion, digits)), ")\n"))
	int <- attr(xd$terms, "intercept")
	if(is.null(int))
		int <- 1
	cat("\n    Scaled Null Deviance:", format(round(xd$null.deviance, digits)), "on", n - 
		int, "degrees of freedom\n")
	cat("Scaled Residual Deviance:", format(round(xd$deviance, digits)), "on", round(
		rdf, digits), "degrees of freedom\n")
#	cat("\nNumber of Fisher Scoring Iterations:", format(trunc(xd$iter)), "\n")
	if(!is.null(correl)) {
		p <- dim(correl)[2]
		if(p > 1) {
			cat("\nCorrelation of Coefficients:\n")
			ll <- lower.tri(correl)
			correl[ll] <- format(round(correl[ll], digits))
			correl[!ll] <- ""
			print(correl[-1,  - p, drop = F], quote = F, digits = digits)
		}
	}
#
#  Overall iteration
#
	cat("\nMinus Twice the Log-Likelihood:", format(round(x$m2loglik, digits)), "\n")
	cat("Number of Alternating Iterations:", format(trunc(x$outer.iter)), "\n")
	invisible(NULL)
}

anova.dglm <- function(object)
{
#  ANOVA for double glm (likelihood ratio tests for mean and dispersion models)
#  GKS  16 Jan 98
#
	response <- as.character(formula(object)[2])
	mterms <- object$terms
	dterms <- object$dispersion.fit$terms
	df <- c(length(mterms),length(dterms))
#
#  Extract null mean model formula
#
	if(df[1]>0) {
		o <- attr(mterms, "offset")
		i <- attr(mterms, "intercept")
		factor.labels <- dimnames(attr(mterms,"factors"))[[1]]
		fm <- paste(response,"~")
		if(!is.null(o)) {
			fm <- paste(fm,factor.labels[o])
			if(i == 0)
				fm <- paste(fm,"-1")
			}
		else {
			if(i == 0)
				fm <- paste(fm,"-1")
			else
				fm <- paste(fm,"1")
			}
		fm <- parse(text=fm)
		mode(fm) <- "call"
		df[1] <- object$rank - i
	}
#
#  Extract null dispersion model formula
#
	if(df[2]>0) {
		o <- attr(dterms, "offset")
		i <- attr(dterms, "intercept")
		factor.labels <- dimnames(attr(dterms,"factors"))[[1]]
		if(!is.null(o)) {
			fd <- paste("~",factor.labels[o])
			if(i == 0)
				fd <- paste(fd,"-1")
			}
		else {
			if(i == 0)
				fd <- "~ -1"
			else
				fd <- "~ 1"
			}
		fd <- parse(text=fd)
		mode(fd) <- "call"
		df[2] <- object$dispersion.fit$rank - i
	}
#
#  Fit null models and store likelihoods
# 
	lik <- rep(object$m2loglik,4)
	names(lik) <- c("Null","Mean","Full","Disp")
	ncall <- object$call
	if(df[2]>0) {	
		ncall["dformula"] <- fd
		lik["Mean"] <- eval(ncall)$m2loglik
		if(df[1]>0) {
			ncall["formula"] <- fm
			lik["Null"] <- eval(ncall)$m2loglik
			ncall["dformula"] <- object$call["dformula"]
			ncall["formula"] <- fm
			lik["Disp"] <- eval(ncall)$m2loglik
			}
		else {
			lik["Null"] <- lik["Mean"]
			lik["Disp"] <- lik["Full"]
			}
		}
	else {
		lik["Mean"] <- lik["Full"]
		if(df[1]>0) {
			ncall["formula"] <- fm
			lik["Null"] <- eval(ncall)$m2loglik
			lik["Disp"] <- lik["Null"]
			}
		else
			lik["Disp"] <- lik["Null"] <- lik["Full"]
		}
	seqdev <- c(lik["Null"] - lik["Mean"], lik["Mean"] - lik["Full"])
	adjdev <- c(lik["Disp"] - lik["Full"], lik["Mean"] - lik["Full"])
#
#  Output table
#
	heading <- c("Analysis of Deviance Table",
		paste("\n",object$family["name"]," double generalized linear model",sep=""),
		paste("\nResponse: ", response,"\n", sep = "") )
	aod <- data.frame(row.names = c("Mean model","Dispersion model"),
		DF = df,
		"Seq Chisq" = seqdev,
		P = 1-pchisq(seqdev,ifelse(df>0,df,NA)),
		"Adj Chisq"=adjdev,
		P = 1-pchisq(adjdev,ifelse(df>0,df,NA)))
	return(as.anova(aod, heading))
}
