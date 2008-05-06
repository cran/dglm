dglm <- function(formula = formula(data),
      dformula = ~1,
      family = gaussian,
      dlink = "log",
      data = sys.parent(),
      subset = NULL,
      weights = NULL,
      contrasts = NULL,
      method = "ml",
      mustart = NULL,
      betastart = NULL,
      etastart = NULL,
      phistart = NULL,
      control = dglm.control(...),
      ykeep = TRUE,
      xkeep = FALSE,
      zkeep = FALSE,
      verbose=FALSE,
      ...)
#
#   Double generalized linear models
#   Gordon Smyth, Walter and Eliza Hall Institute of Medical Research
#   S-Plus version created 8 Dec 1997, last revised 22 Oct 1999.
#
#   Ported to R by Peter Dunn, 22 August 2005
#
{
#  Set up mean submodel: 
#               y   response
#               X   design matrix
#               w   prior weights
#          offset   offset in linear predictor
#          family   response family
#         mustart   starting values for mu (optional)
#       betastart   starting values for coefficients (optional)
#
#   Save call for future reference
   call <- match.call()

#   Get family for mean model
   if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
   if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
#
#   Evaluate the model frame
   mnames <- c("", "formula", "data", "weights", "subset")
   cnames <- names(call)

   cnames <- cnames[match(mnames,cnames,0)]
   mcall <- call[cnames]
   mcall[[1]] <- as.name("model.frame")
   mframe <- eval(mcall, sys.parent())
   mf <- match.call(expand.dots=FALSE)

#   Now extract the glm components
   y <- model.response(mframe, "numeric")
   
   if(is.null(dim(y))) {
      N <- length(y) 
   } else {
      N <- dim(y)[1]
   }
   nobs <- N # Needed for some of the  eval  calls
   mterms <- attr(mframe, "terms")
   X <- model.matrix(mterms, mframe, contrasts)
   
   weights <- model.weights(mframe)
   if(is.null(weights)) weights <- rep(1,N)
   if ( is.null(weights) ) weights <- rep(1, N )
   if (!is.null(weights) && any(weights < 0)) {
        stop("negative weights not allowed")
   }
   
   offset <- model.offset(mframe)
   if ( is.null(offset) ) offset <- rep(0, N )
   if (!is.null(offset) && length(offset) != NROW(y)) {
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(offset), NROW(y)), domain = NA)
   }

#
#  Set up dispersion submodel:
#               Z   design matrix 
#         doffset   offset in linear predictor
#         dfamily   family for unit deviances
#        phistart   starting values for phi (optional)
#
#   Setup dformula with y on left hand side
#   (to ensure that mean and dispersion submodels have same length)
   mcall$formula <- formula
   mcall$formula[3] <- switch(match(length(dformula),c(0,2,3)),
      1,dformula[2],dformula[3])
#
#   Evaluate the model frame and extract components
   mframe <<- eval(mcall, sys.parent())
   dterms <- attr(mframe, "terms")
   Z <- model.matrix(dterms, mframe, contrasts)
   doffset <- model.extract(mframe, offset)
   if ( is.null(doffset) ) doffset <- rep(0, N )
#
#   Parse the dispersion link and evaluate as list
   name.dlink <- substitute(dlink)
   if(is.name(name.dlink)) {
      if(is.character(dlink)) { # link is character variable
         name.dlink <- dlink
      } else {                  # link is name without quotes
         dlink <- name.dlink <- as.character(name.dlink)
      }
   } else {
      if(is.call(name.dlink))  # power link
         name.dlink <- deparse(name.dlink)
   }
#   if(!is.null(dlink$name))   name.dlink <- dlink$name ### WHAT DOES THIS DO??
#
   if ( verbose ) cat("dlink =",dlink,"\n")
#   Construct the dispersion variance function
   if ( family$family=="Tweedie") {
      tweedie.p <- call$family$var.power
   }

   Digamma <- family$family=="Gamma" || (family$family=="Tweedie" && tweedie.p==2)

   if (Digamma) {

      linkinv <- make.link(name.dlink)$linkinv
      linkfun <- make.link(name.dlink)$linkfun
      mu.eta <- make.link(name.dlink)$mu.eta
      valid.eta <- make.link(name.dlink)$valid.eta
      init <- expression({
         if (any(y <= 0)){
            print(y)
            print( any(y<=0) )
           stop("non-positive values not allowed for the DM gamma family")
           }
            n <- rep.int(1, nobs)
            mustart <- y
            })

      
      dfamily <- structure(list(family="Digamma", variance=varfun.digamma, 
            dev.resids=function(y,mu,wt){wt * unitdeviance.digamma(y,mu)}, # gaussian()$dev.resids are resids^2
            aic=function(y, n, mu, wt, dev) NA,  
            link=name.dlink, linkfun=linkfun, linkinv=linkinv, 
            mu.eta=mu.eta, initialize=init, 
            validmu=function(mu) { all(mu>0) }, valideta=valid.eta) )

   } else { # Not digamma family
      eval(substitute(dfamily <- Gamma(link=lk), list(lk=name.dlink ) ))
   }
#
#   Remember if log-link for use with initial values
#
   dlink <- as.character(dfamily$link)
   logdlink <- dlink=="log"
#
#  Match method (ml or reml)
#
   if(!is.null(call$method)) {
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

   if ( is.null(mustart) ) { 
      etastart <- NULL
      eval(family$initialize)  
      mu <- mustart
      mustart <- NULL
   }

   if(!is.null(betastart)) {
      eta <- X %*% betastart
      mu <- family$linkinv(eta+offset)
   } else {
      if(!is.null(mustart)) {
         mu <- mustart
         eta <- family$linkfun(mu)-offset
      } else {

         eta <- lm.fit(X,family$linkfun(mu)-offset,singular.ok=TRUE)$fitted.values
         # Recall:  fitted values are on the linear predictor scale
         mu <- family$linkinv(eta+offset)
      }
   }

#   Careful:   Called dev.resid, but appear to be the deviance residuals squared
   d <- family$dev.resids(y, mu, weights)

    if(!is.null(phistart)) {
       phi <- phistart
       deta <- dfamily$linkfun(phi)-doffset
   } else {
       deta <- lm.fit(Z,dfamily$linkfun(d+(d==0)/6)-doffset,singular.ok=TRUE)$fitted.values
       if(logdlink) deta <- deta + 1.27036
       phi <- dfamily$linkinv(deta+offset)
    }

   zm <- eta + (y - mu) / family$mu.eta(eta)
   wm <- eval(family$variance(mu))*weights/phi

   mfit <- lm.wfit(X,zm,wm,method="qr",singular.ok=TRUE) # ,qr=reml)
   eta <- mfit$fitted.values
   mu <- family$linkinv(eta+offset)
   d <- family$dev.resids(y, mu, weights)

#  Initial (minus twice log) likelihood or adjusted profile likelihood
   const <- dglm.constant(y,family,weights)

   if(Digamma) {
      h <- 2*(lgamma(weights/phi)+(1+log(phi/weights))*weights/phi)
   } else {
      h <- log(phi/weights)
   }
   m2loglik <- const + sum(h + d/phi)

   if(reml)
      m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))
      
   m2loglikold <- m2loglik+1
#  Estimate model by alternate iterations

   epsilon <- control$epsilon
   maxit <- control$maxit
   trace <- control$trace
   iter <- 0
   
        while ( abs(m2loglikold-m2loglik)/(abs(m2loglikold)+1) > epsilon && iter < maxit )
   {
################################
	if ( verbose ) cat("Iteration",iter+1,":\n")
#      dispersion submodel
      hdot <- 1/dfamily$mu.eta( deta )
      if(Digamma) {
         delta <- 2*weights*( log(weights/phi) - digamma(weights/phi) )
         u <- 2*weights^2*( trigamma(weights/phi) - phi/weights)
         fdot <- phi^2 / u * hdot  # p 50
      } else { # In normal and iG cases, eg, the dispersion sub-model is gamma
          delta <- phi
         u <- phi^2 # variance function for disp. model; u(delta)=delta^2 is gamma
         fdot <- hdot
      }
      wd <- 1 / (fdot^2 * u)   # ie Smyth, p 50.  We don't include the factor of 2,
               # as that is the disp. parameter, which enters via
               # the  dispersion=2  argument for the summary.

      if(reml) {
         h <- hat(mfit$qr)
         delta <- delta - phi*h
         wd <- wd - 2*(h/hdot^2/phi^2) + h^2
      }
      
      if(any(wd<0)) {
         wd[wd<0] <- 0
      }
      zd <- deta + (d - delta) * fdot

      dfit <- lm.wfit(Z, zd, wd, method="qr", singular.ok=TRUE)
      deta <- dfit$fitted.values
   
      phi <- dfamily$linkinv(deta+doffset)
      if ( verbose ) cat(" Dispersion sub-model done\n")

################################
#      mean submodel
      zm <- eta + (y - mu) / family$mu.eta(eta)
      fam.wt <- expression( weights * family$variance(mu) ) 
      wm <- eval( fam.wt )/phi
      mfit <- lm.wfit(X,zm,wm,method="qr",singular.ok=TRUE)
      eta <- mfit$fitted.values
      mu <- family$linkinv(eta+offset)
      d <- family$dev.resids(y, mu, weights)

#      overall likelihood
      m2loglikold <- m2loglik
      if(Digamma) {
         h <- 2*(lgamma(weights/phi)+(1+log(phi/weights))*weights/phi)
      } else {
         h <- log(phi/weights)
      }
      m2loglik <- const + sum(h + d/phi)
      if(reml) {
         m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))
      }
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
   if ( verbose ) cat("Model fitted\n")
   mfit$formula <- call$formula
   mfit$call <- call
   mfit$family <- family

   mfit$linear.predictors <- mfit$fitted.values+offset
   mfit$fitted.values <- mu
   mfit$prior.weights <- weights
   mfit$terms <- mterms
   mfit$contrasts <- attr(X, "contrasts")
   intercept <- attr(mterms, "intercept")
   mfit$df.null <- N - sum(weights == 0) - as.integer(intercept)
   mfit$call <- call
   mfit$deviance <- sum(d/phi)
   mfit$aic <- NA
   mfit$null.deviance <- glm.fit(x=X, y=y, weights=weights/phi, offset=offset, family=family)
   if(length(mfit$null.deviance)>1) mfit$null.deviance <- mfit$null.deviance$null.deviance
   if(ykeep) mfit$y <- y
   if(xkeep) mfit$x <- X
   class(mfit) <- c("glm","lm")

#
#  Output for dispersion model:
#  As for glm.object except that prior.weights are not relevant.  Is
#  nested in one output component.
#
   dfit$family <- dfamily
   dfit$prior.weights <- rep(1, N)
   dfit$linear.predictors <- dfit$fitted.values+doffset
   dfit$fitted.values <- phi
   dfit$terms <- dterms
   dfit$aic <- NA
   call$formula <- call$dformula
   call$dformula <- NULL
   call$family <- call(dfamily$family,link=name.dlink)
   dfit$call <- call
   dfit$residuals <- dfamily$dev.resid(d, phi, wt=rep(1/2,N) )
   dfit$deviance <- sum( dfit$residuals  )
if ( verbose ) {
	cat("Hmm2?\n")
	print(doffset) 	
	print(dfamily)
	cat("y=",y,"\n")
	cat("x=",Z,"\n")
	cat("fitted values:",dfit$fitted.values,"\n")
	cat("deviance:",dfit$deviance,"\n")
	cat("residuals:",dfit$residuals,"\n")
	cat("d=",d,"\n")
}
   dfit$null.deviance <- glm.fit(x=Z, y=d, weights=rep(1/2,N), offset=doffset, family=dfamily)

if ( verbose ) cat("Hmm2?\n")
   if(length(dfit$null.deviance)>1) dfit$null.deviance <- dfit$null.deviance$null.deviance
if ( verbose ) cat("Hmm2?\n")
   if(ykeep) dfit$y <- d
   if(zkeep) dfit$z <- Z
   dfit$formula <- as.vector(attr(dterms, "formula"))
   dfit$iter <- iter
   class(dfit) <- c("glm","lm")
   out <- c(mfit, list(dispersion.fit = dfit, iter=iter, method=method, m2loglik=m2loglik))
   class(out) <- c("dglm","glm","lm")
   out
}














######################################################################################
dglm.control <- function(epsilon = 1e-007, maxit = 50, trace = FALSE, ...)
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






dglm.constant <- function(y,family,weights=1)
{
#  Constant term appearing in glm log-likelihood
#  Used by dglm
#  GKS  6 Jan 98, 4 Jul 98, 23 Sep 99.
# "Binomial" chanegd to "binomial": PKD 05 Sep 2006
#
#  Exact cases (in binomial case, exact for phi near 1)
#

   if ( family$family=="Tweedie"){ 
      tweedie.p <- get("tweedie.p", envir = parent.frame() )
   }
   
   const <- switch(family$family[1],
      "Gaussian" = length(y)*log(2*pi),
      "Poisson" = 2*sum(y-y*ifelse(y>0,log(y),0)+lgamma(y+1)),
      "Gamma" = 2*sum(log(y)),
      "Inverse Gaussian" = sum(log(2*pi*y^3)),
      "Tweedie" = switch(match( tweedie.p,c(0,1,2,3),nomatch=0),
         length(y)*log(2*pi),
         2*sum(y-y*ifelse(y>0,log(y),0)+lgamma(y+1)),
         2*sum(log(y)),
         sum(log(2*pi*y^3)) ),
      "binomial" = -2*sum(lgamma(weights+1)-lgamma(weights*y+1)-
      lgamma(weights*(1-y)+1)+weights*(y*ifelse(y>0,log(y),0)+
      (1-y)*ifelse(1-y>0,log(1-y),0)))+sum(log(weights))
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





summary.dglm <- function(object, dispersion=NULL, correlation = FALSE, ...)
{
#  Summarize a double glm
#  GKS  7 Jan 98
   sm <- summary.glm(object,dispersion=dispersion,correlation=correlation, ...)
   sd <- summary.glm(object$dispersion.fit,dispersion=2,correlation=correlation, ...)
   ans <- structure( c(sm, list(dispersion.summary=sd, outer.iter=object$iter,
      m2loglik=object$m2loglik)), class=c("summary.dglm","summary.glm"))
      
   class(ans) <- "summary.dglm"
   return(ans)
}






print.summary.dglm <- function(x, ..., digits = NULL, quote = TRUE, prefix = "", residuals = FALSE)
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
#   cat("\nMEAN MODEL")
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

   cat(paste("(Dispersion Parameters for", x$family$family, 
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
         print(correl[-1,  - p, drop = FALSE], quote = FALSE, digits = digits)
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
   cat(paste("(Dispersion parameter for", xd$family$family, 
      "family taken to be", format(round(xd$dispersion, digits)), ")\n"))
   int <- attr(xd$terms, "intercept")
   if(is.null(int))
      int <- 1
   cat("\n    Scaled Null Deviance:", format(round(xd$null.deviance, digits)), "on", n - 
      int, "degrees of freedom\n")
   cat("Scaled Residual Deviance:", format(round(xd$deviance, digits)), "on", round(
      rdf, digits), "degrees of freedom\n")
#   cat("\nNumber of Fisher Scoring Iterations:", format(trunc(xd$iter)), "\n")
   if(!is.null(correl)) {
      p <- dim(correl)[2]
      if(p > 1) {
         cat("\nCorrelation of Coefficients:\n")
         ll <- lower.tri(correl)
         correl[ll] <- format(round(correl[ll], digits))
         correl[!ll] <- ""
         print(correl[-1,  - p, drop = FALSE], quote = FALSE, digits = digits)
      }
   }
#
#  Overall iteration
#
   cat("\nMinus Twice the Log-Likelihood:", format(round(x$m2loglik, digits)), "\n")
   cat("Number of Alternating Iterations:", format(trunc(x$outer.iter)), "\n")
   invisible(NULL)
}








anova.dglm <- function(object, ..., dispersion = NULL, test = NULL)
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
      paste("\n",object$family$family," double generalized linear model",sep=""),
      paste("\nResponse: ", response,"\n", sep = "") )
   aod <- data.frame(row.names = c("Mean model","Dispersion model"),
      DF = df,
      Seq.Chisq = seqdev,
      Seq.P = 1-pchisq(seqdev,ifelse(df>0,df,NA)),
      Adj.Chisq=adjdev,
      Adj.P = 1-pchisq(adjdev,ifelse(df>0,df,NA)))
   structure(aod, heading=heading, class=c("anova","data.frame") )
   
}
