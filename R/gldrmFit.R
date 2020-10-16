#' Control arguments for \code{gldrm} algorithm
#'
#' This function returns control arguments for the \code{gldrm} algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.

#' @param eps Convergence threshold. The fitting algorithm has converged when the
#' relative change in log-likelihood between iterations is less than \code{eps}.
#' A single iteration consists of a \code{beta} update followed by an \code{f0}
#' update.
#' @param maxiter Maximum number of iterations allowed.
#' @param returnfTiltMatrix Logical. Return nonparametric fitted probabilities for
#' each observation. This is a matrix with nrow equal to the number of
#' observations and ncol equal to the number of unique response values observed.
#' @param returnf0ScoreInfo Logical. If \code{TRUE}, the score and information for
#' \code{log(f0)} are returned as components of the "gldrm" object.
#' @param betaStart Optional vector of starting values for \code{beta}. If the
#' call to gldrm contains a formula, the values of betaStart should correspond to
#' the columns of the model matrix.
#' @param f0Start Optional vector of starting values for \code{f0}. The length
#' of the vector should be the number of unique values in the response, and the
#' vector should correspond to these values sorted in increasing order. The starting
#' values will be scaled to sum to one and tilted to have mean \code{mu0}. All values
#' should be strictly positive.
#' @param print Logical. If \code{TRUE}, the relative change in the log-likelihood
#' will be printed after each iteration.
#'
#' @return Object of S3 class "gldrmControl", which is a list of control arguments.
#'
#' @export
gldrm.control <- function(eps=1e-10, maxiter=100, returnfTiltMatrix=TRUE,
                          returnf0ScoreInfo=FALSE, 	constrainedf0Var=FALSE, print=FALSE,
                          betaStart=NULL, f0Start=NULL)
{
    gldrmControl <- as.list(environment())
    class(gldrmControl) <- "gldrmControl"
    gldrmControl
}

nullspace <- function(x){ # function to calculate nullspace of matrix x
    if (!is.numeric(x)) 
        stop("Argument 'M' must be a numeric matrix.")
    if (is.vector(x)) 
        x <- matrix(c(x), nrow = length(x), ncol = 1)
    qrx <- qr(t(x))
    rnk <- qrx$rank
    if (rnk == ncol(x)) 
        return(NULL)
    ind <- if (rnk == 0) 
        1:ncol(x)
    else -(1:rnk)
    qrQ <- qr.Q(qrx, complete = TRUE)[, ind, drop = FALSE]
    if (length(qrQ) == 0) 
        return(NULL)
    else return(qrQ)
}


#' Main optimization function
#'
#' This function is called by the main \code{gldrm} function.
#'
#' @keywords internal
gldrmFit <- function(x, y, linkfun, linkinv, mu.eta, mu0=NULL, offset=NULL, sampprobs=NULL, effInfo=FALSE, spt= NULL,
					 estprobs = FALSE, fullDataSize=NULL, sampind=NULL, sampgroups=NULL, groups=NULL, xiInfo=NULL,
                     gldrmControl=gldrm.control(), thetaControl=theta.control(),
                     betaControl=beta.control(), f0Control=f0.control())
{
	
	
	
    ## Extract control arguments
    if (class(gldrmControl) != "gldrmControl")
        stop("gldrmControl must be an object of class \'gldrmControl\' returned by
              gldrmControl() function.")
    eps <- gldrmControl$eps
    maxiter <- gldrmControl$maxiter
    returnHess <- gldrmControl$returnHess
    returnfTiltMatrix <- gldrmControl$returnfTiltMatrix
    returnf0ScoreInfo <- gldrmControl$returnf0ScoreInfo
    print <- gldrmControl$print
    betaStart <- gldrmControl$betaStart
    f0Start <- gldrmControl$f0Start
	  constrainedf0Var <- gldrmControl$constrainedf0Var
		
	sampprobs.for.f0start <- sampprobs #it appears that sampprobs gets rewritten below - I am adding this line for the f0start calculation under ODS

    ## Tabulation and summary of responses used in estimating f0
    n <- length(y)
    if(is.null(spt)){
    spt <- sort(unique(y))}  # observed support
    ySptIndex <- match(y, spt)# index of each y value within support
    sptFreq <- table(factor(ySptIndex, levels = 1:length(spt)))
    attributes(sptFreq) <- NULL

    ## Check sampprobs
    if (!is.null(sampprobs)) {
        if (!(is.vector(sampprobs) || is.matrix(sampprobs)) || !is.numeric(sampprobs) || any(sampprobs < 0))
            stop("sampprobs must be a matrix or vector of nonnegative numeric values")

        if (is.vector(sampprobs)) {
            if (length(sampprobs) != length(spt)){
                stop(paste0("sampprobs vector should have length equal to the ",
                            "number of unique observed values in the response OR have length equal to specified spt vector."))
            }
            sampprobs <- matrix(sampprobs, nrow=n, ncol=length(sampprobs), byrow=TRUE) #sampprobs is now a matrix where each row denotes subject and each column denotes a support point, so each entry is the sampling probability for a particular subject at a particular support point
          } else {
             #sampprobs must be a matrix
            if (nrow(sampprobs) != n)
                stop(paste0("sampprobs matrix should have row dimension equal to ",
                            "the number of observations."))
            if (ncol(sampprobs) != length(spt))
                stop(paste0("sampprobs matrix should have column dimension equal ",
                            "to the number of unique observe values in the response OR have length equal to specified spt vector."))
        }
    }

    ## Initialize offset
    if (is.null(offset))
        offset <- rep(0, n)

    ## Initialize mu0 if not provided by user
    if (is.null(mu0)) {
        mu0 <- mean(y)
        # mu0 <- linkinv(0)
        # mu0 <- mean(range(y))
        # mu0 <- mean(spt)
    } else if (mu0<=min(spt) || mu0>=max(spt)) {
        stop(paste0("mu0 must lie within the range of observed values. Choose a different ",
                    "value or set mu0=NULL to use the default value, mean(y)."))
    }

    ## Initialize f0
    if (is.null(f0Start)) {
		if (!is.null(sampprobs)){
			f0 <- sptFreq / n
			#f0 <- ( (sptFreq / n) / sampprobs.for.f0start ) # divide by sampling probs when doing ODS
			#f0 <- f0 / sum(f0) # renormalize
			if (is.null(mu0)) {
			mu0 <- sum(sort(unique(y))*f0)} # recalculate mean

		}
		else {
        f0 <- sptFreq / n
    
		}
        if (mu0 != mean(y))
            f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1, thetaStart=0,
                           thetaControl=thetaControl)$fTilt[, 1]
						  
	} else {
        
        if (length(f0Start) != length(spt))
            stop("Length of f0Start should equal number of unique values in the response.")
        if (any(f0Start <= 0))
            stop("All values in f0Start should be strictly positive.")
        f0 <- f0Start / sum(f0Start)
        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1, thetaStart=0,
                       thetaControl=thetaControl)$fTilt[, 1]
    }

    ## Initialize beta
    ## The starting values returned by lm.fit guarantee that all mu values are
    ## within the support range, even if there is no intercept.
    ## Offset could still create problems.
    lmcoef <- stats::lm.fit(x, linkfun(mu0) - offset)$coef
    if (is.null(betaStart)) {
        beta <- lmcoef
    } else {
        if (length(betaStart) != ncol(x))
            stop("Length of betaStart should equal the number of columns in the model matrix.")
      beta <- betaStart
    }

    ## Drop coefficients if x is not full rank (add NA values back at the end)
    naID <- is.na(lmcoef)
    beta <- beta[!naID]
    x <- x[, !naID, drop=FALSE]
    eta <- c(x %*% beta + offset)
    mu <- linkinv(eta)
    if (ncol(x) >= n)
    stop("gldrm requires n > p.")
    if (any(mu<min(spt) | mu>max(spt)))
    stop("Unable to find beta starting values that do not violate convex hull condition.")
    ## Get initial theta and log likelihood
    th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
                   thetaStart=NULL, thetaControl=thetaControl)
    llik <- th$llik
	iter.scoresNorm <- list(NA)
	llik.iter <- list(NA)
	scoref0.log.T3 <- list(NA)
	
    conv <- FALSE
    iter <- 0
    while (!conv && iter <= maxiter)
    {
        iter <- iter+1
        betaold <- beta
        f0old <- f0
        llikold <- llik

        ## update beta (mu) and theta, with fixed f0:
		
		
        bb <- getBeta(x=x, y=y, spt=spt, ySptIndex=ySptIndex, f0=f0,
                      linkinv=linkinv, mu.eta=mu.eta, offset=offset, sampprobs=sampprobs,
                      betaStart=beta, thStart=th,
                      thetaControl=thetaControl, betaControl=betaControl)
		
        th <- bb$th
        llik <- bb$llik
        mu <- bb$mu
        beta <- bb$beta
		
		dmudeta <- bb$dmudeta
	

        ## update f0 and theta, with fixed beta (mu)
        ff <- getf0(x=x, y=y, spt=spt, ySptIndex=ySptIndex, sptFreq=sptFreq,
                    sampprobs=sampprobs, estprobs=estprobs, groups=groups,
				    effInfo=effInfo, beta=beta, offset=offset, dmudeta=dmudeta, mu=mu, mu0=mu0, f0Start=f0, thStart=th,
                    thetaControl=thetaControl, f0Control=f0Control)
        th <- ff$th
        llik <- ff$llik
        f0 <- ff$f0
		

        ## Check convergence
        del1 <- abs((llik - llikold) / llik)
        if (llik == 0) del1 <- 0
		#conv <- del1 < eps & del2 < eps
		conv <- del1<eps
        if (print) {
            cat("iteration ", iter,
                "\nrelative change in log-likelihood = ", del,
                "  (eps = ", eps, ")\n")
        }
		
    }

	
    ## Final values
    eta <- linkfun(mu)
    dmudeta <- mu.eta(eta)
    llik <- ff$llik
    theta <- th$theta
    bPrime <- th$bPrime
    bPrime2 <- th$bPrime2
    fTilt <- th$fTilt[cbind(ySptIndex, seq_along(ySptIndex))]
	
	
    ## Compute betaHat variance
    if (!is.null(sampprobs)) {
		q <- (th$bPrime2SW / th$bPrime2)
        w <- dmudeta^2 / th$bPrime2 * q
        wSqrt <- sqrt(w)
		#beta.info.check <- matrix(data=0, nrow=ncol(x), ncol=ncol(x))
		#for(i in 1:nrow(x)){
		#	beta.info.check <- beta.info.check+x[i,]%*%t(x[i,])*(dmudeta[i]^2*(th$bPrime2SW[i]/th$bPrime2[i]^2))}
		
		#print("beta info subj is")
		#print(beta.info.check)	
			
			
		
    } else {
        w <- dmudeta^2 / th$bPrime2
        wSqrt <- sqrt(w)
    }
    if (any(wSqrt == Inf)) {
        ## if any weights are infinite, return all standard errors as zero
        varbeta <- matrix(0, nrow=length(beta), ncol=length(beta))
    } else {
        wtdX <- wSqrt * x
        # varbeta <- chol2inv(qr.R(qr(wtdX)))  # not stable
        # varbeta <- solve(crossprod(wtdX))  # not stable
        varbeta <- tcrossprod(backsolve(qr.R(qr(wtdX)), diag(ncol(wtdX))))
		if (effInfo==TRUE){
			if (estprobs==TRUE){
				
				
				infobeta <- crossprod(wtdX)
				
				#print("beta info matrix is")
				#print(infobeta)
				#print("uncorrected SE")
				#print(sqrt(diag(solve(infobeta))))
				infof0 <- ff$info.log
				infocross <- ff$crossinfo.log	
				#print("f0 info matrix is")
				#print(infof0)
				#print("beta f0 cross info matrix is")
				#print(infocross)
		
				infotheta.1 <- cbind(infobeta,infocross)
				infotheta.2 <- cbind(t(infocross),infof0) 	
				infotheta <- rbind(infotheta.1, infotheta.2)	
				#print(solve(infotheta))
		
				#length.betas <- length(beta)
				#length.f0 <- length(f0)
				#supp.vals <- sort(unique(y))
				#g0 <- log(f0)
				#grad.constraint <- matrix(nrow=2, ncol=(length.betas+length.f0))
				#grad.constraint[1,(1:length.betas)] <- 0
				#grad.constraint[2,(1:length.betas)] <- 0
				#grad.constraint[1,((length.betas+1):((length.betas+length.f0)))] <- exp(g0)	
				#grad.constraint[2,((length.betas+1):((length.betas+length.f0)))] <- supp.vals*exp(g0)
				#U <- nullspace(grad.constraint)
				#U1 <- U[1:length.betas,]
				#U2 <- U[(length.betas+1):((length.betas+length.f0)),]
				#tmp <- t(U1)%*%infobeta%*%U1 + t(U2)%*%t(infocross)%*%U1 + t(U1)%*%infocross%*%U2 + t(U2)%*%infof0%*%U2
				#constrained.info.inv.theta <- U%*%solve(tmp,t(U))
				
				thetaxi.crossinfo <- rbind(ff$betaxi.crossinfo, ff$f0xi.crossinfo.log)
				#print("theta xi cross info")
				#print(thetaxi.crossinfo)
				
				#constrained.info.inv<-  constrained.info.inv.theta- constrained.info.inv.theta%*%thetaxi.crossinfo%*%solve(xiInfo)%*%t(thetaxi.crossinfo)%*%constrained.info.inv.theta
				
				#print("without constraints")
				#print(constrained.info.inv)
				
				#print("without estimating xi")
				#print(constrained.info.inv.theta)
				#print("correction term")
				#print(constrained.info.inv.theta%*%thetaxi.crossinfo%*%solve(xiInfo)%*%t(thetaxi.crossinfo)%*%constrained.info.inv.theta)
				
				D_matrix1 <- cbind(infotheta,thetaxi.crossinfo)
				D_matrix2 <- cbind(matrix(0,nrow=ncol(thetaxi.crossinfo),ncol=nrow(thetaxi.crossinfo)), xiInfo)
				D_matrix <- rbind(D_matrix1, D_matrix2)
				
				V_matrix1 <- D_matrix1
				V_matrix2 <- cbind(t(thetaxi.crossinfo),xiInfo)
				V_matrix <- rbind(V_matrix1,V_matrix2)
				
				
				length.betas <- length(beta)
				length.f0 <- length(f0)
				length.xi <- nrow(xiInfo)
				supp.vals <- sort(unique(y))
				g0 <- log(f0)
				grad.constraint <- matrix(nrow=2, ncol=(length.betas+length.f0+length.xi))
				grad.constraint[1,(1:length.betas)] <- 0
				grad.constraint[2,(1:length.betas)] <- 0
				grad.constraint[1,((length.betas+1):((length.betas+length.f0)))] <- exp(g0)	
				grad.constraint[2,((length.betas+1):((length.betas+length.f0)))] <- supp.vals*exp(g0)
				grad.constraint[1,((length.betas+length.f0+1):((length.betas+length.f0+length.xi)))] <- 0
				grad.constraint[2,((length.betas+length.f0+1):((length.betas+length.f0+length.xi)))] <- 0
				U <- nullspace(grad.constraint)
				#U<- round(U)
				#print("nullspace matrix")
				#print(U)
				#	print(dim(U))
				#U1 <- U[1:length.betas,]
				#U2 <- U[(length.betas+1):((length.betas+length.f0)),]
				#tmp <- t(U1)%*%infobeta%*%U1 + t(U2)%*%t(infocross)%*%U1 + t(U1)%*%infocross%*%U2 + t(U2)%*%infof0%*%U2
				
				#svd.Vmat <- svd(V_matrix)	
				#print(solve(V_matrix))
				#pseudo.inv.Vmat <- svd.Vmat$v%*%(diag(1/svd.Vmat$d)*t(svd.Vmat$u))
				#print(pseudo.inv.Vmat)
				#print("xi info")
				#print(xiInfo)
			
				
				
				svd.v_matrix <- svd(V_matrix)
				psuedo.inv_v <- svd.v_matrix$v%*%((1/svd.v_matrix$d) * t(svd.v_matrix$u))
				unconstrained.info <- t(D_matrix)%*%psuedo.inv_v%*%D_matrix
				#print("inverted v matrix")
				#print(ginv(V_matrix))
				
				
				#print("matrix to invert")
				#print(t(U)%*%unconstrained.info%*%U)
				#print(eigen(t(U)%*%unconstrained.info%*%U))
				constrained.info.inv <- U%*%solve(t(U)%*%unconstrained.info%*%U,t(U))
				
				#constrained.D.matrix.inv <- U%*%solve(t(U)%*%D_matrix%*%U,t(U))
				#constrained.V.matrix.inv <- U%*%solve(t(U)%*%V_matrix%*%U,t(U))
				#svd.constrained.V.inv <- svd(constrained.V.matrix.inv)
	
				
				
				#constrained.V.matrix <- svd.constrained.V.inv$v%*%diag(1/svd.constrained.V.inv$d)%*%t(svd.constrained.V.inv$u)
				#constrained.V.matrix <- U%*%(t(U)%*%V_matrix%*%U)%*%t(U)
				#constrained.info.inv <-constrained.D.matrix.inv%*%constrained.V.matrix%*%t(constrained.D.matrix.inv)
				
				#print("with constraints")
				#print(constrained.info.inv)
				
				
				
			}	
				
				else {
			
			
		# if we want effective information which include constraints!!!	
		infobeta <- crossprod(wtdX)
		#print("beta info matrix is")
		#print(infobeta)
		#print("uncorrected SE")
		#print(sqrt(diag(solve(infobeta))))
		infof0 <- ff$info.log
		infocross <- ff$crossinfo.log	
		
		
		infotheta.1 <- cbind(infobeta,infocross)
		infotheta.2 <- cbind(t(infocross),infof0) 	
		infotheta <- rbind(infotheta.1, infotheta.2)	
		#print(solve(infotheta))
		
		length.betas <- length(beta)
		length.f0 <- length(f0)
		supp.vals <- sort(unique(y))
		g0 <- log(f0)
		grad.constraint <- matrix(nrow=2, ncol=(length.betas+length.f0))
		grad.constraint[1,(1:length.betas)] <- 0
		grad.constraint[2,(1:length.betas)] <- 0
		grad.constraint[1,((length.betas+1):((length.betas+length.f0)))] <- exp(g0)	
		grad.constraint[2,((length.betas+1):((length.betas+length.f0)))] <- supp.vals*exp(g0)
		U <- nullspace(grad.constraint)
		U1 <- U[1:length.betas,]
		U2 <- U[(length.betas+1):((length.betas+length.f0)),]

		
		if(length.betas==1){
		  infobeta <- as.double(infobeta)
		  U1 <- t(U1) # this technically not mathematically correct, doing this because R forces U1 to column vector when it should be row
		  tmp <-  (t(U1)%*%U1)*infobeta + t(U2)%*%t(infocross)%*%U1 + t(U1)%*%infocross%*%U2 + t(U2)%*%infof0%*%U2  
		} else{

		tmp <- t(U1)%*%infobeta%*%U1 + t(U2)%*%t(infocross)%*%U1 + t(U1)%*%infocross%*%U2 + t(U2)%*%infof0%*%U2}
		
		constrained.info.inv <- U%*%solve(tmp,t(U))
		
		
		#print("constrained info is")
		#print(constrained.info.inv)
		if(constrainedf0Var==TRUE){
			constrainedf0Var <- constrained.info.inv[((length.betas+1):((length.betas+length.f0))), ((length.betas+1):((length.betas+length.f0)))]
			b.mat <- matrix(NA, nrow=length.f0-1, ncol=length.f0)
			for(i in 1:(length.f0-1)){
				b.mat[i,] <- c(f0[1:i],rep(0,length.f0-i))
			}
			
			b.mat.test <- matrix(NA, nrow=length.f0, ncol=length.f0)
			for(i in 1:length.f0){
				b.mat.test[i,] <- c(f0[1:i],rep(0,length.f0-i))
			}
			
			constrainedf0_cdfVar <- b.mat%*%constrainedf0Var%*%t(b.mat)
			
			b.mat2 <-  matrix(NA, nrow=length.f0-1, ncol=length.f0)
			for(i in 1:(length.f0-1)){
				
				b.mat2[i,] <- c((f0[1:i]/sum(f0[1:i]))+ (f0[1:i]/(1-sum(f0[1:i]))),rep(0,length.f0-i))
			}
			constrainedf0_cdfVar.logit <- b.mat2%*%constrainedf0Var%*%t(b.mat2)	

		}
		else{
			constrainedf0Var<-NULL
			constrainedf0_cdfVar<- NULL}
		}			
		
		#infobeta.eff <- infobeta - infocross%*%solve(infof0,t(infocross))
		#print("information for beta is")
		#print(infobeta)
		#print("correction is")
		#print(infocross%*%solve(infof0,t(infocross)))
		#print("final info is")
		#print(infobeta - infocross%*%solve(infof0,t(infocross)))
		
		varbeta <- constrained.info.inv[1:length.betas,1:length.betas]
		#varbeta <- solve(infobeta.eff)
	}
    }

    ## Compute standard errors
    seBeta <- sqrt(diag(varbeta))
    print(varbeta)
    seEta <- sqrt(pmax(0, apply(x, 1, function(xx) crossprod(xx, varbeta) %*% xx)))
    seMu <- dmudeta * seEta

    ## Add NA values back into beta vector and varbeta if covariate matrix is not full rank
    nBeta <- length(beta) + sum(naID)
    betaTemp <- seBetaTemp <- rep(NA, nBeta)
    betaTemp[!naID] <- beta
    print(seBeta)
    seBetaTemp[!naID] <- seBeta
    beta <- betaTemp
    seBeta <- seBetaTemp
    varbetaTemp <- matrix(NA, nrow=nBeta, ncol=nBeta)
    varbetaTemp[!naID, !naID] <- varbeta
    varbeta <- varbetaTemp

    # Inference vs. null model
    containsIntercept <- any(apply(x, 2, function(xx) all(xx == xx[1])))
    if (!containsIntercept) {
        llikNull <- lr.stat <- lr.df <- lr.pval <- NA
    } else {
        xrank <- ncol(x)  # columns of x have already been dropped, if necessary to make x full rank
        lr.df <- c(max(xrank-1, 1), n-xrank)  # force df[1] >= 1 just in case the full model is a null model
        llikNull <- sum(sptFreq * log(sptFreq / n))
        lr.stat <- 2 * (llik - llikNull) / lr.df[1]
        lr.pval <- 1 - stats::pf(lr.stat, lr.df[1], lr.df[2])
    }

    ## Return gldrm object
    attributes(beta) <- NULL
    attributes(f0) <- NULL

    fit <- list(conv=conv, iter=iter, llik=llik,
                beta=beta, mu=mu, eta=eta, f0=f0, spt=spt, mu0=mu0,
                varbeta=varbeta, seBeta=seBeta, seMu=seMu, seEta=seEta, constrained.info.inv=constrainedf0Var, #constrainedf0_cdfVar=constrainedf0_cdfVar,
				#constrainedf0_cdfVar.logit=constrainedf0_cdfVar.logit, 
				xiInfo=xiInfo,
                theta=theta, bPrime=bPrime, bPrime2=bPrime2, fTilt=fTilt, sampprobs=sampprobs, dmudeta=dmudeta, mu=mu,
                llikNull=llikNull, lr.stat=lr.stat, lr.df=lr.df, lr.pval=lr.pval)

    if (returnfTiltMatrix)
        fit$fTiltMatrix <- t(th$fTilt)

    if (returnf0ScoreInfo) {
        fit$score.logf0 <- ff$score.log
        fit$info.logf0 <- ff$info.log
		fit$info.logcross <- ff$crossinfo.log
    }

    class(fit) <- "gldrm"
    fit
}
