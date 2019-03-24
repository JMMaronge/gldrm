################################################################################
# These functions are only used for the 2 term approximate information
# (method = "approx2")
################################################################################

## Computes inverse of a 2x2 matrix
invert2by2 <- function(m) {
    matrix(c(m[4], -m[2], -m[3], m[1]), nrow=2) /
        (m[1] * m[4] - m[2] * m[3])
}

## Computes (A+BCB')^{-1}, where Ainv is available and B is rank 1 or 2
## Adiag is an indicator of whether A is a diagonal matrix
woodbury <- function(Ainv, Cinv, B) {
    AinvB <- Ainv %*% B
    mid <- Cinv + crossprod(B, AinvB)
    midinv <- invert2by2(mid)
    inv <- Ainv - AinvB %*% tcrossprod(midinv, AinvB)
    inv
}
################################################################################

#' Control arguments for f0 update algorithm
#'
#' This function returns control arguments for the \eqn{f_0} update algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold. The update has converged when the relative
#' change in log-likelihood between iterations is less than \code{eps}.
#' absolute change is less than \code{thesh}.
#' @param maxiter Maximum number of iterations allowed.
#' @param maxhalf Maximum number of half steps allowed per iteration if
#' log-likelihood does not improve between iterations.
#' @param maxlogstep Maximum optimization step size allowed on the
#' \code{log(f0)} scale.
#'
#' @return Object of S3 class "f0Control", which is a list of control arguments.
#'
#' @export
f0.control <- function(eps=1e-10, maxiter=1000, maxhalf=20, maxlogstep=2, trueHess=FALSE)
{
    f0Control <- as.list(environment())
    class(f0Control) <- "f0Control"
    f0Control
}

#' f0 optimization routine
#'
#' @param y Vector of response values.
#' @param spt Vector of unique observed support points in the response.
#' @param ySptIndex Index of each \code{y} value within \code{spt}.
#' @param sptFreq Vector containing frequency of each \code{spt} value.
#' @param sampprobs Optional matrix of sampling probabilities.
#' @param mu Fitted mean for each observation. Only used if \code{sampprobs=NULL}.
#' @param mu0 Mean constraing for f0.
#' @param f0Start Starting f0 values. (Typically the estimate from the previous
#' iteration.)
#' @param thStart Starting theta values. Needs to be a list of values matching
#' the output of the \code{getTheta} function.
#' @param thetaControl A "thetaControl" object returned from the \code{theta.control}
#' function.
#' @param f0Control An "f0Control" object returned from the \code{f0.control}
#' function.
#' trace Logical. If TRUE, then progress is printed to terminal at each iteration.
#'
#' @return A list containing the following:
#' \itemize{
#' \item \code{f0} Updated values.
#' \item \code{llik} Updated log-likelihood.
#' \item \code{th} Updated list returned from the \code{getTheta} function.
#' \item \code{conv} Convergence indicator.
#' \item \code{iter} Number of iterations until convergence.
#' \item \code{nhalf} The number of half steps taken on the last iteration if the
#' initial BFGS update did not improve the log-likelihood.
#' \item \code{score.log} Score function with respect to log(f0) at convergence.
#' \item \code{info.log} Information matrix with respect to log(f0) at convergence.
#' }
#'
#' @keywords internal
#' @export
getf0 <- function(x, y, spt, ySptIndex, sptFreq, sampprobs, effInfo, beta, offset, dmudeta, mu, mu0, f0Start, thStart,
	thetaControl=theta.control(), f0Control=f0.control(), trace=FALSE)
{
    # Initialize nhalf to prevent error when maxiter=0
    nhalf <- 0
	## Extract theta control arguments
	if (class(f0Control) != "f0Control")
    stop("f0Control must be an object of class f0Control returned by f0Control() function.")
	eps <- f0Control$eps
	maxiter <- f0Control$maxiter
	maxhalf <- f0Control$maxhalf
	maxlogstep <- f0Control$maxlogstep
	trueHess <- f0Control$trueHess

	f0 <- f0Start  # assumes sum(f0Start) = 1 and sum(f0Start * spt) = mu0
	th <- thStart
	llik <- th$llik
	score.log <- NULL
	scoreNorm.log.byIter <- NULL #this is for debugging convergence of score, may remove later
	llik.byIter <- th$llik #this is for debugging convergence of score, may remove later
	score.logT3.iter <- vector(length=length(unique(y)))
    conv <- FALSE
    iter <- 0	
    while (!conv && iter<maxiter) {
		
        iter <- iter + 1

        # Score calculation
        score.logOld <- score.log
		
		if (is.null(sampprobs)) {
	        smm <- outer(spt, mu, "-")
	        ymm <- y - mu
	        yeqmu <- which(abs(ymm) < 1e-15)
			
	        fTiltSWSums <- rowSums(th$fTiltSW)
	        smmfTiltSW <- smm * th$fTiltSW
	        ystd <- ymm / th$bPrime2SW
	        ystd[yeqmu] <- 0  # prevent 0/0
	        score.logT1 <- sptFreq
	        score.logT2 <- fTiltSWSums
	        score.logT3 <- c(smmfTiltSW %*% ystd)
	        score.log <- score.logT1 - score.logT2 - score.logT3	
			
	        if (iter == 1) {
	            d1 <- min(fTiltSWSums)  # max inverse diagonal of first information term, on log scale
	            d2 <- max(abs(score.log)) / maxlogstep
	            d <- max(d1, d2)
	            infoinvBFGS.log <- diag(1/d, nrow=length(f0))
				if (trueHess==TRUE){
			    info.logT1 <- diag(fTiltSWSums)
			    info.logT2 <- tcrossprod(th$fTiltSW)
				info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
			    #info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(ystd, each=nrow(smmfTiltSW))) #I don't think this is right -- ystd includes (y_i-\mu_i) 
			    info.log <- info.logT1 - info.logT2 - info.logT3	
				infoinvBFGS.log <- solve(info.log) # not BFGS, using exact hessian, just using this to make code easier		
				}
	        } else {
	            scorestep.log <- score.log - score.logOld
	            f0step.log <- log(f0) - log(f0old)
	            sy <- sum(f0step.log * scorestep.log)
	            yiy <- c(crossprod(scorestep.log, infoinvBFGS.log %*% scorestep.log))
	            iys <- tcrossprod(infoinvBFGS.log %*% scorestep.log, f0step.log)
				if (trueHess==TRUE){ # use true hessian
				    info.logT1 <- diag(fTiltSWSums)
				    info.logT2 <- tcrossprod(th$fTiltSW)
					info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
				    #info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(ystd, each=nrow(smmfTiltSW))) #I don't think this is right -- ystd includes (y_i-\mu_i) 
				    info.log <- info.logT1 - info.logT2 - info.logT3	
					infoinvBFGS.log <- solve(info.log) # not BFGS, using exact hessian, just using this to make code easier	
				} else {
	            infoinvBFGS.log <- infoinvBFGS.log + ((yiy - sy) / sy^2) * tcrossprod(f0step.log) - (1 / sy) * (iys + t(iys))
				}	
	        }
	        logstep <- c(infoinvBFGS.log %*% score.log)

	        # Cap log(f0) step size
	        logstep.max <- max(abs(logstep))
	        if (logstep.max > maxlogstep)
	            logstep <- logstep * (maxlogstep / logstep.max)

	        # Save values from previous iteration
	        f0old <- f0
	        thold <- th
	        llikold <- llik

	        # Take update step
	        f0 <- exp(log(f0) + logstep)
	        # Scale and tilt f0
	        f0 <- f0 / sum(f0)
	        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1,
	                       thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	        # Update theta and likelihood
	        thold <- th
	        llikold <- llik
	        th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
	                       thetaStart=th$theta, thetaControl=thetaControl)
	        llik <- th$llik
	        conv <- abs((llik - llikold) / (llik + 1e-100)) < eps

	        # If log-likelihood does not improve, change step direction to be along gradient
	        # Take half steps until likelihood improves
	        # Continue taking half steps until log likelihood no longer improves
	        nhalf <- 0
	        if (llik<llikold) {
	            llikprev <- -Inf
	            while ((llik<llikold || llik>llikprev) && nhalf<maxhalf) {
	                nhalf <- nhalf + 1

	                # Set previous values
	                llikprev <- llik
	                thprev <- th
	                f0prev <- f0
	                infoinvBFGS.logprev <- infoinvBFGS.log

	                f0 <- exp((log(f0) + log(f0old)) / 2)
	                f0 <- f0 / sum(f0)
	                f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1,
	                               thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	                th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
	                               thetaStart=th$theta, thetaControl=thetaControl)
	                llik <- th$llik
	                infoinvBFGS.log <- infoinvBFGS.log / 2
	            }

	            if (llik < llikprev) {
	                nhalf <- nhalf - 1
	                llik <- llikprev
	                th <- thprev
	                f0 <- f0prev
	                infoinvBFGS.log <- infoinvBFGS.logprev
	            }

	            conv <- abs((llik - llikold) / (llik + 1e-100)) < eps
	        }

	        if (llik < llikold) {
	            f0 <- f0old
	            th <- thold
	            llik <- llikold
	            conv <- TRUE
	        }

	        if (trace) {
	            printout <- paste0("iter ", iter, ": llik=", llik)
	            if (nhalf > 0)
	                printout <- paste0(printout, "; ", nhalf, " half steps")
	            cat(printout, "\n")
	        }
		} else { # code added for ODS by JMM 02/07/19
            #smm <- outer(spt, th$bPrimeSW, "-") # this isn't right, need (s_m-\mu_i) NOT (s_m - \mu_i^*)
			smm <- outer(spt, th$bPrime, "-") # corrected from MW's original code
            ymm <- y - th$bPrimeSW # this is (y_i - \mu_i^*)
            yeqmu <- which(abs(ymm) < 1e-15)
			
	        fTiltSWSums <- rowSums(th$fTiltSW)
	        smmfTiltSW <- smm * th$fTiltSW
	        fTiltSums <- rowSums(th$fTilt)
	        smmfTilt <- smm * th$fTilt
	        ystd <- ymm / th$bPrime2SW # this calculates (y_i - \mu_i^*)/b^*''(\theta_i), not in the score function under ODS
			#ystd <- ymm / th$bPrime2 # this calculates (y_i - \mu_i^*)/b''(\theta_i) [notice we no longer have b^*''(\theta_i)]
	        ystd[yeqmu] <- 0  # prevent 0/0
	        score.logT1 <- sptFreq
	        score.logT2 <- fTiltSWSums
	        #score.logT3 <- c(smmfTiltSW %*% ystd) # this isn't right, using ODS versions, actual score uses SRS versions (for ODS score)
			score.logT3 <- c(smmfTilt %*% ystd)
	        score.log <- score.logT1 - score.logT2 - score.logT3
			
	        if (iter == 1) {
	            d1 <- min(fTiltSWSums)  # max inverse diagonal of first information term, on log scale
	            d2 <- max(abs(score.log)) / maxlogstep
	            d <- max(d1, d2)
	            infoinvBFGS.log <- diag(1/d, nrow=length(f0))
				if (trueHess==TRUE){
					smm <- outer(spt, th$bPrime, "-")
					smmStar <-  outer(spt, th$bPrimeSW, "-")
					smmfTiltSW <- smm * th$fTiltSW
					smmfTilt <- smm * th$fTilt
					smmStarfTiltSW <- smmStar * th$fTiltSW
		
					info.logT1 <- diag(fTiltSums)
					info.logT2 <- tcrossprod(th$fTiltSW)
			    	info.logT3.1 <- tcrossprod(smmStarfTiltSW, smm * th$fTilt * rep(1/th$bPrime2, each=nrow(smmfTiltSW))) # under ODS, 3rd term decomposes into 3 more terms
					info.logT3.2 <- tcrossprod(smmfTilt, smmStarfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
					info.logT3.3 <- tcrossprod(smmfTilt, smmfTilt* rep(th$bPrime2SW/(th$bPrime2)^2, each=nrow(smmfTiltSW)))
					info.logT3 <- info.logT3.1 + info.logT3.2 - info.logT3.3
			    	info.log <- info.logT1 - info.logT2 - info.logT3	
					infoinvBFGS.log <- solve(info.log) # not BFGS, using true hessian, but this makes coding easier
				}
	        } else {
	            scorestep.log <- score.log - score.logOld
	            f0step.log <- log(f0) - log(f0old)
	            sy <- sum(f0step.log * scorestep.log)
	            yiy <- c(crossprod(scorestep.log, infoinvBFGS.log %*% scorestep.log))
	            iys <- tcrossprod(infoinvBFGS.log %*% scorestep.log, f0step.log)
				if (trueHess==TRUE){
					smm <- outer(spt, th$bPrime, "-")
					smmStar <-  outer(spt, th$bPrimeSW, "-")
					smmfTiltSW <- smm * th$fTiltSW
					smmfTilt <- smm * th$fTilt
					smmStarfTiltSW <- smmStar * th$fTiltSW
		
					info.logT1 <- diag(fTiltSums)
					info.logT2 <- tcrossprod(th$fTiltSW)
			    	info.logT3.1 <- tcrossprod(smmStarfTiltSW, smm * th$fTilt * rep(1/th$bPrime2, each=nrow(smmfTiltSW))) # under ODS, 3rd term decomposes into 3 more terms
					info.logT3.2 <- tcrossprod(smmfTilt, smmStarfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
					info.logT3.3 <- tcrossprod(smmfTilt, smmfTilt* rep(th$bPrime2SW/(th$bPrime2)^2, each=nrow(smmfTiltSW)))
					info.logT3 <- info.logT3.1 + info.logT3.2 - info.logT3.3
			    	info.log <- info.logT1 - info.logT2 - info.logT3	
					infoinvBFGS.log <- solve(info.log) # not BFGS, using true hessian, but this makes coding easier
				} else{
	            infoinvBFGS.log <- infoinvBFGS.log + ((yiy - sy) / sy^2) * tcrossprod(f0step.log) - (1 / sy) * (iys + t(iys))
			}
	        }
	        logstep <- c(infoinvBFGS.log %*% score.log)

	        # Cap log(f0) step size
	        logstep.max <- max(abs(logstep))
	        if (logstep.max > maxlogstep)
	            logstep <- logstep * (maxlogstep / logstep.max)

	        # Save values from previous iteration
	        f0old <- f0
	        thold <- th
	        llikold <- llik

	        # Take update step
	        f0 <- exp(log(f0) + logstep)
	        # Scale and tilt f0
	        f0 <- f0 / sum(f0)
	        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1,
	                       thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	        # Update theta and likelihood
	        thold <- th
	        llikold <- llik
	        th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
	                       thetaStart=th$theta, thetaControl=thetaControl)
	        llik <- th$llik
	        conv <- abs((llik - llikold) / (llik + 1e-100)) < eps

	        # If log-likelihood does not improve, change step direction to be along gradient
	        # Take half steps until likelihood improves
	        # Continue taking half steps until log likelihood no longer improves
	        nhalf <- 0
	        if (llik<llikold) {
	            llikprev <- -Inf
	            while ((llik<llikold || llik>llikprev) && nhalf<maxhalf) {
	                nhalf <- nhalf + 1

	                # Set previous values
	                llikprev <- llik
	                thprev <- th
	                f0prev <- f0
	                infoinvBFGS.logprev <- infoinvBFGS.log

	                f0 <- exp((log(f0) + log(f0old)) / 2)
	                f0 <- f0 / sum(f0)
	                f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1,
	                               thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	                th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
	                               thetaStart=th$theta, thetaControl=thetaControl)
	                llik <- th$llik
	                infoinvBFGS.log <- infoinvBFGS.log / 2
	            }

	            if (llik < llikprev) {
	                nhalf <- nhalf - 1
	                llik <- llikprev
	                th <- thprev
	                f0 <- f0prev
	                infoinvBFGS.log <- infoinvBFGS.logprev
	            }

	            conv <- abs((llik - llikold) / (llik + 1e-100)) < eps
	        }

	        if (llik < llikold) {
	            f0 <- f0old
	            th <- thold
	            llik <- llikold
	            conv <- TRUE
	        }

	        if (trace) {
	            printout <- paste0("iter ", iter, ": llik=", llik)
	            if (nhalf > 0)
	                printout <- paste0(printout, "; ", nhalf, " half steps")
	            cat(printout, "\n")
	        }
			
			
		}
	scoreNorm.log.byIter <- c(scoreNorm.log.byIter, sqrt(sum(score.log^2)))	# added for debugging
	llik.byIter <- c(llik.byIter, llik) # added for debugging
	score.logT3.iter <- cbind(score.logT3.iter, score.logT3) #added for debugging
	}			  
				  
    # Final score calculation
    if (is.null(sampprobs)) {
        smm <- outer(spt, mu, "-")
        ymm <- y - mu
        yeqmu <- which(abs(ymm) < 1e-15)
		
        fTiltSWSums <- rowSums(th$fTiltSW)
        smmfTiltSW <- smm * th$fTiltSW
        ystd <- ymm / th$bPrime2SW
        ystd[yeqmu] <- 0  # prevent 0/0
        score.logT1 <- sptFreq
        score.logT2 <- fTiltSWSums
        score.logT3 <- c(smmfTiltSW %*% ystd)
        score.log <- score.logT1 - score.logT2 - score.logT3	
		scoreNorm.log.byIter <- c(scoreNorm.log.byIter, sqrt(sum(score.log^2))) #added for debugging scoref0
		score.logT3.iter <- cbind(score.logT3.iter, score.logT3) #added for debugging
    } else{#smm <- outer(spt, th$bPrimeSW, "-") # this isn't right, need (s_m-\mu_i) NOT (s_m - \mu_i^*)
			smm <- outer(spt, th$bPrime, "-") # corrected from MW's original code
            ymm <- y - th$bPrimeSW # this is (y_i - \mu_i^*)
            yeqmu <- which(abs(ymm) < 1e-15)
			
	        fTiltSWSums <- rowSums(th$fTiltSW)
	        smmfTiltSW <- smm * th$fTiltSW
	        fTiltSums <- rowSums(th$fTilt)
	        smmfTilt <- smm * th$fTilt
	        #ystd <- ymm / th$bPrime2SW # this calculates (y_i - \mu_i^*)/b^*''(\theta_i), not in the score function under ODS
			ystd <- ymm / th$bPrime2 # this calculates (y_i - \mu_i^*)/b''(\theta_i) [notice we no longer have b^*''(\theta_i)]
	        ystd[yeqmu] <- 0  # prevent 0/0
	        score.logT1 <- sptFreq
	        score.logT2 <- fTiltSWSums
	        #score.logT3 <- c(smmfTiltSW %*% ystd) # this isn't right, using ODS versions, actual score uses SRS versions (for ODS score)
			score.logT3 <- c(smmfTilt %*% ystd)
	        score.log <- score.logT1 - score.logT2 - score.logT3
			scoreNorm.log.byIter <- c(scoreNorm.log.byIter, sqrt(sum(score.log^2))) # added for debugging score.f0		
			score.logT3.iter <- cbind(score.logT3.iter, score.logT3) #added for debugging
    }

    # Final info calculation
	if (is.null(sampprobs)) {
    	info.logT1 <- diag(fTiltSWSums)
    	info.logT2 <- tcrossprod(th$fTiltSW)
		info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
    	#info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(ystd, each=nrow(smmfTiltSW))) #I don't think this is right -- ystd includes (y_i-\mu_i) 
    	info.log <- info.logT1 - info.logT2 - info.logT3
	} else{
		smm <- outer(spt, th$bPrime, "-")
		smmStar <-  outer(spt, th$bPrimeSW, "-")
		smmfTiltSW <- smm * th$fTiltSW
		smmfTilt <- smm * th$fTilt
		smmStarfTiltSW <- smmStar * th$fTiltSW
		
		info.logT1 <- diag(fTiltSums)
		info.logT2 <- tcrossprod(th$fTiltSW)
    	info.logT3.1 <- tcrossprod(smmStarfTiltSW, smm * th$fTilt * rep(1/th$bPrime2, each=nrow(smmfTiltSW))) # under ODS, 3rd term decomposes into 3 more terms
		info.logT3.2 <- tcrossprod(smmfTilt, smmStarfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
		info.logT3.3 <- tcrossprod(smmfTilt, smmfTilt* rep(th$bPrime2SW/(th$bPrime2)^2, each=nrow(smmfTiltSW)))
		info.logT3 <- info.logT3.1 + info.logT3.2 - info.logT3.3
    	info.log <- info.logT1 - info.logT2 - info.logT3
	}
	if (effInfo==TRUE){
	q <- th$bPrime2SW/th$bPrime2
	crossinfo.log <-  (t(x)%*%diag(dmudeta*(1/th$bPrime2)))%*%(t(smmStar*th$fTiltSW - (smm*th$fTilt*q)))
	}else{crossinfo.log <- 0}
	
    list(f0=f0, llik=llik, th=th, conv=conv, iter=iter, nhalf=nhalf,
         score.log=score.log, info.log=info.log, crossinfo.log=crossinfo.log, scoreNorm.log.byIter=scoreNorm.log.byIter, llik.byIter=llik.byIter, score.logT3.iter=score.logT3.iter[,-1])# last three listed items added for debugging of f0 score
}			  
				  
				  
				  
				  
				  
				  
				  
				  
 	
