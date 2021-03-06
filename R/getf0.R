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
f0.control <- function(eps=1e-10, maxiter=1000, maxhalf=20, maxlogstep=2, trace=FALSE)
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
getf0 <- function(x, y, spt, ySptIndex, sptFreq, weights=weights, sampprobs, estprobs, groups, effInfo, beta, offset, dmudeta, mu, mu0, f0Start, thStart,
	thetaControl=theta.control(), f0Control=f0.control()	)
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
	trace <- f0Control$trace

	f0 <- f0Start  # assumes sum(f0Start) = 1 and sum(f0Start * spt) = mu0
	th <- thStart
	llik <- th$llik
	score.log <- NULL
	
	


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
	        score.logT2 <- rowSums(t(weights*t(th$fTiltSW)))
	        score.logT3 <- c(smmfTiltSW %*% (weights*ystd))
	        score.log <- score.logT1 - score.logT2 - score.logT3
	        if (iter == 1) {
	            d1 <- min(fTiltSWSums)  # max inverse diagonal of first information term, on log scale
				#d1 <- nrow(x)*diag(nrow(th$fTilt))
	            d2 <- max(abs(score.log)) / maxlogstep
	            d <- max(d1, d2)
	            infoinvBFGS.log <- diag(1/d, nrow=length(f0))

	        } else {
	            scorestep.log <- score.log - score.logOld
	            f0step.log <- log(f0) - log(f0old)
	            sy <- sum(f0step.log * scorestep.log)
	            yiy <- c(crossprod(scorestep.log, infoinvBFGS.log %*% scorestep.log))
	            iys <- tcrossprod(infoinvBFGS.log %*% scorestep.log, f0step.log)
 
	            infoinvBFGS.log <- infoinvBFGS.log + ((yiy - sy) / sy^2) * tcrossprod(f0step.log) - (1 / sy) * (iys + t(iys))
					
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
	        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, weights=1, sampprobs=NULL, ySptIndex=1,
	                       thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	        # Update theta and likelihood
	        thold <- th
	        llikold <- llik
	        th <- getTheta(spt=spt, f0=f0, mu=mu, weights=weights, sampprobs=sampprobs, ySptIndex=ySptIndex,
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
	                f0 <- getTheta(spt=spt, f0=f0, mu=mu0, weights=1, sampprobs=NULL, ySptIndex=1,
	                               thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	                th <- getTheta(spt=spt, f0=f0, mu=mu, weights=weights, sampprobs=sampprobs, ySptIndex=ySptIndex,
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
	        #ystd <- ymm / th$bPrime2SW # this calculates (y_i - \mu_i^*)/b^*''(\theta_i), not in the score function under ODS
			ystd <- ymm / th$bPrime2 # this calculates (y_i - \mu_i^*)/b''(\theta_i) [notice we no longer have b^*''(\theta_i)]
	        ystd[yeqmu] <- 0  # prevent 0/0
			#adding f0star must sum to 1 to code, not sure if this will work - JMM 06/08/19
	        
			
			score.logT1 <- sptFreq
	        score.logT2 <- rowSums(t(weights*t(th$fTiltSW)))
	        #score.logT3 <- c(smmfTiltSW %*% ystd) # this isn't right, using ODS versions, actual score uses SRS versions (for ODS score)
			score.logT3 <- c(smmfTilt %*% (weights*ystd)) # use this one for matrix calcs when return
	        score.log <- score.logT1 - score.logT2 - score.logT3
			#print("matrix calc")
			#print(score.log)
			
			
			#subject calc
			
			#score.tmp <- matrix(data=NA, nrow=nrow(th$fTiltSW), ncol=nrow(x)) # used to debug f0 score, coding score 1 subject at a time
			#for(i in 1:nrow(x)){
			#	for(k in 1:nrow(th$fTiltSW)){
			#		score1 <- ifelse(y[i]==spt[k],1,0) - ifelse(y[i]==spt[k],1,0)*((exp(log(f0[k]))*sampprobs[i,k])/sum(exp(log(f0))*sampprobs[i,]))
			#		score2 <- th$fTiltSW[k,i] - th$fTiltSW[k,i]*((exp(log(f0[k]))*sampprobs[i,k])/sum(exp(log(f0))*sampprobs[i,]))
			#		score3 <- (y[i]- th$bPrimeSW[i])*th$fTilt[k,i]*(spt[k]-th$bPrime[i])*(1/th$bPrime2[i])
			#		score.tmp[k,i] <- score1 - score2 - score3
			#	}
			
			
			#}
			#score.log <- rowSums(score.tmp)
			#print("subject calc")
			#print(score.log)

			
			
			
	        if (iter == 1) {
	            d1 <- min(fTiltSWSums)  # max inverse diagonal of first information term, on log scale
	            d2 <- max(abs(score.log)) / maxlogstep
	            d <- max(d1, d2)
	            infoinvBFGS.log <- diag(1/d, nrow=length(f0))
	        } else {
	            scorestep.log <- score.log - score.logOld
	            f0step.log <- log(f0) - log(f0old)
	            sy <- sum(f0step.log * scorestep.log)
	            yiy <- c(crossprod(scorestep.log, infoinvBFGS.log %*% scorestep.log))
	            iys <- tcrossprod(infoinvBFGS.log %*% scorestep.log, f0step.log)
	            infoinvBFGS.log <- infoinvBFGS.log + ((yiy - sy) / sy^2) * tcrossprod(f0step.log) - (1 / sy) * (iys + t(iys))
			
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
	        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, weights=1, sampprobs=NULL, ySptIndex=1,
	                       thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	        # Update theta and likelihood
	        thold <- th
	        llikold <- llik
	        th <- getTheta(spt=spt, f0=f0, mu=mu, weights=weights, sampprobs=sampprobs,  ySptIndex=ySptIndex,
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
	                f0 <- getTheta(spt=spt, f0=f0, mu=mu0, weights=1, sampprobs=NULL, ySptIndex=1,
	                               thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
	                th <- getTheta(spt=spt, f0=f0, mu=mu, weights=weights, sampprobs=sampprobs, ySptIndex=ySptIndex,
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
        score.logT2 <- rowSums(t(weights*t(th$fTiltSW)))
        score.logT3 <- c(smmfTiltSW %*% (weights*ystd))
        score.log <- score.logT1 - score.logT2 - score.logT3
			
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
	        score.logT2 <- rowSums(t(weights*t(th$fTiltSW)))
	        #score.logT3 <- c(smmfTiltSW %*% ystd) # this isn't right, using ODS versions, actual score uses SRS versions (for ODS score)
			score.logT3 <- c(smmfTilt %*% (weights*ystd)) # use this one when doing matrix calcs
	        score.log <- score.logT1 - score.logT2 - score.logT3
			#print("matrix score")
			#print(score.log)


		
			#score.tmp <- matrix(data=NA, nrow=nrow(th$fTiltSW), ncol=nrow(x)) # used to debug f0 score, coding score 1 subject at a time
			#for(i in 1:nrow(x)){
			#	for(k in 1:nrow(th$fTiltSW)){
			#		score1 <- ifelse(y[i]==spt[k],1,0) - ifelse(y[i]==spt[k],1,0)*((exp(log(f0[k]))*sampprobs[i,k])/sum(exp(log(f0))*sampprobs[i,]))
			#		score2 <- th$fTiltSW[k,i] - th$fTiltSW[k,i]*((exp(log(f0[k]))*sampprobs[i,k])/sum(exp(log(f0))*sampprobs[i,]))
			#		score3 <- (y[i]- th$bPrimeSW[i])*th$fTilt[k,i]*(spt[k]-th$bPrime[i])*(1/th$bPrime2[i])
			#		score.tmp[k,i] <- score1 - score2 - score3
			#	}
			
			
			#}
			#score.log <- rowSums(score.tmp)
				
					
					
					
			#scoreNorm.log.byIter <- c(scoreNorm.log.byIter, sqrt(sum(score.log^2))) # added for debugging score.f0		
			#score.logT3.iter <- cbind(score.logT3.iter, score.logT3) #added for debugging
     }

    # Final info calculation
	if (is.null(sampprobs)) {
    	#info.logT1 <- diag(fTiltSWSums) # from Mike's code, after debugging the ODS info, this shouldn't be right either
		info.logT1 <- diag(rowSums(t(weights*t(th$fTilt))))
		#print("matrix term")
		#print(info.logT1)
    	info.logT2 <- tcrossprod(sqrt(weights)*th$fTiltSW)
		info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(weights*(1/th$bPrime2), each=nrow(smmfTiltSW)))
    	#info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(ystd, each=nrow(smmfTiltSW))) #I don't think this is right -- ystd includes (y_i-\mu_i) 
    	info.log <- info.logT1 - info.logT2 - info.logT3
		#print("matrix f0 info")
		#print(info.log)
		
		
		#info.log2.t1 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt)) #this chunk of code was used to debug the information
		#info.log2.t2 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#info.log2.t3 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#subjs.info.t1 <- list(length=nrow(x))
		#subjs.info.t2 <- list(length=nrow(x))
		#subjs.info.t3 <- list(length=nrow(x))
		#for(i in 1:nrow(x)){
		#	subj.info.t1 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#	subj.info.t2 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#	subj.info.t3 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#	for(k in 1:nrow(th$fTilt)){
		#		for(m in 1:nrow(th$fTilt)){
		#			t1 <- ifelse(k==m,th$fTilt[k,i],0)
		#			t2 <- th$fTiltSW[k,i]*th$fTiltSW[m,i]
		#			t3.1 <- (spt[k]-th$bPrime[i])*(spt[m]-th$bPrime[i])*th$fTilt[k,i]*th$fTilt[m,i]
		#			t3 <- (1/th$bPrime2[i])*(t3.1)
		#			subj.info.t1[k,m] <- t1
		#			subj.info.t2[k,m] <- t2
		#			subj.info.t3[k,m] <- t3
					
					
		#		}
				
		#	}
		#	subjs.info.t1[[i]] <- subj.info.t1
		#	subjs.info.t2[[i]] <- subj.info.t2
		#	subjs.info.t3[[i]] <- subj.info.t3
			
			
		#}
		#k<-2
		#m<-5
		#print("f tilt k")
		#print(th$fTilt[k,i])
		#print("f tilt m")
		#print(th$fTilt[m,i])
		#print("f tilt sw k")
		#print(th$fTiltSW[k,i])
		#print("f tilt sw m")
		#print(th$fTiltSW[m,i])
		#info.log2.t1 <- Reduce('+', subjs.info.t1)
		#info.log2.t2 <- Reduce('+', subjs.info.t2)
		#info.log2.t3 <- Reduce('+', subjs.info.t3)
		#print("subject term")
		#print(info.log2.t1)
		#info.log2 <- info.log2.t1-info.log2.t2-info.log2.t3
		#print("subject f0 info")
		#print(info.log2)
		
	} else{
		smm <- outer(spt, th$bPrime, "-")
		smmStar <-  outer(spt, th$bPrimeSW, "-")
		smmfTiltSW <- smm * th$fTiltSW
		smmfTilt <- smm * th$fTilt
		smmStarfTiltSW <- smmStar * th$fTiltSW
	
		
		
		#info.logT1 <- diag(fTiltSums) # from Mike's Code, after a lot of time, this IS NOT RIGHT
		info.logT1 <- diag(rowSums(th$fTiltSW))
		info.logT2 <- tcrossprod(th$fTiltSW)
    	info.logT3.1 <- tcrossprod(smmStarfTiltSW, smm * th$fTilt * rep(1/th$bPrime2, each=nrow(smmfTiltSW))) # under ODS, 3rd term decomposes into 3 more terms
		info.logT3.2 <- tcrossprod(smmfTilt, smmStarfTiltSW * rep(1/th$bPrime2, each=nrow(smmfTiltSW)))
		info.logT3.3 <- tcrossprod(smmfTilt, smmfTilt* rep(th$bPrime2SW/(th$bPrime2)^2, each=nrow(smmfTiltSW)))
		info.logT3 <- info.logT3.1 + info.logT3.2 - info.logT3.3
		#info.logT3 <- tcrossprod(smmStarfTiltSW, smmStarfTiltSW*rep(1/th$bPrime2SW, each=nrow(smmStarfTiltSW)))
    	info.log <- info.logT1 - info.logT2 - info.logT3
		#print("matrix f0 info")
		#print(info.log)
		
		#info.log2.t1 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt)) #this chunk of code was used to debug the information
		#info.log2.t2 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#info.log2.t3 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#subjs.info.t1 <- list(length=nrow(x))
		#subjs.info.t2 <- list(length=nrow(x))
		#subjs.info.t3 <- list(length=nrow(x))
		#for(i in 1:nrow(x)){
		#	subj.info.t1 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#	subj.info.t2 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#	subj.info.t3 <- matrix(data=NA, nrow=nrow(th$fTilt), ncol=nrow(th$fTilt))
		#	for(k in 1:nrow(th$fTilt)){
		#		for(m in 1:nrow(th$fTilt)){
		#			t1 <- ifelse(k==m,1,0)
		#			t2 <- th$fTiltSW[k,i]*th$fTiltSW[m,i]
		#			t3.1 <- (spt[k]-th$bPrimeSW[i])*(spt[m]-th$bPrime[i])*th$fTiltSW[k,i]*th$fTilt[m,i]
		#			t3.2 <- (spt[m]-th$bPrimeSW[i])*(spt[k]-th$bPrime[i])*th$fTilt[k,i]*th$fTiltSW[m,i]
		#			t3.3 <- (spt[k]-th$bPrime[i])*(spt[m]-th$bPrime[i])*th$fTilt[k,i]*th$fTilt[m,i]*(th$bPrime2SW[i]/th$bPrime2[i])
		#			t3 <- (1/th$bPrime2[i])*(t3.1+t3.2-t3.3)
		#			subj.info.t1[k,m] <- t1
		#			subj.info.t2[k,m] <- t2
		#			subj.info.t3[k,m] <- t3
					
					
		#		}
				
		#	}
		#	subjs.info.t1[[i]] <- subj.info.t1
		#	subjs.info.t2[[i]] <- subj.info.t2
		#	subjs.info.t3[[i]] <- subj.info.t3
			
			
		#}
		#k<-2
		#m<-5
		#print("f tilt k")
		#print(th$fTilt[k,i])
		#print("f tilt m")
		#print(th$fTilt[m,i])
		#print("f tilt sw k")
		#print(th$fTiltSW[k,i])
		#print("f tilt sw m")
		#print(th$fTiltSW[m,i])
		#info.log2.t1 <- Reduce('+', subjs.info.t1)
		#info.log2.t2 <- Reduce('+', subjs.info.t2)
		#info.log2.t3 <- Reduce('+', subjs.info.t3)
		#info.log2 <- info.log2.t1-info.log2.t2-info.log2.t3
		#print("subject f0 info")
		#print(info.log2)
		
		

		
	}
	if (effInfo==TRUE){
		if(is.null(sampprobs)){crossinfo.log <- matrix(0,nrow=ncol(x),ncol=nrow(th$fTilt))}
		else{
	q <- th$bPrime2SW/th$bPrime2
	crossinfo.log <-  t((dmudeta*(1/th$bPrime2))*x)%*%(t(smmStar*th$fTiltSW - (smm*th$fTilt*rep(q,each=nrow(smmfTilt))))) # cross info for beta and f0 on log scale
	#crossinfo.log <- matrix(0,nrow=ncol(x),ncol=nrow(th$fTilt))
		}
	
#	crossinfo.log2 <- matrix(data=NA, nrow=nrow(t(x)), ncol=nrow(th$fTilt))  #	code used for debugging cross information calc
	
#	for(k in 1:nrow(th$fTilt)){
#	subj.info <- matrix(data=NA, nrow=nrow(x), ncol=ncol(x))
#	for(i in 1:nrow(x)){
#	subj.info[i,] <- x[i,]*(dmudeta[i]*(1/th$bPrime2[i]))*((spt[k]-th$bPrimeSW[i])*th$fTiltSW[k,i]-(spt[k]-th$bPrime[i])*th$fTilt[k,i]*(th$bPrime2SW[i]/th$bPrime2[i]))	
#}
#	crossinfo.log2[,k] <- colSums(subj.info)	
#}
#	print("matrix calc is")
#	print(crossinfo.log)
#	print("subject calc is")
#	print(crossinfo.log2)	
	}else{crossinfo.log <- matrix(0,nrow=ncol(x),ncol=nrow(th$fTilt))}
	
	if(estprobs==TRUE){
		#betaxi.crossinfo <- (1/fullDataSize)*(t(x)%*%diag(dmudeta*(1/th$bPrime2)))%*%(t(smmStar*th$fTiltSW)/sampprobs) # cross info for beta and sampling probabilities 
		#print("test1")
				#print(betaxi.crossinfo)
		
				#print(dim(th$fTiltSW))
				#print(dim(sampprobs))
		betaxi.crossinfo <- matrix(data=NA, nrow=ncol(x), ncol=length(unique(groups)))
		betaxi.crossinfo.subj <- list(length=nrow(x))
		#print(dim(th$fTiltSW))
		for(i in 1:nrow(x)){
			tmp <- matrix(data=NA, nrow=ncol(x), ncol=length(unique(groups)))
				for(k in 1:ncol(betaxi.crossinfo)){
					indexing <- which(groups ==k)
					tmp[,k] <-x[i,]*(dmudeta[i]*(1/th$bPrime2[i]))*(sum(spt[indexing]*th$fTiltSW[indexing,i])-th$bPrimeSW[i]*sum(th$fTiltSW[indexing,i]))/sampprobs[i,indexing[1]]
				}	
			betaxi.crossinfo.subj[[i]] <- tmp
		}
		betaxi.crossinfo <-  Reduce('+', betaxi.crossinfo.subj) 	
		#print("test2")
		#print(betaxi.crossinfo)
		
		#print(betaxi.crossinfo)
	#	print("t1")
	#	t1 <- diag(rowSums(th$fTiltSW/t(sampprobs)))
	#	print(t1)
	#	print("t2")
	#	t2 <- tcrossprod(th$fTiltSW,th$fTiltSW/t(sampprobs))
	#	print(t2)
	#	print("t3")
	#	t3 <- tcrossprod((smmfTilt*(smmStar/th$bPrime)), th$fTiltSW/t(sampprobs))
	#	print(t3)
		
		#f0xi.crossinfo.log <- (1/fullDataSize)*(diag(rowSums(th$fTiltSW/t(sampprobs)))-tcrossprod(th$fTiltSW,th$fTiltSW/t(sampprobs))-
		#tcrossprod((smmfTilt*(smmStar/th$bPrime)), th$fTiltSW/t(sampprobs)))
		#print(f0xi.crossinfo.log)
		
		f0xi.crossinfo.log <- matrix(data=NA, nrow=nrow(th$fTiltSW), ncol=length(unique(groups)))
		f0xi.crossinfo.log.subj <- list(length=nrow(x))
		#print(dim(th$fTiltSW))
		for(i in 1:nrow(x)){
			tmp <- matrix(data=NA, nrow=nrow(th$fTiltSW), ncol=length(unique(groups)))
			for(j in 1:nrow(f0xi.crossinfo.log)){
				for(k in 1:ncol(f0xi.crossinfo.log)){
					indexing <- which(groups ==k)
					tmp1 <- ifelse(j%in%indexing,th$fTiltSW[j,i],0)
					tmp2 <- th$fTiltSW[j,i]*sum(th$fTiltSW[indexing,i])
					tmp3 <- ((spt[j]-th$bPrime[i])*th$fTilt[j,i]*sum((spt[indexing]-th$bPrimeSW[i])*th$fTiltSW[indexing,i]))/(th$bPrime2[i])
					tmp[j,k] <- (tmp1-tmp2-tmp3)/sampprobs[i,indexing[1]]
					
				}
			}
			f0xi.crossinfo.log.subj[[i]] <- tmp
		}
		f0xi.crossinfo.log <-  Reduce('+', f0xi.crossinfo.log.subj)
		#print(f0xi.crossinfo.log)
		
		
	}
	else{
		betaxi.crossinfo<- NULL
		f0xi.crossinfo.log <- NULL
		
	}
	
    list(f0=f0, llik=llik, th=th, conv=conv, iter=iter, nhalf=nhalf,
         score.log=score.log, info.log=info.log, crossinfo.log=crossinfo.log, betaxi.crossinfo=betaxi.crossinfo, f0xi.crossinfo.log=f0xi.crossinfo.log)
}			  
				  
				  
				  
				  
				  
				  
				  
				  
 	
