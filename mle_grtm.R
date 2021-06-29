## #~#~#~#~#~#~# ##
##   GRTM::MLE   ##
## #~#~#~#~#~#~# ##

mle_grtm <- function(resp, ipar, tlid, SE=FALSE, DeBug=FALSE,
                     D=1, bscut=NULL, maxitr=200,
                     tol=list(nr=c(1.0e-4, 1.0e-3), loglike=1.0e-6),
                     trunc=list(th=c(-3.5, 3.5), gam=c(-3.0, 3.0))){
  ## Estimate MLE under polytomous random-effect testlet model
  
  
  ## #~#~#~#~#~#~#~#~#~#~#~# ##
  ##   Setting up Variables  ##
  ## #~#~#~#~#~#~#~#~#~#~#~# ##
  
  ## if one person estimation
  if (!is.matrix(resp)){resp <- t(as.matrix(resp))}
  
  nexaminee <- nrow(resp)
  nitem <- ncol(resp)
  tlidu <- unique(tlid) 
  ntl <- length(tlidu)
  
  disc <- ipar[,grep("a", colnames(ipar), value=T)]
  cloc  <- ipar[,grep("b", colnames(ipar), value=T)]
  
  maxK <- dim(cloc)[2]
  if (any(grepl("K_j", colnames(ipar)))){
    K_j <- ipar[,"K_j"]
  } else {
    K_j <- rowSums(1*(!is.na(cloc)))
  }
  
  kmat <- matrix(rep(seq(0, maxK), each=nitem), nitem, maxK+1)
  u_jki <- array(NA, dim=c(nitem, maxK+1, nexaminee),
                 dimnames=list(paste0("j", 1:nitem), paste0("k", 0:maxK),
                               paste0("i", 1:nexaminee)))
  for (i in 1:nexaminee){# i<-2
    u_jki[,,i] <- (matrix(rep(resp[i,], times=maxK+1), nitem, maxK+1)==kmat) * 1
  }
  
  
  
  ## #~#~#~#~#~#~#~# ##
  ##   Initialize    ##
  ## #~#~#~#~#~#~#~# ##  
  
  th_init <- qnorm(rowSums(resp)/sum(K_j))
  th_init[rowSums(resp)==sum(K_j)] <- trunc$th[2]
  th_init[rowSums(resp)==0] <- trunc$th[1]
  
  gam_init <- matrix(NA, nexaminee, ntl) # View(gam_init)
  for (s in 1:ntl){# s<-1
    iset <- which(tlid==tlidu[s])
    tlsc <- rowSums(resp[,iset])
    gam_init[,s] <- qnorm(tlsc / sum(K_j[iset]), mean=0, sd=0.5)
    ## Note. Initial gamma is calculated assuming SD of 0.5.
    gam_init[tlsc==sum(K_j[iset]),s] <- trunc$gam[2]
    gam_init[tlsc==0,s] <- trunc$gam[1]
    
  }
  ## Note. Initial theta & gamma can differ in a single testlet because of different priors imposed.
  ## Note. Initial values had little impact on final MLE. No matter what values the NR started, 
  ##       converged to the same value.  
  
  
  ## #~#~#~#~#~#~#~#~#~# ##
  ##    Newton Raphson   ##
  ## #~#~#~#~#~#~#~#~#~# ##
  
  ## Stopping criteria  
  maxradius <- 2     
  track_itr <- list()
  track_itr$loglike <- matrix(NA, nexaminee, maxitr+1)
  track_itr$loglike[,1] <- 100  
  track_itr$nitr <- rep(NA, nexaminee)
  
  pest <- matrix(NA, nexaminee, ntl+1, dimnames=list(NULL, c("th", paste0("gam" , seq(1:ntl)))))
  pest_nonconv <- rep(NA, nexaminee)
  
  
  ### | Start of for-loop ------------------
  
  for (i in 1:nexaminee){# i <- 1
    
    itr <- 0
    th_curr <- th_init[i]
    gam_curr <- gam_init[i,]    

    while (itr < maxitr){
      
      itr <- itr + 1
      
      if (trunc$th[1] <= th_curr && th_curr <= trunc$th[2]){
        th_try <- th_curr
      } else if (th_curr <= trunc$th[1]){
        th_try <- (trunc$th[1] + th_curr) / 2
      } else if (th_curr >= trunc$th[2]){
        th_try <- (trunc$th[2] + th_curr) / 2
      }
      
      gam_try <- rep(NA, ntl)
      for (s in 1:ntl){ # s<-1
        if (trunc$gam[1] <= gam_curr[s] && gam_curr[s] <= trunc$gam[2]){
          gam_try[s] <- gam_curr[s]
        } else if (gam_curr[s] <= trunc$gam[1]){
          gam_try[s] <- (trunc$gam[1] + gam_curr[s]) / 2
        } else if (gam_curr[s] >= trunc$gam[2]){
          gam_try[s] <- (trunc$gam[2] + gam_curr[s]) / 2
        }
      }
      
      irfs <- icrf_grtm(c(th_try, gam_try), ipar)
      
      ## Derivatives
      lam_j1 <- rowSums(kmat * irfs, na.rm=T) # [nitem x 1]
      lam_j2 <- rowSums(kmat^2 * irfs, na.rm=T)
      
      Lambda  <- matrix(NA, 1, (ntl+1))
      Hessian <- matrix(0, (ntl+1), (ntl+1))
      
      Lambda[1]    <-   sum(D * disc[,1] * u_jki[,,i] * (kmat - lam_j1))
      Hessian[1,1] <- - sum(D^2 * disc[,1]^2 * u_jki[,,i] *(lam_j2 - lam_j1^2))
      for (s in 1:ntl){ # s<-1
        iset <- which(tlid==tlidu[s])
        Lambda[s+1] <- sum(D * disc[iset, 2] * u_jki[iset,,i] * (kmat[iset,] - lam_j1[iset])) 
        Hessian[s+1, s+1] <- - sum(D^2 * disc[iset, 2]^2 * u_jki[iset,,i] * (lam_j2[iset] - lam_j1[iset]^2))
        Hessian[1, s+1] <- Hessian[s+1, 1] <- - sum(D^2 * disc[iset, 1] * disc[iset, 2] * u_jki[iset,,i] 
                                                    * (lam_j2[iset] - lam_j1[iset]^2))
        ## Note. Cross second-order derivatives bw gammas equal 0.
      }
      
      delta <- Lambda %*% MASS::ginv(Hessian)
      
      ## Normalize delta if it exceeds maximum radius.
      if (sqrt(sum(delta^2)) > maxradius){
        delta <- delta * maxradius / abs(sqrt(sum(delta^2)))
      }
      
      th_new  <- th_try - delta[1]
      gam_new <- gam_try - delta[2:(ntl+1)]
      
      track_itr$loglike[i,itr+1] <- sum(u_jki[,,i] * log(irfs), na.rm=T)      
      
      if ( ( abs(delta[1]) < tol$nr[1] && max(abs(delta[2:(ntl+1)])) < tol$nr[2] ) ||
           ( itr > 1 && all(abs(diff(track_itr$loglike[i,(itr-1):(itr+1)])) < tol$loglike ) ) ){  
        ## Convergence: If differences in the estimated values are sufficiently small OR
        ## the log-like is stably small 
        
        th_hat <- th_new
        gam_hat <- gam_new
        break ;
        
      } else {
        
        if (is.numeric(bscut)){
          ## Bisection
          if (abs(th_new - th_init) > bscut){
            th_new <- (th_new + th_init) / 2
          }
          
          if (any( abs(gam_new - gam_init) > bscut)){
            idx <- which(abs(gam_new - gam_init) > bscut)
            gam_new[idx] <- (gam_new[idx] + gam_init[idx]) / 2
          }
          ## Note. Bissection can lead to circulation of the same values
        }
        
        th_curr <- th_new
        gam_curr <- gam_new       

      }   
      
    } # end of while
    
    
    if (itr < maxitr){
      ## If converged within the tolerance criteria
      
      if (abs(th_hat) > (trunc$th[2] + 0.5)){
        pest_nonconv[i] <- 2 # out of bounds
        th_hat[th_hat < (trunc$th[1] - 0.5)] <- trunc$th[1] - 0.5
        th_hat[th_hat > (trunc$th[2] + 0.5)] <- trunc$th[2] + 0.5
      }
      pest[i,] <- c(th_hat, gam_hat)
      
    } else {
      ## If not converged
      
      pest_nonconv[i] <- 1
      th_curr[th_curr < trunc$th[1] - 0.5] <- trunc$th[1] - 0.5
      th_curr[th_curr > trunc$th[2] + 0.5] <- trunc$th[2] + 0.5
      gam_curr[gam_curr < trunc$gam[1] - 0.5] <- trunc$gam[1] - 0.5
      gam_curr[gam_curr > trunc$gam[2] + 0.5] <- trunc$gam[2] + 0.5
      pest[i,] <- c(th_curr, gam_curr)
    }
    
    track_itr$nitr[i] <- itr
  } # end of i (examinee)
  
  ### End of for-loop ------------------ |
  
  
  
  if (isTRUE(SE)){
    
    pest_se <- matrix(NA, nexaminee, ntl+1, dimnames=list(NULL, c("th", paste0("gam" , seq(1:ntl)))))
    for (i in 1:nexaminee){# i <- 1
      irfs <- icrf_grtm(pest[i,], ipar)
      
      lam_j1 <- rowSums(kmat * irfs, na.rm=T)
      lam_j2 <- rowSums(kmat^2 * irfs, na.rm=T)
      
      infmat <- matrix(0, (ntl+1), (ntl+1))
      infmat[1,1] <- sum(D^2 * disc[,1]^2 * (lam_j2 - lam_j1^2))
      for (s in 1:ntl){ # s<-1
        iset <- which(tlid==tlidu[s])
        infmat[s+1, s+1] <- sum(D^2 * disc[iset, 2]^2 * (lam_j2[iset] - lam_j1[iset]^2))
        infmat[1, s+1] <- infmat[s+1, 1] <- sum(D^2 * disc[iset, 1] * disc[iset, 2] * (lam_j2[iset] - lam_j1[iset]^2))
      }
      
      pest_se[i, ] <- sqrt( diag(MASS::ginv(infmat)) )
    }
    
  } else {
    pest_se=NULL
  }
  
  pest_init <- cbind(th_init, gam_init)
  colnames(pest_init) <- colnames(pest)
  
  output <- list(est=pest, 
                 se=pest_se,
                 nonconv=pest_nonconv,
                 init=pest_init,
                 track=track_itr)
  
  return(output)
  
} # end of function

