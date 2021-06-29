## #~#~#~#~#~#~# ##
##   GFTM::MLE   ##
## #~#~#~#~#~#~# ##


mle_gftm <- function(resp, ipar, model="GFTM", SE=FALSE, DeBug=FALSE,
                     D=1, bscut=NULL, maxitr=200, 
                     tol=list(nr=1.0e-4, loglike=1.0e-6),
                     trunc=c(-3.5, 3.5)){
  
  ## #~#~#~#~#~#~#~#~#~#~#~# ##
  ##   Setting up Variables  ##
  ## #~#~#~#~#~#~#~#~#~#~#~# ##
  
  ## if one person estimation
  if (!is.matrix(resp)){resp <- t(as.matrix(resp))}
  
  nexaminee <- nrow(resp)
  nitem <- ncol(resp) 
  
  if (any(grepl("K_j", colnames(ipar)))){
    K_j <- ipar[,"K_j"]
  } else {
    K_j <- rowSums(1*(!is.na(cloc)))
  }
  maxK <- max(K_j)
  cloc <- ipar[,grep("b", colnames(ipar), value=T)[1:maxK]]
  
  
  kmat <- matrix(rep(seq(0, maxK), each=nitem), nitem, maxK+1)
  u_jki <- array(NA, dim=c(nitem, maxK+1, nexaminee))
  for (i in 1:nexaminee){# i<-2
    u_jki[,,i] <- (matrix(rep(resp[i,], times=maxK+1), nitem, maxK+1)==kmat) * 1
  }
  
  
  ## #~#~#~#~#~#~#~# ##
  ##   Initialize    ##
  ## #~#~#~#~#~#~#~# ##  
  
  th_init <- qnorm(rowSums(resp)/sum(K_j))
  th_init[rowSums(resp)==sum(K_j)] <- trunc[2]
  th_init[rowSums(resp)==0] <- trunc[1]
  
  
  ## #~#~#~#~#~#~#~#~#~# ##
  ##    Newton Raphson   ##
  ## #~#~#~#~#~#~#~#~#~# ##
  
  maxradius <- 1

  track_itr <- list()
  track_itr$loglike <- matrix(NA, nexaminee, maxitr+1)
  track_itr$loglike[,1] <- 100
  track_itr$nitr <- rep(NA, nexaminee)
  
  pest <- matrix(NA, nexaminee, 1)
  pest_nonconv <- rep(NA, nexaminee)
  
  ### | Start of for-loop ------------------
  
  for (i in 1:nexaminee){# i<-1
    
    itr <- 0
    th_curr <- th_init[i]
    
    
    while (itr < maxitr){
      
      itr <- itr + 1
      
      if (trunc[1] <= th_curr && th_curr <= trunc[2]){
        th_try <- th_curr
      } else if (th_curr <= trunc[1]){
        th_try <- (trunc[1] + th_curr) / 2
      } else if (th_curr >= trunc[2]){
        th_try <- (trunc[2] + th_curr) / 2
      }
      
      
      irfs <- icrf_gftm(th_try, ipar, model=model)
      
      ## Derivatives
      lam_j1 <- rowSums(kmat * irfs, na.rm=T)
      lam_j2 <- rowSums(kmat^2 * irfs, na.rm=T)
      
      dlogL <- sum(D * ipar[,"a"] * u_jki[,,i] * (kmat - lam_j1))
      d2logL <- - sum(D^2 * ipar[,"a"]^2 * u_jki[,,i] * (lam_j2 - lam_j1^2))
      
      delta <- dlogL / d2logL
      
      ### Normalize delta if it exceeds maximum radius.
      if (sqrt(sum(delta^2)) > maxradius){
        delta <- delta * maxradius / abs(sqrt(sum(delta^2)))
      }
      
      th_new <- th_try - delta
      
      track_itr$loglike[i,itr+1] <- sum(u_jki[,,i] * log(irfs), na.rm=T)      
      
      if ( (abs(delta) < tol$nr) ||
           ( itr > 1 && all(abs(diff(track_itr$loglike[i, (itr-1):(itr+1)])) < tol$loglike) ) ){ 
           
        th_hat <- th_new
        break ;
        
      } else {
        
        if (is.numeric(bscut)){
          if (abs(th_new - th_init) > bscut){
            th_new <- (th_new + th_init) / 2
          }
        }
        
        th_curr <- th_new
      }
      
      rm(th_new)
      
    } # end of while
    
    if (itr < maxitr){
      ## If converged within the tolerance criteria
      if (abs(th_hat) > (trunc[2] + 0.5)){
        pest_nonconv[i] <- 2 # out of bounds
        th_hat[th_hat < (trunc[1] - 0.5)] <- trunc[1] - 0.5
        th_hat[th_hat > (trunc[2] + 0.5)] <- trunc[2] + 0.5
      }
      pest[i] <- th_hat
      
    } else {
      ## If not converged
      pest_nonconv[i] <- 1
      th_curr[th_curr < (trunc[1] - 0.5)] <- trunc[1] - 0.5
      th_curr[th_curr > (trunc[2] + 0.5)] <- trunc[2] + 0.5
      pest[i] <- th_curr
    }
    
    track_itr$nitr[i] <- itr
    
  } # end of i (examinee)
  
  ### End of for-loop ------------------ |
  
  
  if (isTRUE(SE)){
    pest_se <- matrix(NA, nexaminee, 1)
    for (i in 1:nexaminee){# i <- 1
      irfs <- icrf_gftm(pest[i], ipar, model=model)
      lam_j1 <- rowSums(kmat * irfs, na.rm=T)
      lam_j2 <- rowSums(kmat^2 * irfs, na.rm=T)
      inf <- sum( D^2 * ipar[,"a"]^2 * (lam_j2 - lam_j1^2) )
      pest_se[i] <- 1 / sqrt(inf)
    }
  } 
  
  output <- list(est=pest, 
                 se=pest_se,
                 nonconv=pest_nonconv,
                 init=as.matrix(th_init),
                 track=track_itr)
  
  return(output)
  
} # end of function

