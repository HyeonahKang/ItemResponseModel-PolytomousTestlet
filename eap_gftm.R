## #~#~#~#~#~#~# ##
##   GFTM::EAP   ##
## #~#~#~#~#~#~# ##

eap_gftm <- function(resp, ipar, model="GFTM", var_p=NULL,
                     SE=FALSE, trunc=c(-4.5, 4.5), nquad=40){
  
  ## Computes EAP for multiple people & multiple items under FETM
  ## 'ipar': [s, a, gam, b1, ... , bK, K_j]
  
  
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
    K_j <- rowSums(1*(!is.na(ipar[,grep("b", colnames(ipar), value=T)])))
  }
  maxK <- max(K_j)
  
  kmat <- matrix(rep(seq(0, maxK), each=nitem), nitem, maxK+1)
  u_jki <- array(NA, dim=c(nitem, maxK+1, nexaminee),
                 dimnames=list(paste0("j", 1:nitem), paste0("k", 0:maxK),
                               paste0("i", 1:nexaminee)))
  for (i in 1:nexaminee){# i<-2
    u_jki[,,i] <- (matrix(rep(resp[i,], times=maxK+1), nitem, maxK+1)==kmat) * 1
  }
  
  quad <- seq(trunc[1], trunc[2], length.out=nquad)
  mu_p <- 0
  if (is.null(var_p)){var_p <- 1}
  weight <- dnorm(quad, mu_p, sqrt(var_p))
  
  
  ## #~#~#~#~#~# ##
  ##     EAP     ##
  ## #~#~#~#~#~# ##
  
  irfs <- icrf_gftm(quad, ipar, model=model)
  pest <- pest_se <- matrix(NA, nexaminee, 1)
  for (i in 1:nexaminee){
    umat <- (matrix(rep(resp[i,], times=maxK+1), nitem, maxK+1)==kmat) * 1
    like_quad <- apply(irfs^array(rep(umat, nquad), c(nitem, maxK+1, nquad)), 3, prod)
    
    pest[i] <- sum(quad * like_quad * weight) / sum(like_quad * weight)
    
    if (isTRUE(SE)){
      ## Compute conditional SE
      pest_se[i] <- sqrt( sum((quad - pest[i])^2 * like_quad * weight) / sum(like_quad * weight) )
    }
  }
  
  return(list(est=pest, se=pest_se))
  
}

