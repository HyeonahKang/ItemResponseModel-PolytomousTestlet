
icrf_gftm = function(ppar, ipar, model="GFTM"){ 
  ## Computes item category response function for multiple people to multiple items
  ## in an array under the fixed-effect testlet model  
  
  # if (!is.matrix(ppar)){ ppar <- t(as.matrix(ppar)) }
  if (!is.matrix(ipar)){ ipar <- t(as.matrix(ipar)) }
  
  nexaminee <- length(ppar)
  nitem <- nrow(ipar)
  
  cloc <- ipar[,grep("b", colnames(ipar), value=T)]
  if (any(grepl("K_j", colnames(ipar)))){
    ## if 'ipar' contains K_j column
    K_j <- ipar[,"K_j"]
  } else {
    K_j <- rowSums(1*!is.na(cloc))
  }
  maxK <- max(K_j)
  
  if (grepl("FTM", model, ignore.case=T)){
    ikernel <- function(ppar, j, k){exp( k * (ipar[j,"a"] * ppar + ipar[j,"g"]) + sum(cloc[j, 1:k]) ) }
  } else if (grepl("PCM", model, ignore.case=T)){
    ikernel <- function(ppar, j, k){exp( k * ipar[j,"a"] * ppar + sum(cloc[j, 1:k]) )}
  }

  pr <- array(NA, dim=c(nexaminee, maxK+1, nitem), 
              dimnames=list(NULL, paste0("k", 0:maxK), paste0("j", 1:nitem)))
  for (j in 1:nitem){# j<-3
    exps_k <- matrix(0, nexaminee, K_j[j]+1)
    exps_k[,1] <- exp(0)
    for (k in 1:K_j[j]){ # k<-1
      exps_k[,k+1] <- ikernel(ppar, j, k)
    }
    pr[,1:(K_j[j]+1),j] <- exps_k / rowSums(exps_k)
  } # end of j
  
  # if (isFALSE(ikj)){
  pr <- aperm(pr, c(3,2,1))
  # }
  if (nexaminee==1){pr <- pr[,,1]}
  
  return(pr)
  
}

