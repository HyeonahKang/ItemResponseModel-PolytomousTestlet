## #~#~#~#~#~#~#~# ##
##    GFTM-MMLE    ##
## #~#~#~#~#~#~#~# ##

icrf_gftm_mmle = function(ppar, ipar, K_j=NULL, model="GFTM"){ 
  ## Computes item category response function for multiple people to multiple items
  ## in an array under the fixed-effect testlet model
  
  nexaminee <- length(ppar)
  
  if (!is.matrix(ipar)){ipar <- t(as.matrix(ipar))}
  
  nitem <- dim(ipar)[1]
  if (is.null(K_j)){
    if (length(grep("K_j", colnames(ipar))) > 0){
      K_j <- ipar[,"K_j"]
    } else {
      K_j <- rowSums(1*!is.na(ipar[,grep("b", colnames(ipar))]))
    }
  }
  maxK <- max(K_j)
  cloc <- subset(ipar, select=grep("b", colnames(ipar), value=T))
  
  if (grepl("FTM", model, ignore.case=T)){
    ikernel <- function(ppar, j, k){exp( k * (ipar[j,"a"] * ppar + ipar[j,"g"]) + sum(cloc[j, 1:k]) ) }
  } else if (grepl("PCM", model, ignore.case=T)){
    ikernel <- function(ppar, j, k){exp( k * ipar[j,"a"] * ppar + sum(cloc[j, 1:k]) )}
  }
  
  pr <- array(NA, dim=c(nexaminee, maxK+1, nitem),
              dimnames=list(NULL, c(paste0("k", 0:maxK)), c(paste0("j", 1:nitem))))
  for (j in 1:nitem){# j<-1
    exps_k <- matrix(0, nexaminee, K_j[j]+1)
    exps_k[,1] <- exp(0)
    for (k in 1:K_j[j]){ # k<-1
      exps_k[,k+1] <- ikernel(ppar, j, k)
    }
    
    pr[,1:(K_j[j]+1),j] <- exps_k / rowSums(exps_k)
  } # end of j
  
  if (nitem==1){pr<-pr[,,1]}
  
  return(pr)
}


