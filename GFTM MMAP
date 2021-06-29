## #~#~#~#~#~#~# ##
##   GFTM MMAP   ##
## #~#~#~#~#~#~# ##

gftm_mmap <- function(resp, tlid, model="GFTM", mfit=FALSE,
                      prior_i=list(mu_a=1.0, sd_a=0.3),
                      nquad=NULL, qtrunc=c(-4.0, 4.0),
                      itrunc=matrix(c(0.2, 3, -4, 4, -4, 4), nrow=3, ncol=2, byrow=T,
                                    dimnames = list(c("a","g","b"), c("lb","ub"))),
                      bscut=NULL, maxitr=list(em=100, nr=100), tol=list(em=.001, nr=.001, loglike=1e-4)){


  ### #~#~#~#~#~#~#~#~#~#~# ###
  ###   Set up variables    ###
  ### #~#~#~#~#~#~#~#~#~#~# ###

  nexaminee <- dim(resp)[1]
  nitem <- dim(resp)[2]
  tlset <- unique(tlid)
  ntl <- length(tlset)
  K_j <- rep(NA, nitem)
  for (j in 1:nitem){ K_j[j] <- max(resp[,j]) }
  maxK <- max(K_j); ncat <- maxK+1
  name_icol <- c("a", "g", paste("b", 1:maxK, sep=""))

  mu_astr <- log( prior_i$mu_a^2 / sqrt( (prior_i$sd_a)^2 + prior_i$mu_a^2) ) ;  # mean of disc* = log a; pars for random gen
  var_astr <- log( (prior_i$sd_a)^2 / prior_i$mu_a^2 + 1) ;

  ## Glean testlet information ------
  tlinf <- vector(mod="list", length=ntl)
  for (s in 1:ntl){# s<-1

    iidx <- which(tlid==tlset[s])
    tlinf[[s]]$iidx <- iidx
    tlinf[[s]]$nitemtl <- length(iidx)
    tlinf[[s]]$K_j <- K_j[iidx]

    ## Item parameters within the testlet
    name_bjk_s <- NULL
    for (j in 1:length(iidx)){# j<-1
      name_bjk_s <- c(name_bjk_s, paste0("b", iidx[j], ",",  1:K_j[iidx[j]]))
    }
    tlinf[[s]]$name_ipar <- c(paste0("g", s), paste0("a", iidx), name_bjk_s)
    tlinf[[s]]$nipar <- length(tlinf[[s]]$name_ipar)

    ## Item Parameters that need be estimated
    # model <- "FTM"
    if (model=="GFTM"){
      tlinf[[s]]$name_iest <- c(paste0("g", s), paste0("a", iidx), name_bjk_s)
    } else if (model=="FTM"){
      tlinf[[s]]$name_iest <- c(paste0("g", s), name_bjk_s)
    } else if (model=="GPCM"){
      tlinf[[s]]$name_iest <- c(paste0("a", iidx), name_bjk_s)
    } else if (model=="PCM"){
      tlinf[[s]]$name_iest <- name_bjk_s
    }

    tlinf[[s]]$nipar_est <- length(tlinf[[s]]$name_iest)

    mu_i <- rep(0, tlinf[[s]]$nipar_est)
    cov_i <- diag(tlinf[[s]]$nipar_est)
    names(mu_i) <- colnames(cov_i) <- rownames(cov_i) <- tlinf[[s]]$name_iest
    if (grepl("G", model, ignore.case=T)){
      mu_i[grep("a", tlinf[[s]]$name_iest, value=T)] <- mu_astr
      diag(cov_i[grep("a", tlinf[[s]]$name_iest), grep("a", tlinf[[s]]$name_iest)]) <- var_astr
    }
    tlinf[[s]]$mu_i <- mu_i
    tlinf[[s]]$cov_i <- cov_i
    tlinf[[s]]$icov_i <- solve(cov_i)
    rm(mu_i, cov_i)
  } # end of s

  eval(parse(text=paste0("name_ipar <- c(", paste0("tlinf[[",1:ntl,"]]$name_ipar", collapse=","), ")") ))
  eval(parse(text=paste0("name_iest <- c(", paste0("tlinf[[",1:ntl,"]]$name_iest", collapse=","), ")") ))
  nipar <- length(name_ipar)
  niparest <- length(name_iest)


  ### #~#~#~#~#~#~#~# ###
  ###   Initialize    ###
  ### #~#~#~#~#~#~#~# ###

  start_time <- Sys.time()

  imat_init <- matrix(NA, nitem, length(name_icol), dimnames=list(NULL, name_icol))
  ivec_init <- rep(NA, nipar); names(ivec_init) <- name_ipar

  totalsc <- rowSums(resp)
  for (s in 1:ntl){# s<-1
    iidx <- tlinf[[s]]$iidx;
    tlsc <- rowSums(as.matrix(resp[,iidx]))
    ivec_init[paste0("g",s)] <- imat_init[iidx,"g"] <- - qnorm(1 - (mean(tlsc)/sum(K_j[iidx])), mean=0, sd=2)
  }

  for (j in 1:nitem){# j<-1
    ## 1) Slope
    if (K_j[j]==1){
      mean_y1 <- mean(totalsc[resp[,j]==1])
      px <- sum(resp[,j]) / nexaminee
      pcor <- (mean_y1 - mean(totalsc)) / (sd(totalsc)) * sqrt( px / (1-px) )
    } else {
      pcor <- cor(resp[,j], rowSums(resp[,-j]))
    }
    ivec_init[paste0("a",j)] <- imat_init[j,"a"] <- sqrt( pcor^2 / (1 - pcor^2) ) # 1.702 *

    ## 2) Location
    prop_j <- as.vector(table(resp[,j])) / nexaminee
    ivec_init[grep(paste0("b",j,","), name_ipar)] <- imat_init[j, paste0("b",1:K_j[j])] <- - qnorm(cumsum(prop_j)[1:K_j[j]], mean=0, sd=1)
  } # end of j

   if (model=="FTM"){
    ivec_init[grep("a", name_ipar)] <- imat_init[,"a"] <- 1
  } else if (model=="GPCM"){
    ivec_init[grep("g", name_ipar)] <- imat_init[,"g"] <- 0
  } else if (model=="PCM"){
    ivec_init[grep("g", name_ipar)] <- imat_init[,"g"] <- 0
    ivec_init[grep("a", name_ipar)] <- imat_init[,"a"] <- 1
  }

  ivecstr_init <- ivec_init
  imatstr_init <- imat_init
  ivecstr_init[grep("a", name_ipar, value=T)] <- imatstr_init[,"a"] <- log(imat_init[,"a"])


  ### #~#~#~#~#~#~#~#~#~#~#~# ###
  ###   EM Design variables   ###
  ### #~#~#~#~#~#~#~#~#~#~#~# ###

  ### Gaussian quadrature nodes
  if (is.numeric(nquad)){
    quad <- seq(qtrunc[1], qtrunc[2], length.out=nquad)
  } else {
    quad <- seq(qtrunc[1], qtrunc[2], by=0.5)
    nquad <- length(quad)
  }
  quadw <- dnorm(quad) # quadrature weight

  ### Response category coefficients
  k_qkj <- array(NA, dim=c(nquad, maxK+1, nitem),
                 dimnames=list(NULL, paste0("k", 0:maxK), paste0("j", 1:nitem)))
  for (j in 1:nitem){# j<-1
    k_qkj[,1:(K_j[j]+1),j] <- matrix(rep(seq(0, K_j[j]), each=nquad),
                                     nquad, K_j[j]+1, dimnames=list(NULL, paste0("k", 0:K_j[j])))
  }

  ### Response incidence matrix
  u_jki <- array(NA, dim=c(nitem, maxK+1, nexaminee),
                 dimnames=list(paste0("j", 1:nitem), paste0("k", 0:maxK),
                               paste0("i", 1:nexaminee)))
  utmp <- matrix(rep(seq(0, maxK), each=nitem), nitem, maxK+1)
  for (i in 1:nexaminee){# i<-2
    u_jki[,,i] <- (matrix(rep(resp[i,], times=maxK+1), nitem, maxK+1) == utmp) * 1
  }

  revcumsum <- function(x){ rev(cumsum(rev(x))) }


  ### #~#~#~#~#~#~# ###
  ###    EM Cycle   ###
  ### #~#~#~#~#~#~# ###

  cyc <- 1
  delta_em <- rep(100, niparest)
  maxavgradius <- 2 # safeguard in NR
  ivec_curr <- ivec_init;
  imat_curr <- imat_init
  ivecstr_curr <- ivecstr_init;

  track_em <- list()
  track_em$loglike <- rep(NA, maxitr$em+1)

  track_em$delta <- matrix(NA, maxitr$em+1, niparest, dimnames=list(NULL, name_iest)) # delta_em
  track_em$delta_summ <- matrix(NA, maxitr$em+1, 2, dimnames=list(NULL, c("avgabs", "maxabs")))
  track_em$iparvec <- matrix(NA, maxitr$em+1, nipar, dimnames=list(NULL, name_ipar))
  track_em$nitr_nr <- matrix(NA, maxitr$em+1, ntl, dimnames=list(NULL, paste0("s",1:ntl)))
  track_em$iparvec[cyc,] <- as.vector(ivec_curr)
  

  ### ========= Start EM ========= ###
  while (cyc < maxitr$em) {

    ### Compute probability matrix for all items and categories
    p_qkj_curr <- icrf_gftm_mmle(quad, imat_curr, K_j)

    ### Compute likelihood matrix
    like_qi <- matrix(1, nquad, nexaminee)
    for (i in 1:nexaminee){# i<-13
      for (j in 1:nitem){# j<-2
        like_qi[,i] <- like_qi[,i] *
          apply(p_qkj_curr[,,j] ^ matrix(rep(u_jki[j,,i], each=nquad), nquad, ncat), 1, prod)
      }
    }

    ### Compute artificial variables
    num <- like_qi * quadw
    denom <- colSums(num)
    track_em$loglike[cyc+1] <- sum(log( colSums(num) )) ## log-like evaluated at nodes

    tau_q <- rowSums( num / matrix( rep(denom, each=nquad), nquad, nexaminee ) )
    zeta_qkj <- array(0, dim=c(nquad, maxK+1, nitem),
                      dimnames=list(NULL, c(paste0("k", 0:maxK)), c(paste0("j", 1:nitem))))
    for (j in 1:nitem){# j<-2
      for (i in 1:nexaminee){# i<-1
        zeta_qkj[,,j] <- zeta_qkj[,,j] +
          num[,i] * matrix(rep(u_jki[j,,i], each=nquad), nquad, maxK+1) / denom[i]
      }
    }

    ### New item parameters to record
    imat_new <- matrix(NA, nitem, length(name_icol), dimnames=list(NULL, name_icol))
    ivec_new <- rep(NA, length(ivec_curr)); names(ivec_new) <- name_ipar
    ivecstr_new <- rep(NA, length(ivecstr_curr)); names(ivecstr_new) <- name_ipar


    ## Newton-Raphson for one testlet -----------
    for (s in 1:ntl){# s<-1

      ### Local variables
      iidx <- tlinf[[s]]$iidx
      nitemtl <- tlinf[[s]]$nitemtl
      K_j_s <- tlinf[[s]]$K_j
      zeta_qkj_s <- zeta_qkj[,,iidx]
      k_qkj_s <- k_qkj[,,iidx]

      name_ipar_s <- tlinf[[s]]$name_ipar
      name_iest_s <- tlinf[[s]]$name_iest
      nipar_s <- tlinf[[s]]$nipar
      niparest_s <- tlinf[[s]]$nipar_est

      ### Current item par in matrix & vector  forms
      imat_curr_s <- imat_curr[iidx, name_icol]
      ivec_curr_s <- ivec_curr[name_ipar_s]
      ivecstr_curr_s <- ivecstr_curr[name_ipar_s]

      mu_i <- tlinf[[s]]$mu_i
      icov_i <- tlinf[[s]]$icov_i

      itr_nr <- 0
      deltastr_nr <- rep(100, niparest_s)
      conv_nr <- rep(100, 3); 
      ## Track whether NR yielded parameter estimates within the appropriate range
      ## 0= converged; 1= Out of tolerable range from the initial value
      ## 2= Out of lower bound; 3= Out of upper bound

      
      ### Begin NR ----------------------
      while (itr_nr < maxitr$nr) {

        itr_nr <- itr_nr + 1
        ## Reset conv to zero
        conv_nr <- rep(0, 3);

        p_qkj_s <- icrf_gftm_mmle(quad, imat_curr_s, K_j_s)

        ### Compute common coefficients (local variables whose size is defined by testlet)
        lam_qkj0 <- lam_qkj1 <- lam_qkj2 <- array(NA, dim=c(nquad, maxK+1, nitemtl),
                                                  dimnames=list(NULL, c(paste0("k", 0:maxK)), c(paste0("j", iidx))))
        for (j in 1:nitemtl){# j<-1
          cidx <- 1:(K_j_s[j]+1)
          ## Each column computes reverse cumulative sum from m to K_j
          lam_qkj0[,cidx,j] <- t(apply(p_qkj_s[,cidx,j], 1, revcumsum))
          lam_qkj1[,cidx,j] <- t(apply(k_qkj_s[,cidx,j]   * p_qkj_s[,cidx,j], 1, revcumsum))
          lam_qkj2[,cidx,j] <- t(apply(k_qkj_s[,cidx,j]^2 * p_qkj_s[,cidx,j], 1, revcumsum))
          # kmat[,,iidx[1]] * pr_curr_s[,,1]
        }

        lam_q0j1 <- array(apply(lam_qkj1[,"k0",], 2, rep, times=maxK+1), dim=c(nquad, maxK+1, nitemtl))
        lam_q0j2 <- array(apply(lam_qkj2[,"k0",], 2, rep, times=maxK+1), dim=c(nquad, maxK+1, nitemtl))


        ### Compute derivatives
        Lamb <- rep(NA, niparest_s);
        Hess <- matrix(0, niparest_s, niparest_s)
        names(Lamb) <- colnames(Hess) <- rownames(Hess) <- name_iest_s

        name_g <- grep("g", name_iest_s, value=T)
        name_a <- grep("a", name_iest_s, value=T)

        ## b
        for (j in 1:nitemtl){ # j<-1
          for (m in 1:K_j_s[j]){ # m<-1

            name_m <- paste0("b",iidx[j],",",m)

            if (m==K_j_s[j]){
              Lamb[name_m] <- sum( zeta_qkj_s[,(m+1),j] - tau_q * lam_qkj0[,m+1,j] ) -
                sum(icov_i[,name_m] * (ivecstr_curr_s[name_iest_s] - mu_i))
            } else {
              Lamb[name_m] <- sum( rowSums(zeta_qkj_s[,(m+1):(K_j_s[j]+1),j]) - tau_q * lam_qkj0[,m+1,j] ) -
                sum(icov_i[,name_m] * (ivecstr_curr_s[name_iest_s] - mu_i))
            }


            ## H_{b_{m\prime} b_{m}}
            for (mp in 1:K_j_s[j]){# mp<-1

              name_mp <- paste0("b",iidx[j],",",mp)
              # ridx_mp <- grep(paste0("b",iidx[j],",",mp), name_iest_s)
              if (mp <= m){
                Hess[name_m, name_mp] <- Hess[name_mp, name_m] <-
                  - sum( lam_qkj0[,m+1,j] * (1 - lam_qkj0[,mp+1,j]) * tau_q ) - icov_i[name_m, name_mp]
              } else {
                Hess[name_m, name_mp] <- Hess[name_mp, name_m] <-
                  - sum( lam_qkj0[,mp+1,j] * (1 - lam_qkj0[,m+1,j]) * tau_q ) - icov_i[name_m, name_mp]
              }
            }

            ## H_{gb}
            # cidx <- grep("g", name_iest_s)
            #if (length(name_g) > 0){
            Hess[name_m, name_g] <- Hess[name_g, name_m] <-
              - sum( (lam_qkj1[,m+1,j] - lam_qkj1[,1,j] * lam_qkj0[,m+1,j]) * rowSums(zeta_qkj_s[,,j]) ) -
              icov_i[name_g, name_m]
            #}

            ## H_{ab}
            # name_aj <- paste0("a",iidx[j])
            name_aj <- grep(paste0("a",iidx[j]), name_iest_s, value=T)
            # if (length(name_aj) > 0) {
            Hess[name_m, name_aj] <- Hess[name_aj, name_m] <-
              - sum( ivec_curr_s[name_aj] * quad * (lam_qkj1[,m+1,j] - lam_qkj0[,m+1,j] * lam_qkj1[,1,j]) * rowSums(zeta_qkj_s[,,j]) ) -
              icov_i[name_aj, name_m]
            #}

          } # end of m
        } # end of j


        ## L_g, H_{gg}
        Lamb[name_g] <- sum( apply((k_qkj_s - lam_q0j1) * zeta_qkj_s, 1, sum, na.rm=T) ) -
          sum(icov_i[,name_g] * (ivecstr_curr_s[name_iest_s] - mu_i))
        Hess[name_g, name_g] <- - sum(apply((lam_q0j2 - lam_q0j1^2) * zeta_qkj_s, 1, sum, na.rm=T)) -
          icov_i[name_g, name_g]

        if (length(name_a) > 0){
          Lamb[name_a] <- ivec_curr_s[name_a] * apply(quad * (k_qkj_s - lam_q0j1) * zeta_qkj_s, 3, sum, na.rm=T) -
            colSums(icov_i[,name_a] * (ivecstr_curr_s[name_iest_s] - mu_i))

          cff <- lam_qkj1[,"k0",] + matrix(rep(ivec_curr_s[name_a], each=nquad), nquad, nitemtl) *
            quad * ( lam_qkj2[,"k0",] - lam_qkj1[,"k0",]^2 )

          diag(Hess[name_a, name_a]) <- ivec_curr_s[name_a] *
            apply(quad * (k_qkj_s - array(apply(cff, 2, rep, ncat), dim=c(nquad, ncat, nitemtl))) *
                    zeta_qkj_s, 3, sum, na.rm=T) - diag(icov_i[name_a, name_a])

          Hess[name_g, name_a] <- Hess[name_a, name_g] <-
            - ivec_curr_s[name_a] * apply(quad * (lam_q0j2 - lam_q0j1^2) * zeta_qkj_s, 3, sum) - icov_i[name_g, name_a]
        }


        deltastr_nr <- Lamb %*% MASS::ginv(Hess) ;

        ### Normalize delta if it exceeds maximum radius.
        if (sqrt(mean(deltastr_nr^2)) > maxavgradius){
          deltastr_nr <- deltastr_nr * maxavgradius / abs(sqrt(mean(deltastr_nr^2)))
        }

        ivecstr_new_s <- ivecstr_curr_s
        ivecstr_new_s[name_iest_s] <- ivecstr_new_s[name_iest_s] - deltastr_nr;
        names(ivecstr_new_s) <- name_ipar_s

        ivec_new_s <- ivecstr_new_s
        ivec_new_s[grep("a", names(ivec_new_s))] <- exp(ivec_new_s[grep("a", names(ivec_new_s))])


        ### Bisection
        if (is.numeric(bscut)){ # bscut <- 1
          lidx <- which(abs(ivec_new_s[name_iest_s] - ivec_init[name_iest_s]) > bscut)
          if (length(lidx) > 0){
            conv_nr[1] <- 1 ## Out of tolerable range from the initial value
            ivec_new_s[lidx] <- (ivec_new_s[lidx] + ivec_init[name_iest_s[lidx]]) / 2
          }
        }

        ### Truncation
        for (p in 1:dim(itrunc)[1]){# p<-2
          partype <- rownames(itrunc)[p]
          if (any(ivec_new_s[grep(partype, name_iest_s, value=T)] < itrunc[partype,"lb"])){
            conv_nr[2] <- 1 ## Out of lower bound
            lidx <- which(ivec_new_s[grep(partype, name_iest_s, value=T)] < itrunc[partype,"lb"])
            ivec_new_s[grep(partype, name_iest_s, value=T)][lidx] <- itrunc[partype,"lb"]
          }

          if (any(ivec_new_s[grep(partype, name_iest_s, value=T)] > itrunc[partype,"ub"])){
            conv_nr[3] <- 1 ## Exceeded upper bound
            lidx <- which(ivec_new_s[grep(partype, name_iest_s, value=T)] > itrunc[partype, "ub"])
            ivec_new_s[grep(partype, name_iest_s, value=T)][lidx] <- itrunc[partype, "ub"]
          }
        }


        ### New item parameters in a matrix form
        imat_new_s <- matrix(NA, nitemtl, maxK+2, dimnames=list(NULL, name_icol))
        imat_new_s[,"g"] <- ivec_new_s[paste0("g",s)]
        imat_new_s[,"a"] <- ivec_new_s[paste0("a",iidx)]
        for (j in 1:nitemtl){# j<-2
          imat_new_s[j, grep("b", name_icol)[1:K_j_s[j]]] <- ivec_new_s[grep(paste0("b",iidx[j]), name_iest_s, value=T)]
        }


        if ( max(abs(deltastr_nr)) < tol$nr && (sum(conv_nr)==0) ) {
          imat_nr_s <- imat_new_s;
          ivec_nr_s <- ivec_new_s;
          ivecstr_nr_s <- ivecstr_new_s
          break;
        }

        imat_curr_s <- imat_new_s
        ivec_curr_s <- ivec_new_s
        ivecstr_curr_s <- ivecstr_new_s

      } # end of while (NR)


      if (itr_nr < maxitr$nr) {
        imat_new[iidx, ] <- imat_nr_s
        ivec_new[name_ipar_s] <- ivec_nr_s
        ivecstr_new[name_ipar_s] <- ivecstr_nr_s
      } else {
        imat_new[iidx, ] <- imat_curr_s
        ivec_new[name_ipar_s] <- ivec_curr_s
        ivecstr_new[name_ipar_s] <- ivecstr_curr_s
      }

      track_em$nitr_nr[cyc+1, s] <- itr_nr

    } # end of s (testlet)


    delta_em <- ivec_new[name_iest] - ivec_curr[name_iest]

    track_em$delta[cyc+1,] <- delta_em
    track_em$delta_summ[cyc+1,] <- c(mean(abs(delta_em)), max(abs(delta_em)))
    track_em$iparvec[cyc+1,] <- ivec_new 

    if ( (max(abs(delta_em)) < tol$em) ||
         (cyc > 2 && all( abs(diff(track_em$loglike[(cyc+1):(cyc-1)])) < tol$loglike ) ) ) {
      imat_est <- imat_new;
      ivec_est <- ivec_new;
      break;
    }

    ## If not converged,
    cyc <- cyc + 1

    imat_curr <- imat_new
    ivec_curr <- ivec_new
    ivecstr_curr <- ivecstr_new

  } # end of cycle


  ### #~#~#~#~#~#~#~#~#~# ###
  ###    Standard Error   ###
  ### #~#~#~#~#~#~#~#~#~# ###


  if (!exists("imat_est")){imat_est <- imat_curr}

  iest <- cbind(tlid, imat_est, K_j); colnames(iest) <- c("s", name_icol, "K_j") # View(iest)
  iest_se <- matrix(NA, nitem, maxK+4, dimnames=list(NULL, colnames(iest)))
  iest_se[,"s"] <- tlid; iest_se[,"K_j"] <- K_j
  iest_cov <- vector(mode="list", ntl)

  ## Re-evaluate Hessian using updated global variables from EM

  ### Compute probability matrix for all items and categories
  p_qkj_est <- icrf_gftm_mmle(quad, imat_est, K_j) # View(imat_est)

  ### Compute likelihood matrix
  like_qi <- matrix(1, nquad, nexaminee) # View(like_qi)
  for (i in 1:nexaminee){# i<-13
    for (j in 1:nitem){# j<-2
      like_qi[,i] <- like_qi[,i] *
        apply(p_qkj_est[,,j] ^ matrix(rep(u_jki[j,,i], each=nquad), nquad, maxK+1), 1, prod)
    }
  }

  ### Compute artificial variables
  num <- like_qi * quadw
  denom <- colSums(num)

  tau_q <- rowSums( num / matrix( rep(denom, each=nquad), nquad, nexaminee ) )
  zeta_qkj <- array(0, dim=c(nquad, maxK+1, nitem),
                    dimnames=list(NULL, c(paste0("k", 0:maxK)), c(paste0("j", 1:nitem))))
  for (j in 1:nitem){# j<-2
    for (i in 1:nexaminee){# i<-1
      zeta_qkj[,,j] <- zeta_qkj[,,j] +
        num[,i] * matrix(rep(u_jki[j,,i], each=nquad), nquad, maxK+1) / denom[i]
    }
  }


  ### Compute common coefficients
  lam_qkj0 <- lam_qkj1 <- lam_qkj2 <- array(NA, dim=c(nquad, maxK+1, nitem),
                                            dimnames=list(NULL, c(paste0("k", 0:maxK)), c(paste0("j", 1:nitem))))
  for (j in 1:nitem){# j<-1
    cidx <- 1:(K_j[j]+1)
    lam_qkj0[,cidx,j] <- t(apply(p_qkj_est[,cidx,j], 1, revcumsum))
    lam_qkj1[,cidx,j] <- t(apply(k_qkj[,cidx,j]   * p_qkj_est[,cidx,j], 1, revcumsum))
    lam_qkj2[,cidx,j] <- t(apply(k_qkj[,cidx,j]^2 * p_qkj_est[,cidx,j], 1, revcumsum))
  }

  ### Compute derivatives
  lam_q0j1 <- array(apply(lam_qkj1[,"k0",], 2, rep, times=maxK+1), dim=c(nquad, maxK+1, nitem))
  lam_q0j2 <- array(apply(lam_qkj2[,"k0",], 2, rep, times=maxK+1), dim=c(nquad, maxK+1, nitem))


  for (s in 1:ntl){# s<-2

    iidx <- tlinf[[s]]$iidx
    nitemtl <- tlinf[[s]]$nitemtl
    nipartl <- tlinf[[s]]$nipar

    name_ipar_s <- tlinf[[s]]$name_ipar
    name_iest_s <- tlinf[[s]]$name_iest
    nipar_s <- tlinf[[s]]$nipar
    niparest_s <- tlinf[[s]]$nipar_est


    ## Hessian matrix -------------
    Hess_s <- matrix(0, niparest_s, niparest_s, dimnames=list(name_iest_s, name_iest_s))

    name_g <- grep("g", name_iest_s, value=T)
    Hess_s[name_g, name_g] <-
      - sum(apply((lam_q0j2[,,iidx] - lam_q0j1[,,iidx]^2) * zeta_qkj[,,iidx], 3, sum, na.rm=T))

    name_a <- grep("a", name_iest_s, value=T)
    if (length(name_a)>0){
      diag(Hess_s[name_a, name_a]) <- - apply(quad^2 * (lam_q0j2[,,iidx] - lam_q0j1[,,iidx]^2) * zeta_qkj[,,iidx], 3, sum, na.rm=T)
      Hess_s[name_g, name_a] <- Hess_s[name_a, name_g] <-
        - apply(quad * (lam_q0j2[,,iidx] - lam_q0j1[,,iidx]^2) * zeta_qkj[,,iidx], 3, sum)
    }


    for (j in 1:nitemtl){# j<-1
      for (m in 1:K_j[iidx[j]]){ # m<-1

        name_m <- paste0("b", iidx[j],",", m)

        ## H_{b_{mp} b_{m}}
        for (mp in 1:K_j[iidx[j]]){
          name_mp <- paste0("b",iidx[j],",",mp)
          if (mp <= m){
            Hess_s[name_m, name_mp] <- Hess_s[name_mp, name_m] <- - sum( lam_qkj0[,m+1,iidx[j]] * (1 - lam_qkj0[,mp+1,iidx[j]]) * tau_q )
          } else {
            Hess_s[name_m, name_mp] <- Hess_s[name_mp, name_m] <- - sum( lam_qkj0[,mp+1,iidx[j]] * (1 - lam_qkj0[,m+1,iidx[j]]) * tau_q )
          }
        }

        ## H_{gb}
        Hess_s[name_g, name_m] <- Hess_s[name_m, name_g] <-
          - sum( (lam_qkj1[,m+1,iidx[j]] - lam_qkj1[,1,iidx[j]] * lam_qkj0[,m+1,iidx[j]]) * rowSums(zeta_qkj[,,iidx[j]]) )

        ## H_{ab}
        name_aj <- grep(paste0("a",iidx[j]), name_iest_s, value=T)
        Hess_s[name_aj, name_m] <- Hess_s[name_m, name_aj] <- - sum( quad * (lam_qkj1[,m+1,iidx[j]] - lam_qkj0[,m+1,iidx[j]] * lam_qkj1[,1,iidx[j]]) * rowSums(zeta_qkj[,,iidx[j]]) )

      }
    }

    covmat <- MASS::ginv(-Hess_s);
    colnames(covmat) <- rownames(covmat) <- name_iest_s

    iest_cov[[s]] <- covmat
    iest_se_vec <- sqrt( diag(covmat) );
    iest_se[iidx,"g"] <- iest_se_vec[paste0("g",s)]
    iest_se[iidx,"a"] <- iest_se_vec[paste0("a",iidx)]
    for (j in 1:nitemtl){# j<-4
      iest_se[iidx[j], paste0("b", 1:K_j[iidx[j]])] <- iest_se_vec[grep(paste0("b",iidx[j]), name_iest_s)]
    }

  } # end of s (testlet)


  ### #~#~#~#~#~#~#~#~#~# ###
  ###    PPAR Variance    ###
  ### #~#~#~#~#~#~#~#~#~# ###

  ## 'denom' is person-specific; the common denominator across quadratues for the given person
  post_qi <- num / matrix( rep(denom, each=nquad), nquad, nexaminee )

  if (length(grep("G", model))==0){
    var_p <- mean( colSums(quad^2 * post_qi) / colSums(post_qi) )
    npar <- niparest + 1
  } else {
    var_p <- 1
    npar <- niparest
  }

  end_time <- Sys.time()
  est_time <- end_time - start_time


  ### #~#~#~#~#~#~#~#~#~# ###
  ###    Export results   ###
  ### #~#~#~#~#~#~#~#~#~# ###

  if(isTRUE(mfit)){

    loglike <- sum(log( colSums( like_qi * dnorm(quad, mean=0, sd=sqrt(var_p)) ) ))

    relfit <- rep(NA, 8)
    names(relfit) <- c("deviance", "loglike", "AIC", "CAIC", "BIC", "aBIC", "ncyc", "npar")

    relfit[c("deviance", "loglike")] <- c(-2 * loglike, loglike)
    relfit[c("AIC", "CAIC", "BIC", "aBIC")] <- -2*loglike + c(2*npar,
                                                              npar * (log(nexaminee) + 1),
                                                              npar * log(nexaminee),
                                                              npar * log((nexaminee+2)/24))
    relfit[c("ncyc", "npar")] <- c(cyc, npar)
  } else {
    relfit=NULL  }

  output <- list(iest=iest, ise=iest_se, icov=iest_cov, pvar=var_p,
                 mfit=relfit,
                 track=list(
                   cyc=cyc, #niparest=niparest,
                   loglike=track_em$loglike[1:(cyc+1)],
                   delta=track_em$delta[1:(cyc+1),],
                   delta_summ=track_em$delta_summ[1:(cyc+1),],
                   iparvec=track_em$iparvec[1:(cyc+1),],
                   nitr_nr=track_em$nitr_nr[1:(cyc+1),],
                   npar=npar),
                 supp=list(quad=quad,
                           p_jkq=aperm(p_qkj_est, c(3,2,1)),
                           like_iq=t(like_qi),
                           post_iq=t(post_qi),
                           tau_q = tau_q,
                           zeta_qkj=zeta_qkj,
                           est_time=est_time),
                 input=list(K_j=K_j, resp=as.matrix(resp), u_jki=u_jki)
  )
  return(output)

} # end of function
