generate_obs <- function(
    seed = 2023,
    ndays = 123,
    max.day.birth = 100,
    obsdays = seq(35, 90, by = 12),
    daysPerBin = 4,
    obsdays_bin = obsdays / daysPerBin,
    nIntervals_s =  floor(ndays / daysPerBin),
    nIntervals_c = 25,
    cohort_days = (1:nIntervals_c) * daysPerBin - 2,
    nborn = 500,
    mu.leave = 31.5,
    sd.leave = 7,
    mu.moult = 23,
    sd.moult = 5,
    skew.bday = 2,
    mu.bday = 34,
    sd.bday = 12,
    survInt = 4,
    survSlope= 0.1,
    xchi = 4,
    birth.type = c("skewnorm", "dirichlet"),
    obs.type= c('poisson', 'mvn'),
    pow = 0.95,
    pom = 0.95,
    pod= 0.95,
    pcm = 0.91,
    pcw = 1,
    use_pcw = T,  #include or exclude correct classification of whites 
    add_dead = TRUE
) {
  
  library(MASS)
  library(MCMCpack)
  library(sn)
  library(RColorBrewer)
  
  set.seed(seed)
  birth.type <- match.arg(birth.type)
  obs.type <- match.arg(obs.type)
  #------------------------------
  # Birth distribution
  #------------------------------
  if (birth.type == "skewnorm") {
    
    skewnormBDFunction <- function(mu, sigma, alpha, nsteps, days) {
      xvec <- numeric(nsteps)
      for (i in 1:nsteps) {
        xvec[i] <- exp(dsn(days[i], xi = mu, omega = sigma, alpha = alpha, log = TRUE))
      }
      return(xvec / sum(xvec))
    }
    
    Pb <- skewnormBDFunction(mu.bday, sd.bday, skew.bday, nIntervals_c, cohort_days)
    cohort.size <- nborn*Pb
    #cohort.size <- rmultinom(1, nborn, Pb)[, 1]
    #plot(cohort.size)
  } else if (birth.type == "dirichlet") {
    days.per.bin <- c(rep(4, 24), 5)
    rate <- 10
    nintervals <- length(days.per.bin)
    
    alpha.bday <- sapply(1:nintervals, function(i) {
      shape <- dnorm(i * days.per.bin[i], 45, 15) * 100 + 0.1
      rgamma(1, shape * rate, rate)
    })
    
    pb_binned <- rdirichlet(1, alpha.bday)[1, ]
    Pb <- pb_binned / sum(pb_binned)
    cohort.size <- nborn*Pb
    #cohort.size <- rmultinom(1, nborn, Pb)[, 1]
  }
  
  
  #cohort.size <- rmultinom(1, nborn, Pb)[, 1]
  
  #------------------------------
  # survival and state transitions
  #------------------------------
  hazFunction <- function(mu, std, nsteps, dstep) {
    hazprob <- numeric(nsteps)
    hazprob[1] <- pnorm(1 * dstep, mean = mu, sd = std)
    
    for (s in 2:nsteps) {
      upper <- pnorm(s * dstep, mean = mu, sd = std)
      lower <- pnorm((s - 1) * dstep, mean = mu, sd = std)
      denom <- 1 - lower
      prob <- upper - lower
      hazprob[s] <- if (denom < 1e-6) 1.0 else prob / denom
    }
    
    return(hazprob)
  }
  
  # Hazard functions for Leave and Moult
  xzeta = hazFunction(mu.leave, sd.leave, nIntervals_s, daysPerBin)  # Leave hazard function
  xgamma = hazFunction(mu.moult, sd.moult, nIntervals_s, daysPerBin)  # Moult hazard function
  # plot(xzeta,type='l')
  # lines(xgamma,col='red')
  xphi <- plogis(survInt + survSlope * (0:(ndays)))  # Survival 
  xchi <- plogis(xchi)  # Carcass survival
  # plot(xphi,type='l')
  
  
  
  if (add_dead) {
    matarray <- array(0, dim = c(nIntervals_s, 5, 5)) 
    for (i in 1:nIntervals_s) {
      matarray[i, 1, 1] <- (xphi[(i-1)*daysPerBin+1]^daysPerBin * (1 - xgamma[i]) * (1 - xzeta[i]))  # W -> W #survive up to the current interval day
      matarray[i, 2, 1] <- (xphi[(i-1)*daysPerBin+1]^daysPerBin * xgamma[i] * (1 - xzeta[i]))        # W -> M
      matarray[i, 3, 1] <- (xphi[(i-1)*daysPerBin+1]^daysPerBin * xzeta[i])                             # W -> LA
      matarray[i, 4, 1] <- (1 - xphi[(i-1)*daysPerBin+1]^daysPerBin)                                    # W -> D
      
      matarray[i, 2, 2] <- (xphi[(i-1)*daysPerBin+1]^daysPerBin * (1 - xzeta[i]))  # M -> M
      matarray[i, 3, 2] <- (xphi[(i-1)*daysPerBin+1]^daysPerBin * xzeta[i])        # M -> LA
      matarray[i, 4, 2] <- (1 - xphi[(i-1)*daysPerBin+1]^daysPerBin)               # M -> D
      
      matarray[i, 4, 4] <- xchi^daysPerBin  # D -> D #survive as carcass for all days in interval
      matarray[i, 5, 4] <- (1 - xchi^daysPerBin)  # D -> LD
      
      matarray[i, 3, 3] <- 1  # LA is absorbing
      matarray[i, 5, 5] <- 1  # LD absorbing
    }
  } else {
    # If add_dead is F,  3 states
    matarray <- array(0, dim = c(nIntervals_s, 3, 3))  
    for (i in 1:nIntervals_s) {
      matarray[i, 1, 1] <- (1 - xgamma[i]) * (1 - xzeta[i])  # W -> W
      matarray[i, 2, 1] <- (xgamma[i] * (1 - xzeta[i]))      # W -> M
      matarray[i, 3, 1] <- xzeta[i]                          # W -> LA
      
      matarray[i, 2, 2] <- (1 - xzeta[i])                    # M -> M
      matarray[i, 3, 2] <- xzeta[i]                          # M -> LA
      
      matarray[i, 3, 3] <- 1  # LA is absorbing
    }
  }
  
  
  #------------------------------
  # Place cohorts over time
  #------------------------------
  stageN <- array(0, dim = c(nIntervals_s, ifelse(add_dead, 5, 3)))  
  for (c in 1:nIntervals_c) {
    cohort_vec <- rep(0, ifelse(add_dead, 5, 3))
    cohort_vec[1] <- cohort.size[c]  # all newly born cohorts start in stage 1
    
    for (s in c:nIntervals_s) {
      age <- s - c + 1  # age since birth (s-c) (1-based) so +1 is next age
      stageN[s, ] <- stageN[s, ] + cohort_vec           # track the cohort progression through stages at time s
      cohort_vec <- matarray[age, , ] %*% cohort_vec  # update the cohort's stage for the next interval
      
    }
  }
  
  
  #------------------------------
  # Observation process
  #------------------------------
  if (obs.type == "poisson") {
    obs_trans_mat <- array(0, dim = c(ifelse(add_dead, 5, 3), ifelse(add_dead, 4, 3)))  # 5 states for add_dead, 3 states for no dead
    
    # transition probabilities (adjust for 3 or 2 states as needed)
    if (use_pcw) {
      obs_trans_mat[1,1] <- pow * pcw                     # white -> white obs
      obs_trans_mat[1,2] <- pow * (1 - pcw)               # white -> moulting obs (misclassified)
    } else {
      obs_trans_mat[1,1] <- pow                           # white -> white obs
      obs_trans_mat[1,2] <- 0                             # white -> moulting obs (misclassified)
    }
    
    if (add_dead) {
      obs_trans_mat[1,4] <- 1 - pow                       # white -> dead obs
      obs_trans_mat[1,3] <- 0                             # white -> unobserved
      
      obs_trans_mat[2,1] <- pom * (1 - pcm)               # moulting -> white obs (misclassified)
      obs_trans_mat[2,2] <- pom * pcm                     # moulting -> moulting obs
      obs_trans_mat[2,3] <- 0                             # moulting -> dead obs
      obs_trans_mat[2,4] <- 1 - pom                       # moulting -> unobserved
      
      obs_trans_mat[3,1] <- 0                             # leave -> white obs
      obs_trans_mat[3,2] <- 0                             # leave -> moulting obs
      obs_trans_mat[3,3] <- 0                             # leave -> dead obs
      obs_trans_mat[3,4] <- 1                             # leave -> unobserved
      
      obs_trans_mat[4,1] <- 0                             # dead (obs) -> white obs
      obs_trans_mat[4,2] <- 0                             # dead (obs) -> moulting obs
      obs_trans_mat[4,3] <- pod                           # dead (obs) -> dead obs
      obs_trans_mat[4,4] <- 1 - pod                       # dead (obs) -> unobserved
      
      obs_trans_mat[5,1] <- 0                             # dead (unobs) -> white obs
      obs_trans_mat[5,2] <- 0                             # dead (unobs) -> moulting obs
      obs_trans_mat[5,3] <- 0                             # dead (unobs) -> dead obs
      obs_trans_mat[5,4] <- 1                             # dead (unobs) -> unobserved
    } else {
      # for 3 states (white, moulting, leave), just adjust to leave out dead states
      
      
      obs_trans_mat[1, 3] <- 1 - pow
      obs_trans_mat[2, 1] <- pom * (1 - pcm)
      obs_trans_mat[2, 2] <- pom * pcm
      obs_trans_mat[2, 3] <- 1 - pom
      obs_trans_mat[3, 3] <- 1
      
    }
    
    if (add_dead) {
      lambda_obs_vector <- array(0, dim = c(length(obsdays_bin), 4))  # 4 stages (white, moulted, left, dead(unobs))
      obs_vector <- matrix(0, nrow = length(obsdays_bin), ncol = 3)  # Obs classes: W , M, D
      
    } else {
      lambda_obs_vector <- array(0, dim = c(length(obsdays_bin), 3))  # 3 stages white,moulted,left 
      obs_vector <- matrix(0, nrow = length(obsdays_bin), ncol = 2)  # Obs classes: W & M
    }
    
    # Loop through each observation day
    for (obs in 1:length(obsdays)) {
      # obs=1
      if (add_dead) {
        #  Obs classes: W , M, D
        lambda_obs_vector[obs, ] <- stageN[obsdays_bin[obs], ] %*% obs_trans_mat
        obs_vector[obs, 1] <- rpois(1, lambda_obs_vector[obs, 1])  # observed white (W)
        obs_vector[obs, 2] <- rpois(1, lambda_obs_vector[obs, 2])  # observed moulted (M)
        obs_vector[obs, 3] <- rpois(1, lambda_obs_vector[obs, 3])  # observed dead (D)
      } else {
        # Obs classes: W & M
        lambda_obs_vector[obs, ] <- stageN[obsdays_bin[obs], ] %*% obs_trans_mat
        obs_vector[obs, 1] <- rpois(1, lambda_obs_vector[obs, 1])  # observed White (W)
        obs_vector[obs, 2] <- rpois(1, lambda_obs_vector[obs, 2])  # observed Moulted (M)
      }
    }
    
    
  } else if (obs.type == "mvn") {
    
    Ntot <- numeric(length(obsdays))
    Pw <- numeric(length(obsdays))
    Pm <- numeric(length(obsdays))
    
    Psi_w <- numeric(length(obsdays))
    Psi_m <- numeric(length(obsdays))
    
    if (add_dead) {
      Pd <- numeric(length(obsdays))
      Psi_d <- numeric(length(obsdays))
      
      var_white <- numeric(length(obsdays))
      var_moult <- numeric(length(obsdays))
      var_dead <- numeric(length(obsdays))
      
      covar_wm <- numeric(length(obsdays))
      covar_dw <- numeric(length(obsdays))
      covar_dm <- numeric(length(obsdays))
      
      cov_matrix <- array(0, dim = c(length(obsdays), 3, 3))
      mean_vector <- matrix(0, nrow = length(obsdays), ncol = 3)
      obs_vector <- matrix(0, nrow = length(obsdays), ncol = 3)
      
      for (obs in 1:length(obsdays)) {
        Ntot[obs] <- stageN[obsdays_bin[obs], 1] + stageN[obsdays_bin[obs], 2] + stageN[obsdays_bin[obs], 4]
        
        Pw[obs] <- stageN[obsdays_bin[obs], 1] / Ntot[obs]
        Pm[obs] <- stageN[obsdays_bin[obs], 2] / Ntot[obs]
        Pd[obs] <- stageN[obsdays_bin[obs], 4] / Ntot[obs]
        
        if (use_pcw) {
          Psi_w[obs] <- Pw[obs] * pow * pcw + Pm[obs] * pom * (1 - pcm)
          Psi_m[obs] <- Pm[obs] * pom * pcm + Pw[obs] * pow * (1 - pcw)
        } else {
          Psi_w[obs] <- Pw[obs] * pow + Pm[obs] * pom * (1 - pcm)
          Psi_m[obs] <- Pm[obs] * pom * pcm
        }
        Psi_d[obs] <- Pd[obs] * pod
        
        var_white[obs] <- Ntot[obs] * Psi_w[obs] * (1 - Psi_w[obs])
        var_moult[obs] <- Ntot[obs] * Psi_m[obs] * (1 - Psi_m[obs])
        var_dead[obs]  <- Ntot[obs] * Psi_d[obs] * (1 - Psi_d[obs])
        
        covar_wm[obs] <- -Ntot[obs] * Psi_w[obs] * Psi_m[obs]
        covar_dw[obs] <- -Ntot[obs] * Psi_w[obs] * Psi_d[obs]
        covar_dm[obs] <- -Ntot[obs] * Psi_m[obs] * Psi_d[obs]
        
        cov_matrix[obs, , ] <- matrix(c(
          var_white[obs], covar_wm[obs], covar_dw[obs],
          covar_wm[obs], var_moult[obs], covar_dm[obs],
          covar_dw[obs], covar_dm[obs], var_dead[obs]
        ), nrow = 3, byrow = TRUE)
        
        mean_vector[obs, ] <- c(
          Ntot[obs] * Psi_w[obs],
          Ntot[obs] * Psi_m[obs],
          Ntot[obs] * Psi_d[obs]
        )
        
        obs_vector[obs, ] <- MASS::mvrnorm(1, mean_vector[obs, ], cov_matrix[obs, , ])
      }
      
    } else {
      var_white <- numeric(length(obsdays))
      var_moult <- numeric(length(obsdays))
      covar_wm <- numeric(length(obsdays))
      
      cov_matrix <- array(0, dim = c(length(obsdays), 2, 2))
      mean_vector <- matrix(0, nrow = length(obsdays), ncol = 2)
      obs_vector <- matrix(0, nrow = length(obsdays), ncol = 2)
      
      for (obs in 1:length(obsdays)) {
        Ntot[obs] <- stageN[obsdays_bin[obs], 1] + stageN[obsdays_bin[obs], 2]
        
        Pw[obs] <- stageN[obsdays_bin[obs], 1] / Ntot[obs]
        Pm[obs] <- stageN[obsdays_bin[obs], 2] / Ntot[obs]
        
        if (use_pcw) {
          Psi_w[obs] <- Pw[obs] * pow * pcw + Pm[obs] * pom * (1 - pcm)
          Psi_m[obs] <- Pm[obs] * pom * pcm + Pw[obs] * pow * (1 - pcw)
        } else {
          Psi_w[obs] <- Pw[obs] * pow + Pm[obs] * pom * (1 - pcm)
          Psi_m[obs] <- Pm[obs] * pom * pcm
        }
        
        var_white[obs] <- Ntot[obs] * Psi_w[obs] * (1 - Psi_w[obs])
        var_moult[obs] <- Ntot[obs] * Psi_m[obs] * (1 - Psi_m[obs])
        covar_wm[obs]  <- -Ntot[obs] * Psi_w[obs] * Psi_m[obs]
        
        cov_matrix[obs, , ] <- matrix(c(
          var_white[obs], covar_wm[obs],
          covar_wm[obs], var_moult[obs]
        ), nrow = 2, byrow = TRUE)
        
        mean_vector[obs, ] <- c(
          Ntot[obs] * Psi_w[obs],
          Ntot[obs] * Psi_m[obs]
        )
        
        obs_vector[obs, ] <- MASS::mvrnorm(1, mean_vector[obs, ], cov_matrix[obs, , ])
      }
    }
  }
  return(list(
    Pb = Pb,
    cohort.size = cohort.size,
    obs_vector = obs_vector,
    stageN=stageN,
    #lambda_obs_vector = lambda_obs_vector,
    obsdays = obsdays
  ))
  
}


