// greysealpups_Deads.stan
// Full model with explicit mortality and carcass persistence.


functions {
    
    vector skewnormBDFunction(real mu, real sigma, real alpha, int nsteps, vector days){
    vector[nsteps] xvec;
    for (i in 1:nsteps) {
      xvec[i] = exp(skew_normal_lpdf(days[i] | mu, sigma, alpha));
    }
    return xvec / sum(xvec);
  }

  vector hazFunction(real mu, real std, int nsteps, int dstep){
    real prob;
    vector[nsteps] hazprob;
    hazprob[1] = normal_cdf((1*dstep)|mu,std);
    for (s in 2:nsteps){
      prob = normal_cdf((s*dstep)|mu,std) - normal_cdf(((s-1)*dstep)|mu,std);
      hazprob[s] = (1-normal_cdf(((s-1)*dstep)|mu,std))<0.000001 ? 1.0 : prob/(1-normal_cdf(((s-1)*dstep)|mu,std));
    }
    return(hazprob);
  }
  
  matrix transMatrix(real xphi, real xgamma, real xzeta, real xchi, int stages, int dstep){
    matrix[stages,stages] matarray = rep_matrix(0,stages,stages); 
    matarray[1,1] = xphi^dstep * (1-xgamma) * (1-xzeta);
    matarray[2,1] = xphi^dstep * xgamma * (1-xzeta);
    matarray[3,1] = xphi^dstep * xzeta;
    matarray[4,1] = (1-xphi^dstep);
    matarray[2,2] = xphi^dstep * (1-xzeta);
    matarray[3,2] = xphi^dstep * xzeta;
    matarray[4,2] = (1-xphi^dstep);
    matarray[4,4] = xchi^dstep;
    matarray[5,4] = 1-xchi^dstep;
    matarray[3,3] = 1;
    matarray[5,5] = 1;
    return(matarray);
  }
  
} // end of functions block


data {
  int <lower=0> nstates;
  int <lower=0> maxDayBirth;
  int <lower=0> ndays;
  int <lower=0> nobsstates;
  int <lower=0> nsurveys;
  array[nsurveys, nobsstates] int obs_vector;
  array[nsurveys] int <lower=0> obsdays;
  
  // Detection and classification rates (passed from R)
  real<lower=0,upper=1> pobs_w;    // detection probability: whitecoats
  real<lower=0,upper=1> pobs_m;    // detection probability: moulted
  real<lower=0,upper=1> pobs_d;    // detection probability: dead
  real<lower=0,upper=1> pcm;       // probability of correctly classifying moulted
  real<lower=0,upper=1> pcw;       // probability of correctly classifying whitecoat
}
 
transformed data {
  int <lower=0> daysPerBin = 4;
  int nIntervals_s = ndays / daysPerBin;
  
  array[nsurveys] int <lower=0> obsdays_bin; 
  for (z in 1:nsurveys) obsdays_bin[z] = obsdays[z] / daysPerBin;
  
  // Cohort bins
  int <lower=0> nIntervals_c = 25;
  vector[nIntervals_c] cohort_days;
  for (z in 1:nIntervals_c) {
    cohort_days[z] = z * daysPerBin - 2.0;
  } 
  
  // Survival
  real <lower=0> survInt = 4;
  real <lower=0> survSlope = 0.1;  
  vector <lower=0,upper=1>[ndays] xphi;
  for (a in 1:ndays) {
    xphi[a] = inv_logit(survInt + survSlope * (a-1));
  }
  
  // Carcass persistence
  real <lower=0> pcarcSurvInt = 4;
  vector <lower=0,upper=1>[ndays] xchi;
  for (a in 1:ndays) {
    xchi[a] = inv_logit(pcarcSurvInt);
  }
         
  // Moult process
  real <lower=0> muMoult = 23;  
  real <lower=0> sdMoult = 5; 
  vector <lower=0.0,upper=1>[ndays %/% daysPerBin] xgamma = hazFunction(muMoult, sdMoult, nIntervals_s, daysPerBin);

  // Leave process
  real <lower=0> muLeave = 31.5;
  real <lower=0> sdLeave = 7;
  
} // end of transformed data block

parameters {
  real nborn_raw;
  real muBday_raw;
  real<lower=0> sdBday_raw;
  real<lower=0> alphaBday_raw;
}

transformed parameters {
  real nborn = exp(5.0 + 2.0 * nborn_raw);
  real muBday = exp(3.8 + 0.15 * muBday_raw);
  real sdBday = sdBday_raw * 6;
  real alphaBday = alphaBday_raw * 4;

  // Birth process
  vector[nIntervals_c] Pb;
  Pb = skewnormBDFunction(muBday, sdBday, alphaBday, nIntervals_c, cohort_days);  

  // Leave process
  vector[nIntervals_s] xzeta;
  xzeta = hazFunction(muLeave, sdLeave, nIntervals_s, daysPerBin);

  // Distribute births across cohorts
  vector[nIntervals_c] cohortSizeBin; 
  cohortSizeBin = nborn * Pb;

  // Transition matrices (one per age interval)
  array[nIntervals_s] matrix[nstates, nstates] matarray;
  for (i in 1:nIntervals_s) {
    matarray[i] = transMatrix(xphi[((i-1)*daysPerBin)+1], xgamma[i], xzeta[i], xchi[((i-1)*daysPerBin)+1], nstates, daysPerBin);
  }

  // Accumulate stage counts across cohorts and time
  vector[nstates] tempCohortVec;
  array[nIntervals_s] vector[nstates] stageN;
  for (z in 1:nIntervals_s) stageN[z] = rep_vector(0, nstates);

  for (c in 1:nIntervals_c) {
    tempCohortVec = rep_vector(0, nstates);
    tempCohortVec[1] = cohortSizeBin[c];
    stageN[c] += tempCohortVec;
    for (a in 1:(nIntervals_s-c)) {
      tempCohortVec = matarray[a] * tempCohortVec;
      stageN[c+a] += tempCohortVec;
    }
  }

  // Observation matrix
  matrix[nstates, nobsstates+1] obsMatrix = rep_matrix(0, nstates, nobsstates+1); 
  obsMatrix[1,1] = pobs_w * pcw;
  obsMatrix[1,2] = pobs_w * (1 - pcw);
  obsMatrix[1,4] = 1 - pobs_w;
  obsMatrix[2,1] = pobs_m * (1 - pcm);
  obsMatrix[2,2] = pobs_m * pcm;
  obsMatrix[2,4] = 1 - pobs_m;
  obsMatrix[3,4] = 1;
  obsMatrix[4,3] = pobs_d;
  obsMatrix[4,4] = 1 - pobs_d;
  obsMatrix[5,4] = 1;

  array[nsurveys] row_vector[nobsstates + 1] lambda_obs_vector;
  for (obs in 1:nsurveys) {
    lambda_obs_vector[obs] = stageN[obsdays_bin[obs]]' * obsMatrix;
  }
}

model {
  // Priors
  nborn_raw ~ std_normal();  
  muBday_raw ~ normal(0, 1);
  sdBday_raw ~ gamma(2, 1);
  alphaBday_raw ~ normal(0.8, 0.5); 

  // Likelihood
  for (obs in 1:nsurveys) {
    obs_vector[obs, 1] ~ poisson(lambda_obs_vector[obs, 1]);
    obs_vector[obs, 2] ~ poisson(lambda_obs_vector[obs, 2]);
    obs_vector[obs, 3] ~ poisson(lambda_obs_vector[obs, 3]);
  }
}

generated quantities {
  real pups_leave_before_moult = 0;

  {
    vector[nstates] tempCohortCalcVec;
    vector[nstates] nextCohortCalcVec;

    for (c in 1:nIntervals_c) {
      tempCohortCalcVec = rep_vector(0, nstates);
      tempCohortCalcVec[1] = cohortSizeBin[c];
      for (a in 1:(nIntervals_s - c)) {
        matrix[nstates, nstates] tmat = matarray[a];
        real from_W = tempCohortCalcVec[1];
        real to_LA_from_W = from_W * tmat[3, 1];
        if (tempCohortCalcVec[2] < 1e-3) {
          pups_leave_before_moult += to_LA_from_W;
        }
        nextCohortCalcVec = tmat * tempCohortCalcVec;
        tempCohortCalcVec = nextCohortCalcVec;
      }
    }
  }
  
  real prop_leave_before_moult = pups_leave_before_moult / nborn;
  
  array[nsurveys, nobsstates] int yrep;
  for (obs in 1:nsurveys) {
    for (j in 1:nobsstates) {
      yrep[obs, j] = poisson_rng(lambda_obs_vector[obs, j]);
    }
  }
}

