
functions{
    
    vector skewnormBDFunction(real mu, real sigma, real alpha, int nsteps, int dstep, vector days){
    vector[nsteps] xvec;
    for (i in 1:nsteps) {
      xvec[i] = exp(skew_normal_lpdf(days[i] | mu, sigma, alpha));
    }
    return xvec / sum(xvec);
  }


  vector hazFunction(real mu, real std, int nsteps, int dstep){
    real prob;         // temporary structure- probability mass of transitioning within each time interval
    vector[nsteps] hazprob;      // hazard for each bin (prob of transitioning within this time interval, conditional on being in the interval at the start of the time bin)
    hazprob[1] = normal_cdf((1*dstep)|mu,std);     // prob of remaining in state over first dstep days of life
    for (s in 2:nsteps){     // can we vectorize this?
      prob = normal_cdf((s*dstep)|mu,std) - normal_cdf(((s-1)*dstep)|mu,std);                // prob of remaining in state over this time bin, unconditional
      hazprob[s] = (1-normal_cdf(((s-1)*dstep)|mu,std))<0.000001 ? 1.0 : prob/(1-normal_cdf(((s-1)*dstep)|mu,std));         // prob of remaining in state over this time bin, conditional on being in the state at the start of the period
    }
    return(hazprob);
  }
  
  

  matrix transMatrix(real xgamma, real xzeta, int stages){
    matrix[stages, stages] matarray = rep_matrix(0, stages, stages);
    matarray[1, 1] =  (1 - xgamma) * (1 - xzeta);  // W to W
    matarray[2, 1] = xgamma * (1 - xzeta);        // W to M
    matarray[3, 1] =  xzeta;                         // W to LA
    matarray[2, 2] =  (1 - xzeta);                  // M to M
    matarray[3, 2] =  xzeta;                         // M to LA
    matarray[3, 3] = 1;                                         // LA is absorbing
    return matarray;
  }
}

data {
  int <lower=0> nstates;                            // number of pup states (3 states)
  int <lower=0> maxDayBirth;                        // last day of potential births
  int <lower=0> ndays;                              // number of days in the season
  int <lower=0> nobsstates;                         // number of observable states
  int <lower=0> nsurveys;                          // number of surveys
  array[nsurveys, nobsstates] int obs_vector;      // MAIN DATA: counts of pups in each class for each survey
  array[nsurveys] int <lower=0> obsdays;           // date (days since start of season) at which each survey is conducted
}

transformed data {
  int <lower=0> daysPerBin = 4; // DEFINE TIME BINS
  
  int nIntervals_s = ndays %/% daysPerBin;
  
  array[nsurveys] int <lower=0> obsdays_bin; 
  for (z in 1:nsurveys) obsdays_bin[z] = obsdays[z] / daysPerBin;

  // Define Cohort Bins
  int <lower=0> nIntervals_c = 25;
  vector[nIntervals_c] cohort_days;
  for (z in 1:nIntervals_c) {
    cohort_days[z] = z * daysPerBin - 2.0;
  }

  
  // MOULT
  real <lower=0> muMoult = 23;  
  real <lower=0> sdMoult = 5; 
  vector <lower=0.0,upper=1>  [ndays%/%daysPerBin] xgamma = hazFunction(muMoult,sdMoult,nIntervals_s,daysPerBin);
  
  // LEAVE
  real <lower=0> muLeave = 31.5;
  real <lower=0> sdLeave = 7;

  // Detection
  //real pobs_d = 0.95;
  real pobs_w = 0.95;
  real pobs_m = 0.95;
  real pcm = 0.91;
  real pcw=1;
}

parameters {
  real nborn_raw; // "superpopulation": total number born across season
  real muBday_raw; 
  real<lower=0>sdBday_raw;
  real<lower=0> alphaBday_raw;
  //real muLeave_raw; // mean leave date
  //real <lower=0.0> sdLeave_raw;    // standard deviation of leave date
  
   //obs params
//    real <lower=0.8,upper=1> pobs_d;      // prob of observing deads
//    real <lower=0.8,upper=1> pobs_w;      // prob of observing whites
//    real <lower=0.8,upper=1> pobs_m;      // prob of observing moulted
      //real <lower=0.7,upper=1> pcm;         // prob of misclassifying moulted
 //real <lower=0,upper=1> pcw

}

transformed parameters {
  real nborn = exp(5.0 + 2.0*nborn_raw);   // changed to log scale
  //real muBday = 50.0 + muBday_raw * 10.0;  // shift and scale
  //real muBday = exp(4.2 + 0.15*muBday_raw);   // changed to log scale
  real muBday = exp(3.8 + 0.15 * muBday_raw);
  //real sdBday = 1.0 + sdBday_raw * 4;    //0.5+sdBday_raw*5.0;
  real sdBday = sdBday_raw * 6;  // rescale
  real alphaBday = alphaBday_raw * 4;
  //real muLeave = 30.0 + muLeave_raw * 2.5; // shift and scale
  //real sdLeave = 0.5 + sdLeave_raw*5.0;  //shift and rescale (might need to force a little distance from zero)
}

model {
  // NBORN
  nborn_raw ~ std_normal();  

  // BIRTHDAY
  muBday_raw ~ normal(0, 1);
  //muBday_raw ~ std_normal(); 
  sdBday_raw ~ gamma(2, 1);  // set lower bound of prior to contrain to positive values 
  //sdBday_raw ~ gamma(0.1, 0.1); // seems to work okay- plenty of room to explore? 
  //alphaBday_raw ~ std_normal(); // skewness of birth distribution
  alphaBday_raw ~ normal(0.8, 0.5); 
  // LEAVING SUCCESSFULLY
 // muLeave_raw ~ student_t(3, 0, 2);//~ normal(0, 2)//~ std_normal(); 
  //sdLeave_raw ~ gamma(1,1);    // was uniform 1,10
 
    
  //   // DETECTION
  // pobs_d ~  beta(45.0,5.0);     // strong priors on observation process model
  // pobs_w ~ beta(45.0,5.0);
  // pobs_m ~ beta(45.0,5.0);
  // pcm ~ beta(35.0,5.0);
   //pcm ~ beta(91, 9);  
   
   
  // Birth process
  vector[nIntervals_c] Pb;
  Pb = skewnormBDFunction(muBday, sdBday, alphaBday, nIntervals_c, daysPerBin, cohort_days);

  // Leave process
  vector[nIntervals_s] xzeta;
  xzeta = hazFunction(muLeave, sdLeave, nIntervals_s, daysPerBin);

  // Distribute cohorts through time
  vector[nIntervals_c] cohortSizeBin;
  cohortSizeBin = nborn * Pb;  // Distribute nborn across cohorts
   #print("cohort size is ",cohortSizeBin);

 // transition matrix (one matrix defined per day of cohort age)
    array[nIntervals_s] matrix[nstates,nstates] matarray;

 for (i in 1:nIntervals_s) {
  matarray[i] = transMatrix(xgamma[i], xzeta[i], nstates);
}
  // array[nIntervals_s] matrix[nstates, nstates] matarray;
  // for (i in 1:nIntervals_s) {
  //   matarray[i] = transMatrix(xphi[((i - 1) * daysPerBin) + 1], xgamma[i], xzeta[i], xchi[((i - 1) * daysPerBin) + 1], nstates, daysPerBin);
  // }
   #print("matarray15 =",matarray[15]);
  //place cohorts over time
  vector[nstates] tempCohortVec;
  array[nIntervals_s] vector[nstates] stageN;
  for (z in 1:nIntervals_s) stageN[z] = rep_vector(0, nstates);

  for (c in 1:nIntervals_c) {  // loop through cohorts
    tempCohortVec = rep_vector(0, nstates);
    tempCohortVec[1] = cohortSizeBin[c];  // Initialize the cohort at age 0
    stageN[c] += tempCohortVec;           // add to the cumulative matrix of pups in each class
    for (a in 1:(nIntervals_s - c)) {      // loop through age of cohort
      tempCohortVec = matarray[a] * tempCohortVec;
      stageN[c + a] += tempCohortVec;    // add to the cumulative matrix of pups in each class
    }
  }
  
  // print("stageN interval 10 = ",stageN[10]);
  //=====================
  // Observation matrix
  //=====================

  matrix[nstates, nobsstates+1] obsMatrix = rep_matrix(0, nstates, nobsstates+1);
 //  white transitions
  obsMatrix[1,1] = pobs_w *pcw;                  // white -> observed white
  obsMatrix[1,2] = pobs_w * (1 - pcw);    // white -> observed moulted (misclassified)
  obsMatrix[1,3] = 1 - pobs_w;             // white -> unobserved

  //moulted transitions
  obsMatrix[2,1] = pobs_m * (1 - pcm);  // moulted -> observed white (misclassified)
  obsMatrix[2,2] = pobs_m * pcm;        // moulted -> observed moulted (correct classification)
  obsMatrix[2,3] = 1 - pobs_m;          // moulted -> unobserved

  obsMatrix[3,3] = 1;                   // left -> unobserved (stays in unobserved)


 array[nsurveys] row_vector[nobsstates+1] lambda_obs_vector;
for(obs in 1:nsurveys){
  lambda_obs_vector[obs] = stageN[obsdays_bin[obs]]' * obsMatrix ;
  
  // print("lambda_obs_vector_nsurvs =",lambda_obs_vector[obs]);
  // print("log prob before pois = ", target());
  
  obs_vector[obs,1] ~ poisson(lambda_obs_vector[obs,1] ) ; # data node- observed whites
  obs_vector[obs,2] ~ poisson(lambda_obs_vector[obs,2] ) ; # data node- observed moulteds
  
  // print("log prob after pois = ", target());

    }

 }


/*generated quantities {
 vector[nIntervals_c] Pb_out;
  Pb_out = skewnormBDFunction(muBday, sdBday, alphaBday, nIntervals_c, daysPerBin, cohort_days);
  }*/
  



