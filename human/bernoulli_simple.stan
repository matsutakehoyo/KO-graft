data {
  int<lower=1> n_obs;
  int<lower=1> n_sex;
  int<lower=1> n_grp;
  int<lower=1> n_stm;
  int<lower=1> n_smp;
  int<lower=1, upper=n_sex> obs2sex[n_obs];
  int<lower=1, upper=n_grp> obs2grp[n_obs];
  int<lower=1, upper=n_stm> obs2stm[n_obs];
  int<lower=1, upper=n_smp> obs2smp[n_obs]; 
  real obs2spt[n_obs];
  int<lower=0, upper=1> response[n_obs];
  real<lower=0, upper=1> p_non_zero; 
}


parameters {
  // overall mean
  real a0; 
  // main effects
  vector[n_sex] a_sex;
  vector[n_grp] a_grp;
  vector[n_stm] a_stm;
  vector[n_smp] a_smp;
  real <lower=0> smp_s; 
  // covariates
  real spt;
  //guessing parameter for Binomial robustness
  real<lower=0,upper=1> guess;    
}

transformed parameters {

}

model {
  real theta[n_obs];
  // for robust estimate
  guess ~ beta(1,9);
  // overall mean
  a0 ~ normal(logit(p_non_zero), 1); 
  // main effects
  // a_grp ~ normal(0,1);
  a_grp ~ student_t(3,0,1);
  a_stm ~ student_t(3,0,1);
  a_smp ~ normal(0,smp_s);
  smp_s ~ student_t(3,0,1);

  // covariates
  spt ~ student_t(3,0,1);


  for ( i in 1:n_obs ) {
    theta[i] =.5*guess + (1-guess)*inv_logit(
      a0 + 
      a_grp[obs2grp[i]] +  
      a_sex[obs2sex[i]] + 
      a_smp[obs2smp[i]] + 
      a_stm[obs2stm[i]] + 
      spt * obs2spt[i] 
    );
  }
  // likelihood: Bernoulli
  for (i in 1:n_obs) {
    response[i] ~ bernoulli(theta[i]);
  }
}

generated quantities {
  real t0;
  vector[n_grp] t_grp;
  vector[n_sex] t_sex; 
  vector[n_smp] t_smp; 
  vector[n_stm] t_stm; 
  real m[n_sex, n_grp, n_stm, n_smp]; // cell means for b/l
  // estimated mean 
  real y_pred[n_sex, n_grp, n_stm];

  // convert unconstrained predicotrs to constrained (sum-to-zero)
  // --------------------------------------------------------------------------------
  // overall mean
  t0 = 0.0; 
  // cell means
  for (x in 1:n_sex){
    for(g in 1:n_grp){
      for(s in 1:n_stm){
        y_pred[x,g,s] = a0 + a_sex[x] + a_grp[g] + a_stm[s];
        for(a in 1:n_smp){
          m[x,g,s,a] = y_pred[x,g,s] + a_smp[a];
          t0 += m[x,g,s,a];
        }
      }
    }
  }

  t0 = t0 / num_elements(m);  
  // main effects
  for(x in 1:n_sex) t_sex[x] = mean(to_array_1d(m[x,:,:,:])) - t0;
  for(g in 1:n_grp) t_grp[g] = mean(to_array_1d(m[:,g,:,:])) - t0;
  for(s in 1:n_stm) t_stm[s] = mean(to_array_1d(m[:,:,s,:])) - t0; 
  for(a in 1:n_smp) t_smp[a] = mean(to_array_1d(m[:,:,:,a])) - t0; 
  
}

