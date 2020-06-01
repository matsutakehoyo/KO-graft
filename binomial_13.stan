data {
  int<lower=1> n_obs; 
  int<lower=1> n_mus; 
  int<lower=1> n_mus_tp; 
  int<lower=1> n_grp; 
  int<lower=1> n_sex;
  int<lower=1, upper=n_mus> obs2mus[n_obs]; 
  int<lower=1, upper=n_grp> obs2grp[n_obs]; 
  int<lower=1, upper=n_sex> obs2sex[n_obs]; 
  int<lower=0, upper=1> tp[n_obs];
  int y_age[n_obs]; 
  int y_age_tp[n_obs];
  int y_iti[n_obs]; 
  int y_trl[n_obs];
  int<lower=0, upper=30> y[n_obs]; 
  real y_mean;
  real iti_mean;
  real trl_mean; 
  real age_tp_mean;  
  real age_mean;
}

parameters {
  real a0;
  real<lower=0> mus_s;
  vector[n_mus] a_mus;
  vector[n_grp] a_grp;
  vector[n_sex] a_sex;
  matrix[n_grp, n_sex] a_grp_sex;
  real iti_a; 
  real iti_b;
  real trl; 
  // real<lower=0, upper=1> guess; 
}

transformed parameters{

}

model {
  vector[n_obs] alpha; 

  a0 ~ normal(logit(y_mean), 1);
  mus_s ~ normal(0,1);
  a_mus ~ normal(0, mus_s);
  a_grp ~ normal(0, 1);
  a_sex ~ normal(0,1);
  for (g in 1:n_grp){
    a_grp_sex[g] ~ normal(0,1);
  }

  for (i in 1:n_obs){
    alpha[i] = inv_logit(//.5*guess + (1-guess)* 
      (a0 +
       a_grp[obs2grp[i]] +
       a_sex[obs2sex[i]] +
       a_mus[obs2mus[i]] +
       a_grp_sex[obs2grp[i], obs2sex[i]] +
       // tp[i] * age_tp * (y_age_tp[i] - age_tp_mean) +
       // age * (y_age[i] - age_mean) +
       trl * (y_trl[i] - trl_mean) +
       (iti_a*(y_iti[i])^2 + iti_b * y_iti[i] ) 
      ) 
    );        
  }


  y ~ binomial(30, alpha);  
}

generated quantities {
  // for loo 
  // vector[n_obs] log_lik;
  real b0; 
  vector[n_grp] b_grp; 
  vector[n_mus] b_mus; 
  vector[n_sex] b_sex;
  matrix[n_grp, n_sex] b_grp_sex;
  matrix[n_grp, n_sex] m_grp_sex;
  real y_pred[n_grp, n_sex]; // grp, sex
  real m[n_grp, n_sex, n_mus]; 
  // for loo
  // for (i in 1:n_obs){
  //   log_lik[i] = binomial_logit_lpmf(y[i] | 30, alpha[i]);
  // }

  // -------------------------------------------------------------------------------   
  // convert a0 to b0 mean of all deflections, and each effec to sum-to-zero
  b0 = 0.0;
  for (g in 1:n_grp){
    for (s in 1:2){
      y_pred[g,s] = a0 + a_grp[g] + a_sex[s] + a_grp_sex[g,s] ;   
        
      for (a in 1:n_mus){
        m[g,s,a] = a0 + a_grp[g] + a_sex[s] + a_grp_sex[g,s] + a_mus[a];   
        b0 += m[g,s,a];     
      }
    }
  }
  b0 = b0 / num_elements(m); 
  
  // main effects
  for (g in 1:n_grp) b_grp[g] = mean(to_array_1d(m[g,:,:])) - b0; 
  for (a in 1:n_mus) b_mus[a] = mean(to_array_1d(m[:,:,a])) - b0; 
  for (s in 1:2) b_sex[s] = mean(to_array_1d(m[:,s,:])) - b0;
  
  // interactions: line x cond
  for (g in 1:n_grp){
    for (s in 1:n_sex){
      m_grp_sex[g,s] = mean(to_array_1d(m[g,s,:]));
      b_grp_sex[g,s] = m_grp_sex[g,s] - (b_grp[g] + b_sex[s] + b0);
    }
  }
}

