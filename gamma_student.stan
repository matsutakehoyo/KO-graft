data{
  int<lower=1> n_obs;
  int<lower=1> n_grp;
  int<lower=1, upper=n_grp> obs2grp[n_obs];
  real obs2age[n_obs]; 
  int<lower=1> n_sex; 
  int<lower=1, upper=n_sex> obs2sex[n_obs]; 
  real y[n_obs];          
  real mean_y; 
  real mean_t; 
}

parameters{
  real m_age; 
  real mu0;  
  vector[n_grp] mu_grp; 
  vector[n_sex] mu_sex; 
  real<lower=0> va;
  real t_age; 
  real th0; 
  vector[n_grp] th_grp; 
  vector[n_sex] th_sex; 
}

model{
  real mu[n_obs]; 
  real th[n_obs]; 

  mu0 ~ normal(log(mean_y), 1); 
  mu_grp[1:4] ~ normal(0,1); 
  mu_sex ~ normal(0,1); 
  m_age ~ normal(0,1);   
  
  th0 ~ normal(logit(mean_t),1); 
  th_grp[1:4] ~ student_t(3,0,1); 
  th_sex ~ student_t(3,0,1); 
  t_age ~ student_t(3,0,1); 

  for (i in 1:n_obs) {
    mu[i] = exp(
      mu0 + 
      mu_grp[obs2grp[i]] + 
      mu_sex[obs2sex[i]] +
      m_age * obs2age[i]
    );
    
    th[i] = inv_logit(
      th0 + 
      th_grp[obs2grp[i]] + 
      th_sex[obs2sex[i]] + 
      t_age * obs2age[i]
    );

    if (y[i] == 0.0 ) {
      target += bernoulli_lpmf(0 | th[i]);
    } else {
      target += bernoulli_lpmf(1 | th[i]) +
                gamma_lpdf(y[i] | 1/va, 1/va/mu[i]);
    }
  }
}

generated quantities{
  real m0;
  real m_grp[n_grp]; 
  real m_sex[n_sex]; 
  real t0;
  real t_grp[n_grp]; 
  real t_sex[n_sex]; 
  real m[n_grp, n_sex]; 
  real t[n_grp, n_sex]; 

  m0=0.0; 
  t0=0.0; 
  for (g in 1:n_grp){
    for (s in 1:n_sex){
      m[g,s] = mu0 + mu_grp[g] + mu_sex[s];
      t[g,s] = th0 + th_grp[g] + th_sex[s];
      m0 += m[g,s]; 
      t0 += t[g,s]; 
    }
  }
  m0 = m0 / num_elements(m); 
  t0 = t0 / num_elements(t); 

  for(g in 1:n_grp) {
    m_grp[g] = mean(to_array_1d(m[g,:])) - m0; 
    t_grp[g] = mean(to_array_1d(t[g,:])) - t0; 
  }

  for(s in 1:n_sex) {
    m_sex[s] = mean(to_array_1d(m[:,s])) - m0; 
    t_sex[s] = mean(to_array_1d(t[:,s])) - t0; 
  }

}

