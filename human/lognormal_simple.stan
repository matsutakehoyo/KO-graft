data{
  int<lower=1> n_obs;
  int<lower=1> n_sex; 
  int<lower=1> n_grp;
  int<lower=1> n_stm;
  int<lower=1> n_cnd;
  int<lower=1> n_smp;
  int<lower=1, upper=n_sex> obs2sex[n_obs];
  int<lower=1, upper=n_grp> obs2grp[n_obs];
  int<lower=1, upper=n_stm> obs2stm[n_obs];
  int<lower=1, upper=n_cnd> obs2cnd[n_obs];
  int<lower=1, upper=n_smp> obs2smp[n_obs]; 
  real logmean; 
  real logsd; 
  vector<lower=0>[n_obs] y;
}

parameters{
  real a0; 
  vector[n_sex] a_sex;
  vector[n_grp] a_grp;
  vector[n_stm] a_stm;
  vector[n_cnd] a_cnd;
  vector[n_smp] a_smp;
  real <lower=0> m_smp_s; 
  real<lower=0> sdlog;
}

model{
  vector[n_obs] mu; 
  vector[n_obs] sigma; 

  a0 ~ normal(logmean, 1); 
  a_grp ~ normal(0,1);
  a_sex ~ normal(0,1); 
  a_stm ~ normal(0,1);
  a_cnd ~ normal(0,1);
  a_smp ~ normal(0,m_smp_s);
  m_smp_s ~ normal(0,1);  

  for (i in 1:n_obs) {
    mu[i] = a0 +
      a_sex[obs2sex[i]] +
      a_grp[obs2grp[i]] + 
      a_stm[obs2stm[i]] +
      a_cnd[obs2cnd[i]] +
      a_smp[obs2smp[i]];

    y[i] ~ lognormal(mu[i], sdlog);   
  }
  
}

generated quantities {
  // mu 
  real m0;
  vector[n_sex] m_sex; 
  vector[n_grp] m_grp;
  vector[n_stm] m_stm; 
  vector[n_cnd] m_cnd; 
  vector[n_smp] m_smp; 
  real m[n_sex, n_grp, n_stm, n_cnd, n_smp]; // cell means for b/l
  real y_pred_m[n_sex, n_grp, n_stm, n_cnd];
  
  // overall mean
  m0 = 0.0; 
  // cell means

  for (x in 1:n_sex){
    for(g in 1:n_grp){
      for(s in 1:n_stm){
        for(c in 1:n_cnd){
          y_pred_m[x,g,s,c] = 
            a0 +
            a_sex[x] + 
            a_grp[g] + 
            a_stm[s] + 
            a_cnd[c];
          for(a in 1:n_smp){
            m[x,g,s,c,a] = 
              y_pred_m[x,g,s,c] + a_smp[a];
            m0 += m[x,g,s,c,a];
          }
        }
      }
    }
  }


  m0 = m0 / num_elements(m);  
  // main effects
  for(x in 1:n_sex) m_sex[x] = mean(to_array_1d(m[x,:,:,:,:])) - m0;
  for(g in 1:n_grp) m_grp[g] = mean(to_array_1d(m[:,g,:,:,:])) - m0;
  for(s in 1:n_stm) m_stm[s] = mean(to_array_1d(m[:,:,s,:,:])) - m0;
  for(c in 1:n_cnd) m_cnd[c] = mean(to_array_1d(m[:,:,:,c,:])) - m0; 
  for(a in 1:n_smp) m_smp[a] = mean(to_array_1d(m[:,:,:,:,a])) - m0;

}
