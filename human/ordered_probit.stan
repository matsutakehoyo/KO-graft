data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1,upper=K> y[N];
  int<lower=1> n_grp; 
  int<lower=1, upper=n_grp> grp[N]; 
  int<lower=1> n_sex;
  int<lower=1, upper=n_sex> sex[N];
  int<lower=1> n_aml; 
  int<lower=1, upper=n_aml> aml[N]; 
  real c_lower;  // 1 + 0.5
  real c_upper;  // K - 0.5
}

parameters {
  // a is for mean and b for sd
  real a0; 
  real a_grp[n_grp];
  real a_sex[n_sex];
  real a_aml[n_aml];
  real<lower=0> a_aml_s; 
  real<lower=0> sigma;
  simplex[K-2] c_prop;

}

transformed parameters {
  ordered[K-1] c;  // cutoff: ordered vector with fixed upper and lower bounds
  c[1] = c_lower;

  for (i in 2:(K-1)) {
    c[i] = c[i-1] + c_prop[i-1] * (c_upper - c_lower);
  }

}

model {
  vector[K] theta[N];  // simplex[K] theta[N];  // probabilities
  real mu[N]; 

  for (i in 1:N){
    mu[i] = a0 + a_grp[grp[i]] + a_sex[sex[i]] + a_aml[aml[i]];
  }

  for (i in 1:N){
    theta[i,1] = Phi((c[1] - mu[i]) / sigma);
    theta[i,K] = 1 - Phi((c[K-1] - mu[i]) / sigma);

    for (k in 2:(K-1)) {
      theta[i,k] = Phi((c[k] - mu[i]) / sigma) - Phi((c[k-1] - mu[i]) / sigma);
    }           
  }

  a0 ~ normal((1 + K) / 2.0, K);
  a_grp ~ normal(0,1); 
  a_sex ~ normal(0,1); 
  a_aml ~ normal(0,a_aml_s); 
  a_aml_s ~ normal(0,1); 
  sigma ~ normal(0,1); 
  for (n in 1:N){
    y[n] ~ categorical(theta[n]); 
  }
}

generated quantities{
  // mu 
  real m0=0.0;
  vector[n_grp] m_grp;
  vector[n_sex] m_sex;
  vector[n_aml] m_aml;
  real m[n_grp, n_sex, n_aml]; 
  real theta_pred[K, n_grp, n_sex, n_aml]; 
    
  // convert unconstrained predicotrs to constrained (sum-to-zero)
  // cell means
  for(g in 1:n_grp){
    for(x in 1:n_sex){
      for(a in 1:n_aml){
        m[g,x,a] = a0 + a_grp[g] + a_sex[x] + a_aml[a];
        m0 += m[g,x,a];
        theta_pred[1,g,x,a] = Phi((c[1] - m[g,x,a]) / sigma);
        theta_pred[K,g,x,a] = 1 - Phi((c[K-1] - m[g,x,a]) / sigma);
          for (i in 2:(K-1)) {
            theta_pred[i,g,x,a] = Phi((c[i] - m[g,x,a]) / sigma) - Phi((c[i-1] - m[g,x,a]) / sigma);
          } 
      }
    }
  }
  m0 = m0 / num_elements(m);  
  // main effects
  for(g in 1:n_grp) m_grp[g] = mean(to_array_1d(m[g,:,:])) - m0;
  for(x in 1:n_sex) m_sex[x] = mean(to_array_1d(m[:,x,:])) - m0;
  for(a in 1:n_aml) m_aml[a] = mean(to_array_1d(m[:,:,a])) - m0;
}

