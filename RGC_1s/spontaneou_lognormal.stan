data{
  int<lower=1> n_obs;
  int<lower=1> n_res; 
  int<lower=1> n_lin;
  int<lower=1> n_cnd;
  int<lower=1> n_stm;
  int<lower=1> n_tpl;
  int<lower=1> n_tim;
  int<lower=1> n_l7g; 
  int<lower=1> n_dmg;
  int<lower=1> n_mus;
  int<lower=1, upper=n_lin> obs2lin[n_obs];
  int<lower=1, upper=n_res> obs2res[n_obs];
  int<lower=1, upper=n_cnd> obs2cnd[n_obs];
  int<lower=1, upper=n_stm> obs2stm[n_obs];
  int<lower=1, upper=n_tpl> obs2tpl[n_obs];
  int<lower=1, upper=n_tim> obs2tim[n_obs];
  int<lower=1, upper=n_l7g> obs2l7g[n_obs];   
  int<lower=1, upper=n_dmg> obs2dmg[n_obs];
  int<lower=1, upper=n_mus> obs2mus[n_obs]; 
  real logmean; 
  real logsd; 
  vector<lower=0>[n_obs] y;
}

parameters{
  real a0; 
  vector[n_lin] a_lin;
  vector[n_res] a_res;
  vector[n_cnd] a_cnd;
  vector[n_stm] a_stm;
  vector[n_tim] a_tim; 
  vector[n_l7g] a_l7g;
  vector[n_dmg] a_dmg;
  vector[n_mus] a_mus;
  real <lower=0> m_mus_s; 
  
  // interactions: 
  matrix[n_lin, n_res] a_lin_res;
  matrix[n_lin, n_tim] a_lin_tim;
  matrix[n_lin, n_l7g] a_lin_l7g; 
  real<lower=0> sdlog;
}

model{
  vector[n_obs] mu; 

  a0 ~ normal(logmean, .5); 
  a_lin ~ normal(0,1);
  a_res ~ normal(0,1);
  a_cnd ~ normal(0,1);
  a_stm ~ normal(0,1);
  a_tim ~ normal(0,1); 
  a_l7g ~ normal(0,1); 
  a_dmg ~ normal(0,1); 
  a_mus ~ normal(0,m_mus_s);
  m_mus_s ~ normal(0,1);  

  // interactions
  for(l in 1:3){  //index line
      a_lin_tim[l] ~ normal(0,1); 
      a_lin_l7g[l] ~ normal(0,1); 
      a_lin_res[l] ~ normal(0,1); 
  }


  for (i in 1:n_obs) {
    mu[i] = a0 +
      a_lin[obs2lin[i]] + a_res[obs2res[i]] + a_cnd[obs2cnd[i]] +
      a_stm[obs2stm[i]] + a_tim[obs2tim[i]] + a_l7g[obs2l7g[i]] + 
      a_dmg[obs2dmg[i]] + a_mus[obs2mus[i]] + 
      a_lin_tim[obs2lin[i], obs2tim[i]] + 
      a_lin_l7g[obs2lin[i], obs2l7g[i]] + 
      a_lin_res[obs2lin[i], obs2res[i]];
    y[i] ~ lognormal(mu[i], sdlog);   
  }
  
}

generated quantities {
  // mu 
  real m0;
  vector[n_lin] m_lin;
  vector[n_res] m_res;
  vector[n_cnd] m_cnd; 
  vector[n_stm] m_stm; 
  vector[n_tim] m_tim;
  vector[n_l7g] m_l7g;  
  vector[n_dmg] m_dmg; 
  vector[n_mus] m_mus; 
  matrix[n_lin, n_tim] m_lin_tim;
  matrix[n_lin, n_l7g] m_lin_l7g;  
  matrix[n_lin, n_res] m_lin_res;  
  real m[n_lin, n_res, n_cnd, n_stm, n_tim, n_l7g, n_dmg, n_mus]; 
  real y_pred_m[n_lin, n_res, n_cnd, n_stm, n_tim, n_l7g];


  // convert unconstrained predicotrs to constrained (sum-to-zero)
  // --------------------------------------------------------------------------------
  // overall mean
  m0 = 0.0; 
  // cell means
  for(l in 1:n_lin){
    for(r in 1:n_res){
      for(c in 1:n_cnd){
        for(i in 1:n_stm){
          for(a in 1:n_tim){
            for(g in 1:n_l7g){
              y_pred_m[l,r,c,i,a,g] = a0 +
                a_lin[l] + a_res[r] + a_cnd[c] + a_stm[i] + a_tim[a] + a_l7g[g] +
                a_lin_tim[l,a] + a_lin_l7g[l,g] + a_lin_res[l,r];
              for (d in 1:n_dmg){
                for (u in 1:n_mus){
                  m[l,r,c,i,a,g,d,u] = y_pred_m[l,r,c,i,a,g] +
                    a_dmg[d] + a_mus[u]; 
                  m0 += m[l,r,c,i,a,g,d,u];
                }
              }
            }
          }
        }     
      }
    }
  }


  m0 = m0 / num_elements(m);  
  // main effects
  for(l in 1:n_lin) m_lin[l] = mean(to_array_1d(m[l,:,:,:,:,:,:,:])) - m0;
  for(r in 1:n_res) m_res[r] = mean(to_array_1d(m[:,r,:,:,:,:,:,:])) - m0;
  for(c in 1:n_cnd) m_cnd[c] = mean(to_array_1d(m[:,:,c,:,:,:,:,:])) - m0;
  for(i in 1:n_stm) m_stm[i] = mean(to_array_1d(m[:,:,:,i,:,:,:,:])) - m0;
  for(a in 1:n_tim) m_tim[a] = mean(to_array_1d(m[:,:,:,:,a,:,:,:])) - m0;
  for(g in 1:n_l7g) m_l7g[g] = mean(to_array_1d(m[:,:,:,:,:,g,:,:])) - m0;
  for(d in 1:n_dmg) m_dmg[d] = mean(to_array_1d(m[:,:,:,:,:,:,d,:])) - m0;
  for(u in 1:n_mus) m_mus[u] = mean(to_array_1d(m[:,:,:,:,:,:,:,u])) - m0;

  // interactions: 
  for (l in 1:n_lin){
    for (r in 1:n_res)
      m_lin_res[l,r] = mean(to_array_1d(m[l,r,:,:,:,:,:,:])) - (m_lin[l] + m_res[r] + m0);
    for (a in 1:n_tim)
      m_lin_tim[l,a] = mean(to_array_1d(m[l,:,:,:,a,:,:,:])) - (m_lin[l] + m_tim[a] + m0);
    for (g in 1:n_l7g)
      m_lin_l7g[l,g] = mean(to_array_1d(m[l,:,:,:,:,g,:,:])) - (m_lin[l] + m_l7g[g] + m0);
  }
}



