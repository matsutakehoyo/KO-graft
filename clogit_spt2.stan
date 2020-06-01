data{
  int<lower=1> n_obs;
  int<lower=1> n_cat; 
  int<lower=1> n_res; 
  int<lower=1> n_mus; 
  int<lower=1> n_lin;
  int<lower=1> n_cnd;
  int<lower=1> n_tpl;
  int<lower=1> n_tim;
  int<lower=1> n_stm;
  int<lower=1> n_dmg;
  int<lower=1> n_l7g;
  int<lower=1, upper=n_lin> obs2lin[n_obs];
  int<lower=1, upper=n_res> obs2res[n_obs];
  int<lower=1, upper=n_mus> obs2mus[n_obs];
  int<lower=1, upper=n_cnd> obs2cnd[n_obs];
  int<lower=1, upper=n_tpl> obs2tpl[n_obs];
  int<lower=1, upper=n_stm> obs2stm[n_obs];
  int<lower=1, upper=n_dmg> obs2dmg[n_obs]; 
  int<lower=1, upper=n_tim> obs2tim[n_obs]; 
  int<lower=1, upper=n_l7g> obs2l7g[n_obs];  
  int<lower=1, upper=n_cat> response[n_obs];
  real prior[n_cat-1]; 
}

parameters{
  // vector<lower=0,upper=1>[n_cat-1] guess; 
  vector[n_cat-1] a0;
  vector[n_cat-1] a_lin[n_lin];
  vector[n_cat-1] a_res[n_res];
  vector[n_cat-1] a_tim[n_tim];
  vector[n_cat-1] a_l7g[n_l7g];
  vector[n_cat-1] a_tpl[n_tpl];
  vector[n_cat-1] a_stm[n_stm];
  vector[n_cat-1] a_cnd[n_cnd];
  vector[n_cat-1] a_dmg[n_dmg];
  vector[n_cat-1] a_mus[n_mus];
  real<lower=0> mus_s;
  vector[n_cat-1] a_lin_tim[n_lin,n_tim];
  vector[n_cat-1] a_lin_l7g[n_lin,n_l7g];
  vector[n_cat-1] a_lin_tpl[n_lin,n_tpl];
  vector[n_cat-1] a_lin_stm[n_lin,n_stm];
  vector[n_cat-1] a_lin_cnd[n_lin,n_stm];
  vector[n_cat-1] a_tpl_res[n_tpl,n_res];
}

transformed parameters{

}

model{
  vector[n_cat] mu[n_obs];
  vector[n_cat-1] phi[n_obs];

  // priors
  // guess ~ beta(1,18); 
  a0 ~ normal(logit(prior),1); 
  for (l in 1:n_lin) a_lin[l] ~ normal(0, 1);
  for (r in 1:n_res) a_res[r] ~ normal(0, 1);
  for (a in 1:n_tim) a_tim[a] ~ normal(0, 1);
  for (g in 1:n_l7g) a_l7g[g] ~ normal(0, 1);
  for (t in 1:n_tpl) a_tpl[t] ~ normal(0, 1);
  for (s in 1:n_stm) a_stm[s] ~ normal(0, 1);
  for (c in 1:n_cnd) a_cnd[c] ~ normal(0, 1);
  for (d in 1:n_dmg) a_dmg[d] ~ normal(0, 1);
  for (u in 1:n_mus) a_mus[u] ~ normal(0, mus_s);
  mus_s ~ normal(0,1);
  // interactions
  for(l in 1:3){  //index line
    for(i in 1:3){ //index for cnd tpl stm
      a_lin_cnd[l,i] ~ normal(0,1); 
      a_lin_tpl[l,i] ~ normal(0,1); 
      a_lin_stm[l,i] ~ normal(0,1);       
    }
    for (i in 1:2) {
      a_tpl_res[l,i] ~ normal(0,1); 
      a_lin_tim[l,i] ~ normal(0,1); 
      a_lin_l7g[l,i] ~ normal(0,1); 
    }
  }

  for (n in 1:n_obs){
    for(k in 1:(n_cat-1)){
      phi[n,k] =  inv_logit( a0[k] + 
        a_lin[obs2lin[n],k] + a_res[obs2res[n],k] + a_tim[obs2tim[n],k] + 
        a_l7g[obs2l7g[n],k] + a_tpl[obs2tpl[n],k] + a_stm[obs2stm[n],k] + 
        a_cnd[obs2cnd[n],k] + a_dmg[obs2dmg[n],k] + a_mus[obs2mus[n],k] +
        a_lin_tim[obs2lin[n],obs2tim[n],k] +
        a_lin_l7g[obs2lin[n],obs2l7g[n],k] +
        a_lin_tpl[obs2lin[n],obs2tpl[n],k] +
        a_lin_stm[obs2lin[n],obs2stm[n],k] +
        a_lin_cnd[obs2lin[n],obs2cnd[n],k] +
        a_tpl_res[obs2tpl[n],obs2res[n],k] ); 
    }

    mu[n,1] = phi[n,1]; 
    mu[n,2] = phi[n,2] * (1-phi[n,1]);
    mu[n,3] = phi[n,3] * (1-phi[n,2]) * (1-phi[n,1]);
    mu[n,4] = (1-phi[n,3]) * (1-phi[n,2]) * (1-phi[n,1]);
    response[n] ~ categorical(mu[n]);
  }
}

generated quantities{
  vector[n_cat-1] b0;
  vector[n_cat-1] b_lin[n_lin];
  vector[n_cat-1] b_res[n_res];  
  vector[n_cat-1] b_tim[n_tim];
  vector[n_cat-1] b_l7g[n_l7g];
  vector[n_cat-1] b_cnd[n_cnd];
  vector[n_cat-1] b_stm[n_stm];
  vector[n_cat-1] b_tpl[n_tpl];
  vector[n_cat-1] b_dmg[n_dmg];
  vector[n_cat-1] b_mus[n_mus];
  vector[n_cat-1] m_lin_cnd[n_lin, n_cnd];
  vector[n_cat-1] m_lin_tpl[n_lin, n_tpl];
  vector[n_cat-1] m_lin_tim[n_lin, n_tim];
  vector[n_cat-1] m_lin_stm[n_lin, n_stm];
  vector[n_cat-1] m_lin_l7g[n_lin, n_l7g];
  vector[n_cat-1] m_tpl_res[n_tpl, n_res];
  vector[n_cat-1] b_lin_cnd[n_lin, n_cnd];
  vector[n_cat-1] b_lin_tpl[n_lin, n_tpl];
  vector[n_cat-1] b_lin_stm[n_lin, n_stm];
  vector[n_cat-1] b_lin_tim[n_lin, n_tim];
  vector[n_cat-1] b_lin_l7g[n_lin, n_l7g];
  vector[n_cat-1] b_tpl_res[n_tpl, n_res];
  vector[n_cat-1] m[n_lin, n_res, n_cnd, n_stm, n_tpl, n_tim, n_l7g, n_dmg, n_mus];  //cell means
  vector[n_cat-1] y_pred[n_lin, n_res, n_cnd, n_stm, n_tpl, n_tim, n_l7g];

  b0 = rep_vector(0, n_cat-1); 

  for(l in 1:n_lin){
    for(r in 1:n_res){
      for(c in 1:n_cnd){
        for(s in 1:n_stm){
          for(t in 1:n_tpl){
            for(a in 1:n_tim){
              for(g in 1:n_l7g){
                y_pred[l,r,c,s,t,a,g] = a0 +
                  a_lin[l] + a_res[r] + a_cnd[c] + 
                  a_stm[s] + a_tpl[t] + a_tim[a] + a_l7g[g] + 
                  a_lin_cnd[l,c] + a_lin_tpl[l,t] +
                  a_lin_tim[l,a] + a_lin_stm[l,s] +
                  a_lin_l7g[l,g] + a_tpl_res[t,r];
                for (d in 1:n_dmg){
                  for (u in 1:n_mus){
                    m[l,r,c,s,t,a,g,d,u] = y_pred[l,r,c,s,t,a,g] +
                      a_dmg[d] + a_mus[u];
                    b0 += m[l,r,c,s,t,a,g,d,u];
                  }
                }
              }
            }
          }     
        }
      }
    }
  }
  b0 = b0 / num_elements(m) * (n_cat-1); 
  // main effects
  for(k in 1:(n_cat-1)){
    for(l in 1:n_lin) b_lin[l,k] = mean(to_array_1d(m[l,:,:,:,:,:,:,:,:,k])) - b0[k];
    for(r in 1:n_res) b_res[r,k] = mean(to_array_1d(m[:,r,:,:,:,:,:,:,:,k])) - b0[k];
    for(c in 1:n_cnd) b_cnd[c,k] = mean(to_array_1d(m[:,:,c,:,:,:,:,:,:,k])) - b0[k];
    for(s in 1:n_stm) b_stm[s,k] = mean(to_array_1d(m[:,:,:,s,:,:,:,:,:,k])) - b0[k];
    for(t in 1:n_tpl) b_tpl[t,k] = mean(to_array_1d(m[:,:,:,:,t,:,:,:,:,k])) - b0[k];
    for(a in 1:n_tim) b_tim[a,k] = mean(to_array_1d(m[:,:,:,:,:,a,:,:,:,k])) - b0[k];
    for(g in 1:n_l7g) b_l7g[g,k] = mean(to_array_1d(m[:,:,:,:,:,:,g,:,:,k])) - b0[k];
    for(d in 1:n_dmg) b_dmg[d,k] = mean(to_array_1d(m[:,:,:,:,:,:,:,d,:,k])) - b0[k];
    for(u in 1:n_mus) b_mus[u,k] = mean(to_array_1d(m[:,:,:,:,:,:,:,:,u,k])) - b0[k];
  }

  // interactions
  for (l in 1:n_lin){
    for (c in 1:n_cnd){
      for (k in 1:(n_cat-1)) m_lin_cnd[l,c,k] = mean(to_array_1d(m[l,:,c,:,:,:,:,:,:,k]));
      b_lin_cnd[l,c] = m_lin_cnd[l,c] - (b_lin[l] + b_cnd[c] + b0);
    }
    // interactions: line x stim
    for (s in 1:n_stm){
      for (k in 1:(n_cat-1)) m_lin_stm[l,s,k] = mean(to_array_1d(m[l,:,:,s,:,:,:,:,:,k]));
      b_lin_stm[l,s] = m_lin_stm[l,s] - (b_lin[l] + b_stm[s] + b0);
    }
    // interactions: line x topo
    for (t in 1:n_tpl){
      for (k in 1:(n_cat-1)) m_lin_tpl[l,t,k] = mean(to_array_1d(m[l,:,:,:,t,:,:,:,:,k]));
      b_lin_tpl[l,t] = m_lin_tpl[l,t] - (b_lin[l] + b_tpl[t] + b0);
    }
    // interactions: line x tim
    for (a in 1:n_tim){
      for (k in 1:(n_cat-1)) m_lin_tim[l,a,k] = mean(to_array_1d(m[l,:,:,:,:,a,:,:,:,k]));
      b_lin_tim[l,a] = m_lin_tim[l,a] - (b_lin[l] + b_tim[a] + b0);
    }
    // interactions: line x tim
    for (g in 1:n_l7g){
      for (k in 1:(n_cat-1)) m_lin_l7g[l,g,k] = mean(to_array_1d(m[l,:,:,:,:,:,g,:,:,k]));
      b_lin_l7g[l,g] = m_lin_l7g[l,g] - (b_lin[l] + b_l7g[g] + b0);
    }
  }
  for (t in 1:n_tpl){
    for (r in 1:n_res){
      for (k in 1:(n_cat-1)) m_tpl_res[t,r,k] = mean(to_array_1d(m[:,r,:,:,t,:,:,:,:,k]));
      b_tpl_res[t,r] = m_tpl_res[t,r] - (b_tpl[t] + b_res[r] + b0);
    }
  }

}


