data{
  int<lower=1> n_obs;
  int<lower=1> n_cat; 
  int<lower=1> n_mus; 
  int<lower=1> n_lin;
  int<lower=1> n_tpl;
  int<lower=1> n_tim;
  int<lower=1> n_stm;
  int<lower=1> n_dmg;
  int<lower=1> n_l7g;
  int<lower=1, upper=n_lin> obs2lin[n_obs];
  int<lower=1, upper=n_mus> obs2mus[n_obs];
  int<lower=1, upper=n_tpl> obs2tpl[n_obs];
  int<lower=1, upper=n_stm> obs2stm[n_obs];
  int<lower=1, upper=n_dmg> obs2dmg[n_obs]; 
  int<lower=1, upper=n_tim> obs2tim[n_obs]; 
  int<lower=1, upper=n_l7g> obs2l7g[n_obs];  
  real obs2spt[n_obs];
  int<lower=1, upper=n_cat> response[n_obs];
  real prior[n_cat-1]; 
}

parameters{
  vector<lower=0,upper=1>[n_cat-1] guess; 
  vector[n_cat-1] a0;
  vector[n_cat-1] a_lin[n_lin];
  vector[n_cat-1] a_tim[n_tim];
  vector[n_cat-1] a_l7g[n_l7g];
  vector[n_cat-1] a_tpl[n_tpl];
  vector[n_cat-1] a_stm[n_stm];
  vector[n_cat-1] a_dmg[n_dmg];
  vector[n_cat-1] a_mus[n_mus];
  real<lower=0> mus_s;
  vector[n_cat-1] a_lin_tim[n_lin,n_tim];
  vector[n_cat-1] a_lin_l7g[n_lin,n_l7g];
  vector[n_cat-1] a_lin_tpl[n_lin,n_tpl];
  vector[n_cat-1] a_lin_stm[n_lin,n_stm];
  vector[n_cat-1] spt;
}

transformed parameters{

}

model{
  vector[n_cat] mu[n_obs];
  vector[n_cat-1] phi[n_obs];

  // priors
  guess ~ beta(1,9);
  
  a0 ~ normal(logit(prior),1); 
  // main effects
  for (l in 1:n_lin) a_lin[l] ~ student_t(3,0,1);
  for (a in 1:n_tim) a_tim[a] ~ student_t(3,0,1);
  for (g in 1:n_l7g) a_l7g[g] ~ student_t(3,0,1);
  for (t in 1:n_tpl) a_tpl[t] ~ student_t(3,0,1);
  for (s in 1:n_stm) a_stm[s] ~ student_t(3,0,1);
  for (d in 1:n_dmg) a_dmg[d] ~ student_t(3,0,1);
  for (u in 1:n_mus) a_mus[u] ~ normal(0, mus_s);
  mus_s ~ gamma(1.64,0.32); //mode=2, sd=4
  spt ~ student_t(3,0,1);

  // interactions
  for(l in 1:3){  //index line
    for(i in 1:3){ //index for cnd tpl stm
      a_lin_tpl[l,i] ~ student_t(3,0,1);
      a_lin_stm[l,i] ~ student_t(3,0,1);
    }
    for (i in 1:n_tim) {
      a_lin_tim[l,i] ~ student_t(3,0,1);
      a_lin_l7g[l,i] ~ student_t(3,0,1);
    }
  }

  for (n in 1:n_obs){
    for(k in 1:(n_cat-1)){
      phi[n,k] = .5*guess[k] + (1-guess[k]) * inv_logit( a0[k] + a_lin[obs2lin[n],k] + 
        a_tim[obs2tim[n],k] + a_l7g[obs2l7g[n],k] +
        a_tpl[obs2tpl[n],k] + a_stm[obs2stm[n],k] + 
        a_dmg[obs2dmg[n],k] + a_mus[obs2mus[n],k] + 
        a_lin_tim[obs2lin[n],obs2tim[n],k] +
        a_lin_l7g[obs2lin[n],obs2l7g[n],k] +
        a_lin_tpl[obs2lin[n],obs2tpl[n],k] +
        a_lin_stm[obs2lin[n],obs2stm[n],k] +
        spt[k] * obs2spt[n] ); 
    }

    mu[n,1] = phi[n,1]; // silent
    mu[n,2] = phi[n,2] * (1-phi[n,1]); //onxoff
    mu[n,3] = phi[n,3] * (1-phi[n,2]) * (1-phi[n,1]); //on 
    mu[n,4] = (1-phi[n,3]) * (1-phi[n,2]) * (1-phi[n,1]); //adapted_on
    response[n] ~ categorical(mu[n]);
  }
}

generated quantities{
  vector[n_cat-1] b0;
  vector[n_cat-1] b_lin[n_lin];
  vector[n_cat-1] b_tim[n_tim];
  vector[n_cat-1] b_l7g[n_l7g];
  vector[n_cat-1] b_stm[n_stm];
  vector[n_cat-1] b_tpl[n_tpl];
  vector[n_cat-1] b_dmg[n_dmg];
  vector[n_cat-1] b_mus[n_mus];
  vector[n_cat-1] m_lin_tpl[n_lin, n_tpl];
  vector[n_cat-1] m_lin_tim[n_lin, n_tim];
  vector[n_cat-1] m_lin_stm[n_lin, n_stm];
  vector[n_cat-1] m_lin_l7g[n_lin, n_l7g];
  vector[n_cat-1] b_lin_tpl[n_lin, n_tpl];
  vector[n_cat-1] b_lin_stm[n_lin, n_stm];
  vector[n_cat-1] b_lin_tim[n_lin, n_tim];
  vector[n_cat-1] b_lin_l7g[n_lin, n_l7g];
  vector[n_cat-1] m[n_lin, n_stm, n_tpl, n_tim, n_l7g, n_dmg, n_mus];  //cell means
  vector[n_cat-1] y_pred[n_lin, n_stm, n_tpl, n_tim, n_l7g];

  b0 = rep_vector(0, n_cat-1); 

  for(l in 1:n_lin){
    for(s in 1:n_stm){
      for(t in 1:n_tpl){
        for(a in 1:n_tim){
          for(g in 1:n_l7g){
            y_pred[l,s,t,a,g] = a0 + a_lin[l] + 
              a_stm[s] + a_tpl[t] + 
              a_tim[a] + a_l7g[g] + 
              a_lin_tpl[l,t] + a_lin_tim[l,a] + 
              a_lin_stm[l,s] + a_lin_l7g[l,g];
            for (d in 1:n_dmg){
              for (u in 1:n_mus){
                m[l,s,t,a,g,d,u] = a0 + a_lin[l] + 
                  a_tim[a] + a_l7g[g] + a_stm[s] + 
                  a_tpl[t] + a_dmg[d] + a_mus[u] +
                  a_lin_tpl[l,t] + a_lin_tim[l,a] + 
                  a_lin_stm[l,s] + a_lin_l7g[l,g];
                b0 += m[l,s,t,a,g,d,u];
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
    for(l in 1:n_lin) b_lin[l,k] = mean(to_array_1d(m[l,:,:,:,:,:,:,k])) - b0[k];
    for(s in 1:n_stm) b_stm[s,k] = mean(to_array_1d(m[:,s,:,:,:,:,:,k])) - b0[k];
    for(t in 1:n_tpl) b_tpl[t,k] = mean(to_array_1d(m[:,:,t,:,:,:,:,k])) - b0[k];
    for(a in 1:n_tim) b_tim[a,k] = mean(to_array_1d(m[:,:,:,a,:,:,:,k])) - b0[k];
    for(g in 1:n_l7g) b_l7g[g,k] = mean(to_array_1d(m[:,:,:,:,g,:,:,k])) - b0[k];
    for(d in 1:n_dmg) b_dmg[d,k] = mean(to_array_1d(m[:,:,:,:,:,d,:,k])) - b0[k];
    for(u in 1:n_mus) b_mus[u,k] = mean(to_array_1d(m[:,:,:,:,:,:,u,k])) - b0[k];
  }

  for (l in 1:n_lin){
    // interactions: line x tim
    for (a in 1:n_tim){
      for (k in 1:(n_cat-1)) m_lin_tim[l,a,k] = mean(to_array_1d(m[l,:,:,a,:,:,:,k]));
      b_lin_tim[l,a] = m_lin_tim[l,a] - (b_lin[l] + b_tim[a] + b0);
    }
    // interactions: line x l7g
    for (g in 1:n_l7g){
      for (k in 1:(n_cat-1)) m_lin_l7g[l,g,k] = mean(to_array_1d(m[l,:,:,:,g,:,:,k]));
      b_lin_l7g[l,g] = m_lin_l7g[l,g] - (b_lin[l] + b_l7g[g] + b0);
    }
    // interactions: line x topo
    for (t in 1:n_tpl){
      for (k in 1:(n_cat-1)) m_lin_tpl[l,t,k] = mean(to_array_1d(m[l,:,t,:,:,:,:,k]));
      b_lin_tpl[l,t] = m_lin_tpl[l,t] - (b_lin[l] + b_tpl[t] + b0);
    }
    // interactions: line x stim
    for (s in 1:n_stm){
      for (k in 1:(n_cat-1)) m_lin_stm[l,s,k] = mean(to_array_1d(m[l,s,:,:,:,:,:,k]));
      b_lin_stm[l,s] = m_lin_stm[l,s] - (b_lin[l] + b_stm[s] + b0);
    }
  }

}


