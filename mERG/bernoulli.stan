// Hierarchical normal crossed model
data {
	int<lower=1> n_obs;
	int<lower=1> n_mus;
	int<lower=1> n_lin;
	int<lower=1> n_cnd;
	int<lower=1> n_stm;
	int<lower=1> n_tpl;
	int<lower=1> n_tim;
	int<lower=1> n_dmg;
	int<lower=1> n_l7g;
	int<lower=1, upper=n_lin> obs2lin[n_obs];
	int<lower=1, upper=n_mus> obs2mus[n_obs];
	int<lower=1, upper=n_stm> obs2stm[n_obs];
	int<lower=1, upper=n_cnd> obs2cnd[n_obs];
	int<lower=1, upper=n_tpl> obs2tpl[n_obs];
	int<lower=1, upper=n_dmg> obs2dmg[n_obs];
	int<lower=1, upper=n_l7g> obs2l7g[n_obs];
	int<lower=1, upper=n_tim> obs2tim[n_obs];
	int<lower=0, upper=1> y[n_obs];
	real<lower=0> p_non_zero; 
}

parameters {
	real<lower=0,upper=1> guess;
	// main effects
	real a0; 
	vector[n_lin] a_lin;
	vector[n_mus] a_mus;
	real<lower=0> mus_s; 
	vector[n_stm] a_stm;
	vector[n_tim] a_tim;
	vector[n_cnd] a_cnd;
	vector[n_tpl] a_tpl;
	vector[n_dmg] a_dmg;
	vector[n_l7g] a_l7g;
	// interactions
	matrix[n_lin, n_cnd] a_lin_cnd;
	matrix[n_lin, n_tpl] a_lin_tpl;
	matrix[n_lin, n_stm] a_lin_stm;
	matrix[n_lin, n_tim] a_lin_tim;
	matrix[n_lin, n_l7g] a_lin_l7g;
}

transformed parameters {
	vector[n_obs] theta;

	for (i in 1:n_obs) {
		theta[i] = .5*guess + (1-guess) *  inv_logit( 
			a0 + 
			a_lin[obs2lin[i]] + 
			a_mus[obs2mus[i]] +
			a_stm[obs2stm[i]] +
			a_cnd[obs2cnd[i]] +
			a_tpl[obs2tpl[i]] +
			a_dmg[obs2dmg[i]] +
			a_tim[obs2tim[i]] +
			a_l7g[obs2l7g[i]] +
			a_lin_cnd[obs2lin[i], obs2cnd[i]] +
			a_lin_tpl[obs2lin[i], obs2tpl[i]] +
			a_lin_stm[obs2lin[i], obs2stm[i]] +
			a_lin_tim[obs2lin[i], obs2tim[i]] +
			a_lin_l7g[obs2lin[i], obs2l7g[i]]
			);
	}
}

model {
	a0 ~ normal(logit(p_non_zero), 1);
	a_lin ~ student_t(3,0,1);
	a_mus ~ normal(0,mus_s);
	mus_s ~ gamma(1.64,0.32); //mode=2, sd=4
	a_stm ~ student_t(3,0,1);
	a_cnd ~ student_t(3,0,1);
	a_tpl ~ student_t(3,0,1);
	a_dmg ~ student_t(3,0,1);
	a_tim ~ student_t(3,0,1);
	a_l7g ~ student_t(3,0,1);
 	// interactions
 	for(i1 in 1:n_lin) { //index of lin
 		for(i2 in 1:3) { //index of cnd == tpl ==3
			a_lin_cnd[i1, i2] ~ student_t(3,0,1); 
			a_lin_tpl[i1, i2] ~ student_t(3,0,1);  			
 		}
 	}
 	for(i1 in 1:n_lin){ //# of lin
 		for(i2 in 1:2){ //# of tim==stm==2
			a_lin_stm[i1, i2] ~ student_t(3,0,1); 
			a_lin_tim[i1, i2] ~ student_t(3,0,1); 
			a_lin_l7g[i1, i2] ~ student_t(3,0,1); 
 		}
 	}

	// for robust estmation
	guess ~ beta(1,9);
	y ~ bernoulli(theta);	
}

generated quantities {
  // for loo 
  // vector[n_obs] log_lik;
	// change a to deflaction from overall mean: t
	real t0; 
	// main effects
	vector[n_lin] t_lin;
	vector[n_mus] t_mus; 
	vector[n_stm] t_stm;
	vector[n_cnd] t_cnd;
	vector[n_tpl] t_tpl; 
	vector[n_dmg] t_dmg; 
	vector[n_l7g] t_l7g; 
	vector[n_tim] t_tim; 
	// interactions
	matrix[n_lin, n_cnd] m_lin_cnd;
	matrix[n_lin, n_tpl] m_lin_tpl;
	matrix[n_lin, n_stm] m_lin_stm;
	matrix[n_lin, n_tim] m_lin_tim;
	matrix[n_lin, n_l7g] m_lin_l7g;
	matrix[n_lin, n_cnd] t_lin_cnd;
	matrix[n_lin, n_tpl] t_lin_tpl;
	matrix[n_lin, n_stm] t_lin_stm;
	matrix[n_lin, n_tim] t_lin_tim;
	matrix[n_lin, n_l7g] t_lin_l7g;
	// cell means
	real m[n_lin, n_mus, n_cnd, n_tpl, n_dmg, n_stm, n_tim, n_l7g]; // cell means for m
	// estimated mean 
	real y_pred[n_lin, n_cnd, n_tpl, n_stm, n_tim, n_l7g];
	// to calculate difference 
	vector[3] diff_lin_t; 
	real diff_stm_t; 
	vector[3] diff_cnd_t; 
	vector[3] diff_tpl_t; 
	real diff_dmg_t; 
	real diff_tim_t; 
	real diff_l7g_t;
	// --------------------------------------------------------------------------------
	// for loo
	// for (n in 1:n_obs) {
	// 	log_lik[n] = bernoulli_logit_lpmf( y[n] | theta[n] );			
	// }
	// --------------------------------------------------------------------------------
	// convert unconstrained predicotrs to constrained (sum-to-zero)
	// calculate sum to zero parameters
	// overall mean
	t0 = 0.0; 
	// cell means
	for(l in 1:n_lin){
		for(r in 1:n_mus){
			for(c in 1:n_cnd){
				for(t in 1:n_tpl){
					for(d in 1:n_dmg){
						for(s in 1:n_stm){
							for(a in 1:n_tim){
								for(j in 1:n_l7g){
									m[l,r,c,t,d,s,a,j] = a0 +
										a_lin[l] +
										a_mus[r] +
										a_cnd[c] +
										a_tpl[t] +
										a_dmg[d] + 
										a_stm[s] +
										a_tim[a] +
										a_l7g[j] +
										a_lin_cnd[l,c] +
										a_lin_tpl[l,t] +
										a_lin_stm[l,s] +
										a_lin_tim[l,a] +
										a_lin_l7g[l,j];
									t0 += m[l,r,c,t,d,s,a,j];
								}
							}
						}
					}
				}
			}
		}
	}

	for(l in 1:n_lin){
		for(c in 1:n_cnd){
			for(t in 1:n_tpl){
				for(s in 1:n_stm){
					for(a in 1:n_tim){
						for(j in 1:n_l7g){
							y_pred[l,c,t,s,a,j] = a0 +
								a_lin[l] +
								a_cnd[c] +
								a_tpl[t] +
								a_stm[s] +
								a_tim[a] +
								a_l7g[j] +
								a_lin_cnd[l,c] +
								a_lin_tpl[l,t] +
								a_lin_stm[l,s] +
								a_lin_tim[l,a] +
								a_lin_l7g[l,j];
						}
					}
				}
			}
		}
	}
	t0 = t0 / num_elements(m);
	// main effects
	for(l in 1:n_lin) t_lin[l] = mean(to_array_1d(m[l,:,:,:,:,:,:,:])) - t0; 
	for(r in 1:n_mus)	t_mus[r] = mean(to_array_1d(m[:,r,:,:,:,:,:,:])) - t0; 
	for(c in 1:n_cnd) t_cnd[c] = mean(to_array_1d(m[:,:,c,:,:,:,:,:])) - t0; 
	for(t in 1:n_tpl)	t_tpl[t] = mean(to_array_1d(m[:,:,:,t,:,:,:,:])) - t0; 
	for(d in 1:n_dmg)	t_dmg[d] = mean(to_array_1d(m[:,:,:,:,d,:,:,:])) - t0; 
	for(s in 1:n_stm) t_stm[s] = mean(to_array_1d(m[:,:,:,:,:,s,:,:])) - t0; 
	for(a in 1:n_tim) t_tim[a] = mean(to_array_1d(m[:,:,:,:,:,:,a,:])) - t0; 
	for(j in 1:n_l7g)	t_l7g[j] = mean(to_array_1d(m[:,:,:,:,:,:,:,j])) - t0; 
	
	// interactions
	for (l in 1:n_lin){
		for (c in 1:n_cnd){
			m_lin_cnd[l,c] = mean(to_array_1d(m[l,:,c,:,:,:,:,:]));
			t_lin_cnd[l,c] = m_lin_cnd[l,c] - (t_lin[l] + t_cnd[c] + t0);			
		}
	}
	for (l in 1:n_lin){
		for (t in 1:n_tpl){
			m_lin_tpl[l,t] = mean(to_array_1d(m[l,:,:,t,:,:,:,:]));
			t_lin_tpl[l,t] = m_lin_tpl[l,t] - (t_lin[l] + t_tpl[t] + t0);
		}
	}
	for (l in 1:n_lin){
		for (s in 1:n_stm){
			m_lin_stm[l,s] = mean(to_array_1d(m[l,:,:,:,:,s,:,:]));
			t_lin_stm[l,s] = m_lin_stm[l,s] - (t_lin[l] + t_stm[s] + t0);
		}
	}
	for (l in 1:n_lin){
		for (a in 1:n_tim){
			m_lin_tim[l,a] = mean(to_array_1d(m[l,:,:,:,:,:,a,:]));
			t_lin_tim[l,a] = m_lin_tim[l,a] - (t_lin[l] + t_tim[a] + t0);			
		}
	}
	for (l in 1:n_lin){
		for (j in 1:n_l7g){
			m_lin_l7g[l,j] = mean(to_array_1d(m[l,:,:,:,:,:,:,j]));
			t_lin_l7g[l,j] = m_lin_l7g[l,j] - (t_lin[l] + t_l7g[j] + t0);			
		}
	}


	// --------------------------------------------------------------------------------
	// differences
	diff_lin_t[1] = t_lin[1] - t_lin[2];
	diff_lin_t[2] = t_lin[1] - t_lin[3];
	diff_lin_t[3] = t_lin[2] - t_lin[3];
	diff_stm_t = t_stm[1] - t_stm[2];
	diff_tim_t = t_tim[1] - t_tim[2];
	diff_cnd_t[1] = t_cnd[1] - t_cnd[2];
	diff_cnd_t[2] = t_cnd[1] - t_cnd[3];
	diff_cnd_t[3] = t_cnd[2] - t_cnd[3];
	diff_tpl_t[1] = t_tpl[1] - t_tpl[2];
	diff_tpl_t[2] = t_tpl[1] - t_tpl[3];
	diff_tpl_t[3] = t_tpl[2] - t_tpl[3];
	diff_dmg_t = t_dmg[1] - t_dmg[2];
	diff_l7g_t = t_l7g[1] - t_l7g[2];
}




