data {
	int<lower=1> n_obs;
	int<lower=1> n_lin;
	int<lower=1> n_tim;
	int<lower=1> n_l7g;
	int<lower=1> n_tpl;
	int<lower=1> n_stm;
	int<lower=1> n_cnd;
	int<lower=1> n_dmg;
	int<lower=1> n_mus;
	int<lower=1, upper=n_lin> obs2lin[n_obs];
	int<lower=1, upper=n_tim> obs2tim[n_obs];
	int<lower=1, upper=n_l7g> obs2l7g[n_obs];
	int<lower=1, upper=n_tpl> obs2tpl[n_obs];
	int<lower=1, upper=n_stm> obs2stm[n_obs];
	int<lower=1, upper=n_cnd> obs2cnd[n_obs];
	int<lower=1, upper=n_dmg> obs2dmg[n_obs];
	int<lower=1, upper=n_mus> obs2mus[n_obs];	
	// int<lower=1, upper=2> obs2bwv[n_obs];
	real obs2bwv[n_obs];
	real obs2spt[n_obs];
	int<lower=0, upper=1> response[n_obs];
	real<lower=0, upper=1> p_non_zero; 
}


parameters {
	// overall mean
	real a0; 
	// main effects
	vector[n_mus] a_mus;
	real<lower=0> mus_s; 
	vector[n_lin] a_lin;
	vector[n_stm] a_stm;
	vector[n_cnd] a_cnd;
	vector[n_tpl] a_tpl;
	vector[n_dmg] a_dmg;
	vector[n_tim] a_tim;
	vector[n_l7g] a_l7g;
	// covariates
	real bwv;
	real spt;
	// interactions
	matrix[n_lin,n_cnd] a_lin_cnd;
	matrix[n_lin,n_tpl] a_lin_tpl;
	matrix[n_lin,n_stm] a_lin_stm;
	matrix[n_lin,n_tim] a_lin_tim;
	matrix[n_lin,n_l7g] a_lin_l7g;
	//guessing parameter for Binomial robustness
	real<lower=0,upper=1> guess_theta;		
}

transformed parameters {

}

model {
	real theta[n_obs];
	// for robust estimate
	guess_theta ~ beta(1,9);
	// overall mean
	a0 ~ normal(logit(p_non_zero), 1); 
	// main effects
	a_mus ~ normal(0,mus_s);	
	mus_s ~ gamma(1.64,0.32); //mode=2, sd=4
	a_lin ~ student_t(3,0,1);
	a_cnd ~ student_t(3,0,1);
	a_stm ~ student_t(3,0,1);
	a_tpl ~ student_t(3,0,1);
	a_dmg ~ student_t(3,0,1);
	a_tim ~ student_t(3,0,1);
	a_l7g ~ student_t(3,0,1);
	// covariates
	bwv ~ student_t(3,0,1);
	spt ~ student_t(3,0,1);
	// interactions
	for(i1 in 1:3){	//index line
		for(i2 in 1:3){	//index cond and topo
			a_lin_cnd[i1,i2] ~ student_t(3,0,1);
			a_lin_tpl[i1,i2] ~ student_t(3,0,1);
		}
	}
	for(i1 in 1:3){	//index line
		for(i2 in 1:2){	//index time and stim
			a_lin_stm[i1,i2] ~ student_t(3,0,1);
			a_lin_tim[i1,i2] ~ student_t(3,0,1);
			a_lin_l7g[i1,i2] ~ student_t(3,0,1);
		}
	}

	for ( i in 1:n_obs ) {
		theta[i] =.5*guess_theta + (1-guess_theta)*inv_logit(
			a0 + 
			a_mus[obs2mus[i]] + 
			a_lin[obs2lin[i]] + 
			a_stm[obs2stm[i]] +
			a_cnd[obs2cnd[i]] +
			a_tpl[obs2tpl[i]] +
			a_dmg[obs2dmg[i]] +
			a_tim[obs2tim[i]] +
			a_l7g[obs2l7g[i]] +
			a_lin_cnd[obs2lin[i],obs2cnd[i]] +
			a_lin_stm[obs2lin[i],obs2stm[i]] +
			a_lin_tpl[obs2lin[i],obs2tpl[i]] +
			a_lin_tim[obs2lin[i],obs2tim[i]] +
			a_lin_l7g[obs2lin[i],obs2l7g[i]] +
			spt * obs2spt[i] +
			bwv * obs2bwv[i] 
		);
	}
	// likelihood: Bernoulli
	for (i in 1:n_obs) {
		response[i] ~ bernoulli(theta[i]);
	}
}

generated quantities {
	real t0;
	vector[n_lin] t_lin;
	vector[n_mus] t_mus; 
	vector[n_tim] t_tim;
	vector[n_stm] t_stm; 
	vector[n_cnd] t_cnd; 
	vector[n_tpl] t_tpl; 
	vector[n_l7g] t_l7g; 
	vector[n_dmg] t_dmg; 
	matrix[n_lin, n_cnd] m_lin_cnd;
	matrix[n_lin, n_stm] m_lin_stm;
	matrix[n_lin, n_tpl] m_lin_tpl;
	matrix[n_lin, n_tim] m_lin_tim;
	matrix[n_lin, n_l7g] m_lin_l7g;
	matrix[n_lin, n_cnd] t_lin_cnd;
	matrix[n_lin, n_stm] t_lin_stm;
	matrix[n_lin, n_tpl] t_lin_tpl;
	matrix[n_lin, n_tim] t_lin_tim;
	matrix[n_lin, n_l7g] t_lin_l7g;
	real m[n_lin, n_tim, n_l7g, n_tpl, n_stm, n_cnd, n_dmg, n_mus]; // cell means for b/l
	// estimated mean 
	real y_pred[n_lin, n_tim, n_l7g, n_tpl, n_stm, n_cnd];

	// convert unconstrained predicotrs to constrained (sum-to-zero)
	// --------------------------------------------------------------------------------
	// overall mean
	t0 = 0.0; 
	// cell means
	for(l in 1:n_lin){
		for(a in 1:n_tim){
			for(g in 1:n_l7g){
				for(t in 1:n_tpl){
					for(s in 1:n_stm){
						for(c in 1:n_cnd){
							y_pred[l,a,g,t,s,c] = a0 +
								a_lin[l] + a_tim[a] + a_l7g[g] + a_tpl[t] + a_stm[s] + a_cnd[c] +
								a_lin_l7g[l,g] + a_lin_stm[l,s] + a_lin_tpl[l,t] + a_lin_cnd[l,c];
							for(d in 1:n_dmg){
								for(u in 1:n_mus){
									m[l,a,g,t,s,c,d,u] = a0 +
										a_lin[l] + a_tim[a] + a_l7g[g] + a_tpl[t] + a_stm[s] + a_cnd[c] +
										a_lin_l7g[l,g] + a_lin_stm[l,s] + a_lin_tpl[l,t] + a_lin_cnd[l,c] + 
										a_dmg[d] + a_mus[u];
									t0 += m[l,a,g,t,s,c,d,u];
								}
							}
						}
					}
				}
			}
		}
	}

	t0 = t0 / num_elements(m);	
	// main effects
	for(l in 1:n_lin)	t_lin[l] = mean(to_array_1d(m[l,:,:,:,:,:,:,:])) - t0;
	for(a in 1:n_tim)	t_tim[a] = mean(to_array_1d(m[:,a,:,:,:,:,:,:])) - t0; 
	for(g in 1:n_l7g) t_l7g[g] = mean(to_array_1d(m[:,:,g,:,:,:,:,:])) - t0; 
	for(t in 1:n_tpl)	t_tpl[t] = mean(to_array_1d(m[:,:,:,t,:,:,:,:])) - t0; 
	for(s in 1:n_stm)	t_stm[s] = mean(to_array_1d(m[:,:,:,:,s,:,:,:])) - t0; 
	for(c in 1:n_cnd) t_cnd[c] = mean(to_array_1d(m[:,:,:,:,:,c,:,:])) - t0; 
	for(d in 1:n_dmg)	t_dmg[d] = mean(to_array_1d(m[:,:,:,:,:,:,d,:])) - t0; 
	for(u in 1:n_mus)	t_mus[u] = mean(to_array_1d(m[:,:,:,:,:,:,:,u])) - t0; 		
	
	// interactions
	for (l in 1:n_lin){
		for (a in 1:n_tim){
			m_lin_tim[l,a] = mean(to_array_1d(m[l,a,:,:,:,:,:,:]));
			t_lin_tim[l,a] = m_lin_tim[l,a] - (t_lin[l] + t_tim[a] + t0);
		}
		for (g in 1:n_l7g){
			m_lin_l7g[l,g] = mean(to_array_1d(m[l,:,g,:,:,:,:,:]));
			t_lin_l7g[l,g] = m_lin_l7g[l,g] - (t_lin[l] + t_l7g[g] + t0);
		}
		for (t in 1:n_tpl){
			m_lin_tpl[l,t] = mean(to_array_1d(m[l,:,:,t,:,:,:,:]));
			t_lin_tpl[l,t] = m_lin_tpl[l,t] - (t_lin[l] + t_tpl[t] + t0);
		}

		for (s in 1:n_stm){
			m_lin_stm[l,s] = mean(to_array_1d(m[l,:,:,:,s,:,:,:]));
			t_lin_stm[l,s] = m_lin_stm[l,s] - (t_lin[l] + t_stm[s] + t0);
		}
		for (c in 1:n_cnd){
			m_lin_cnd[l,c] = mean(to_array_1d(m[l,:,:,:,:,c,:,:]));
			t_lin_cnd[l,c] = m_lin_cnd[l,c] - (t_lin[l] + t_cnd[c] + t0);
		}

	}
}

