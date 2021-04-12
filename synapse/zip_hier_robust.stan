// Zero-Inflated Hierarchical Poisson model
data {
	int<lower=1> n_obs;
	int<lower=1> n_mus;
	int<lower=1> n_lin;
	int<lower=1> n_sex;
	int<lower=1, upper=n_lin> obs2lin[n_obs];
	int<lower=1, upper=n_mus> obs2mus[n_obs];
	int<lower=1, upper=n_sex> obs2sex[n_obs];
	int<lower=0>  y[n_obs];
	real<lower=0>  y_mean; // mean of non-zero counts
	real<lower=0>  p_non_zero; // mean of non-zero counts

}

parameters {
	real a0; 
	real b0; 
	vector[n_lin] a_lin;
	vector[n_lin] b_lin;
	vector[n_mus] a_mus;
	vector[n_mus] b_mus;
	vector[n_sex] a_sex;
	vector[n_sex] b_sex;
	real<lower=0> a_mus_s;
	real<lower=0> b_mus_s;
	matrix[n_lin,n_sex] a_lin_sex;
	matrix[n_lin,n_sex] b_lin_sex;
	real<lower=0,upper=1> guess; 		//guessing parameter for Bernoulli robustness
}

transformed parameters { 
}

model {
  vector[n_obs] theta;
  vector[n_obs] lambda;


	guess ~ beta(1,9);
	a0 ~ normal(logit(p_non_zero), 1);
	b0 ~ normal(log(y_mean), 1); 
	a_mus_s ~ normal(0,1);
	b_mus_s ~ normal(0,1);
	a_mus ~ normal(0, a_mus_s);
	b_mus ~ normal(0, b_mus_s);	
	a_lin ~ normal(0,1);
	b_lin ~ normal(0,1);
	a_sex ~ normal(0, 1);
	b_sex ~ normal(0, 1);

	// interactions
	for(l in 1:n_lin){	//index line
		for(s in 1:n_sex){	//index sex
			a_lin_sex[l,s] ~ normal(0,1); 
			b_lin_sex[l,s] ~ normal(0,1); 
		}
	}

  for ( i in 1:n_obs ) {
    theta[i] = .5*guess + (1-guess)*inv_logit(
			a0 + 
			a_lin[obs2lin[i]] + 
			a_mus[obs2mus[i]] + 
			a_sex[obs2sex[i]] +
			a_lin_sex[obs2lin[i],obs2sex[i]]
		);
    
    lambda[i] = exp(
			b0 + 
			b_lin[obs2lin[i]] +
			b_mus[obs2mus[i]] + 
			b_sex[obs2sex[i]] +
			b_lin_sex[obs2lin[i],obs2sex[i]]
		);
  }

	// likelihood: Zero_Inflated Poisson
	for (i in 1:n_obs) {
		if (y[i] == 0) {
			target += log_sum_exp(bernoulli_lpmf(0 | theta[i]),
			                      bernoulli_lpmf(1 | theta[i]) + poisson_lpmf(y[i] | lambda[i]));
		} else {
			target += bernoulli_lpmf(1 | theta[i]) + poisson_lpmf(y[i] | lambda[i]);
		}
	}
}

generated quantities {
	real t0; 
	real l0; 
	vector[n_lin] t_lin;
	vector[n_lin] l_lin;
	vector[n_sex] t_sex;
	vector[n_sex] l_sex;	
	vector[n_mus] t_mus;
	vector[n_mus] l_mus;
	matrix[n_lin, n_sex] t_lin_sex;
	matrix[n_lin, n_sex] l_lin_sex;
	matrix[n_lin, n_sex] m_lin_sex;
	matrix[n_lin, n_sex] n_lin_sex;
	real m[n_lin, n_sex, n_mus];
	real n[n_lin, n_sex, n_mus];
	real t_pred[n_lin, n_sex];
	real l_pred[n_lin, n_sex];

	t0 = 0.0;
	l0 = 0.0;
	for (l in 1:n_lin){
		for (s in 1:2){
			t_pred[l,s] = a0 + a_lin[l] + a_sex[s] + a_lin_sex[l,s];
			l_pred[l,s] = b0 + b_lin[l] + b_sex[s] + b_lin_sex[l,s];
			for (u in 1:n_mus){
					m[l,s,u] = a0 + a_lin[l] + a_mus[u] + a_sex[s] + a_lin_sex[l,s];
					n[l,s,u] = b0 + b_lin[l] + b_mus[u] + b_sex[s] + b_lin_sex[l,s]; 
					t0 += m[l,s,u];
					l0 += n[l,s,u];
			}
		}
	}
	t0 = t0 / num_elements(m);
	l0 = l0 / num_elements(n);

	for (l in 1:n_lin) {
		t_lin[l] = mean(to_array_1d(m[l,:,:])) - t0; 
		l_lin[l] = mean(to_array_1d(n[l,:,:])) - l0; 
	} 
	for (s in 1:n_sex) {
		t_sex[s] = mean(to_array_1d(m[:,s,:])) - t0; 
		l_sex[s] = mean(to_array_1d(n[:,s,:])) - l0;
	} 
	for (u in 1:n_mus) {
		t_mus[u] = mean(to_array_1d(m[:,:,u])) - t0; 	
		l_mus[u] = mean(to_array_1d(n[:,:,u])) - l0; 	
	}

	for (l in 1:n_lin){
		for (s in 1:n_sex){
			m_lin_sex[l,s] = mean(to_array_1d(m[l,s,:]));
			n_lin_sex[l,s] = mean(to_array_1d(n[l,s,:]));
			t_lin_sex[l,s] = m_lin_sex[l,s] - (t_lin[l] + t_sex[s] + t0);
			l_lin_sex[l,s] = n_lin_sex[l,s] - (l_lin[l] + l_sex[s] + l0);

		}
	}

}

