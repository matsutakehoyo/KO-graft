data {
	int<lower=1> n_obs; 
	int<lower=1> n_aml;
	int<lower=1> n_lin;
	int<lower=1> obs2aml[n_obs];       //convert aml to lin id
  int<lower=1> obs2lin[n_obs];        //convert aml to lin id
	int<lower=0> positive[n_obs];
	int<lower=0> total[n_obs];
  real mean_y; 
  real sd_y;
}

parameters {
  real a0; 
  vector[n_aml] a_aml; 
  real<lower=0> aml_s;
  vector[n_lin] a_lin; 
  real <lower=0,upper=1> guess; 
  }

transformed parameters{


}

model {
  vector[n_obs] theta;
  // priors
  a0 ~ normal(logit(mean_y),1);
  a_aml ~ normal(0,aml_s);
  aml_s ~ normal(0,1);
  a_lin ~ normal(0,1);
  // a_lin_s ~ normal(0,1);
  guess ~ beta(1,18);

  for ( i in 1:n_obs ) {
    theta[i] = .5*guess + (1-guess)*inv_logit(
    // theta[i] = inv_logit(
        a0 + a_aml[obs2aml[i]] + a_lin[obs2lin[i]]);
  }  
  
  // likelihood
  for (i in 1:n_obs) {  
    positive[i] ~ binomial(total[i], theta[i]);  
  }
}

generated quantities {
  real b0; 
  vector[n_aml] b_aml; 
  vector[n_lin] b_lin; 
  vector[n_lin] y_pred; // expected value in lins ignoring (assuming mean) other predictors

  {
    real m[n_lin, n_aml]; // cell means
    b0 = 0.0;
    for (l in 1:n_lin){
      y_pred[l] = a0 + a_lin[l]; 
      for (a in 1:n_aml){
          m[l,a] = a0 + a_aml[a] + a_lin[l];
          b0 += m[l,a];
      }
    }
    b0 = b0 / num_elements(m); 
    
    // main effects
    for (l in 1:n_lin) b_lin[l] = mean(to_array_1d(m[l,:])) - b0; 
    for (a in 1:n_aml) b_aml[a] = mean(to_array_1d(m[:,a])) - b0; 
    // -------------------------------------------------------------------------------   
  }

}
