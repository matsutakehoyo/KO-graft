data {
	int<lower=1> n_obs; 
	int<lower=1> n_mus;
	int<lower=1> n_lin;
	int<lower=1> obs2mus[n_obs];       //convert mus to lin id
  int<lower=1> obs2lin[n_obs];        //convert mus to lin id
  int<lower=1, upper=2> obs2sex[n_obs];
	int<lower=0> positive[n_obs];
	int<lower=0> total[n_obs];
  real mean_y; 
  real sd_y;
}

parameters {
  real a0; 
  vector[n_mus] a_mus; 
  real<lower=0> mus_s;
  vector[n_lin] a_lin; 
  vector[2] a_sex; 
  real <lower=0,upper=1> guess; 
  }

transformed parameters{


}

model {
  vector[n_obs] theta;
  // priors
  a0 ~ normal(logit(mean_y),1);
  a_mus ~ normal(0,mus_s);
  mus_s ~ normal(0,1);
  a_lin ~ normal(0,1);
  a_sex ~ normal(0,1); 
  guess ~ beta(1,9);

  for ( i in 1:n_obs ) {
    theta[i] = .5*guess + (1-guess)*inv_logit(
        a0 + a_mus[obs2mus[i]] + a_lin[obs2lin[i]] + a_sex[obs2sex[i]]);
  }  
  
  // likelihood
  for (i in 1:n_obs) {  
    positive[i] ~ binomial(total[i], theta[i]);  
  }
}

generated quantities {
  real b0; 
  vector[n_mus] b_mus; 
  vector[n_lin] b_lin; 
  vector[2] b_sex;
  vector[n_lin] y_pred; // expected value in lins ignoring (assuming mean) other predictors

  {
    real m[n_lin, n_mus, 2]; // cell means
    b0 = 0.0;
    for (l in 1:n_lin){
      y_pred[l] = a0 + a_lin[l]; 
      for (a in 1:n_mus){
        for (s in 1:2){
          m[l,a,s] = a0 + a_mus[a] + a_lin[l] + a_sex[s];
          b0 += m[l,a,s];
        }
      }
    }
    b0 = b0 / num_elements(m); 
    
    // main effects
    for (l in 1:n_lin) b_lin[l] = mean(to_array_1d(m[l,:,:])) - b0; 
    for (a in 1:n_mus) b_mus[a] = mean(to_array_1d(m[:,a,:])) - b0; 
    for (s in 1:2) b_sex[s] = mean(to_array_1d(m[:,:,s])) - b0;
    // -------------------------------------------------------------------------------   
  }

}
