functions{
  //modified Gompertz growth curve
  real gompertz(real A, real mu, real lambda, real time){
    return A * exp(-exp(mu * exp(1) / A * (lambda - time) + 1));
  }
}

data{
  int<lower=1> n_obs;     //number of 3d retina
  int<lower=1> n_line;       //number of cell lines
  int<lower=0> x[n_obs];     //time points
  real<lower=0> y[n_obs];     //rel fuorescence value
  int<lower=1, upper=n_line> obs2line[n_obs];
}

parameters{
  real a0; 
  real m0; 
  real l0; 
  // // line effects
  vector<lower=0>[n_line] a_l;
  vector<lower=0>[n_line] m_l;
  vector<lower=0>[n_line] l_l;
  real<lower=0> rate;
}

transformed parameters{

}

model{
  real mu[n_obs];
  //priors
  a0 ~ normal(11, 2); 
  m0 ~ normal(0.7, 0.5); 
  l0 ~ normal(15, 2); 
  
  rate ~ cauchy(0, 1);

  a_l ~ normal(a0, 2);
  m_l ~ normal(m0, 2);
  l_l ~ normal(l0, 2);
  
  for (i in 1:n_obs){
    mu[i] = gompertz(a_l[obs2line[i]], m_l[obs2line[i]], l_l[obs2line[i]], x[i]);
  }  

 
  for (i in 1:n_obs){
    y[i] ~ gamma(mu[i]*rate, rate);
  }

}


